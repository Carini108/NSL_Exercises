#include "mpi.h"
#include "parallel.h"

using namespace std;

int main(int argc, char* argv[]) { 
    
    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    // INITIALIZE MPI ENVIRONMENT
    
    int size, rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    cout<<" I am node "<<rank<<" of the "<<size<<" you invoked"<<endl;

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    // READ RUN CARD

    Input();

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    // GET CAPITALS FROM FILE AND DISPLAY CITIES

    string filenameCAPITALS = "American_capitals.dat";
    const vector<Capital> cities = readCapitalsFromFile(filenameCAPITALS);
    int number_of_cities = cities.size();    
    cout << "File 'American_capitals.dat' read. The salesman has to visit the following "<< number_of_cities << " cities:" << endl;
    /*for (const City& city : cities) { printCity(city); } // use iterators*/
    cout << "===================================" << endl;
    
    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    // INITIALIZE POPULATION

    vector<vector<int>> Population = initializePopulation(population_size, number_of_cities, rnd);
    if (checkPopulation(Population,number_of_cities)) { cout << "Population properly initialized" << endl; } 
    else { cout << "Error during initialization! Exiting... " << endl; exit(0); }
    /* cout << "Starting population:" << endl;
    printPopulation(Population); */
    cout << "===================================" << endl;

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    // INITIALIZE VECTORS TO STORE INFO TO BE PRINTED THEN

    vector<double> bestPath; 
    vector<double> bestHalf;
    
    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    // ADVANCE GENERATIONS

    cout << "Commencing evolution process:" << endl;
    for (unsigned int gen = 1; gen<=number_of_generations; gen++ ) {

        //////////////////////////////// user update: show step
        
        if (gen%50==0) { cout << "Generation " << gen << "/" << number_of_generations << endl; }

        //////////////////////////////// mutations

        for (auto& individual : Population) { // use iterators to mutate w/ a certain probability 
            if (rnd.Rannyu()<P_pair_permutation) 
            { GenMut_PairPerm(individual,rnd,number_of_cities);             }
            if (rnd.Rannyu()<P_shift)
            { GenMut_Shift(individual,rnd,number_of_cities);                }
            if (rnd.Rannyu()<P_inversion) 
            { GenMut_Inversion(individual,rnd,number_of_cities);            }
            if (rnd.Rannyu()<P_group_permutation) 
            { GenMut_GroupsPermutation(individual, rnd, number_of_cities);  }
        }

        //////////////////////////////// population check
        
        if (!checkPopulation(Population, number_of_cities)){
            cout << "MUTATION WARNING: something went wrong at step "<< gen << endl;
            exit(0);
        }

        //////////////////////////////// prepare to kill some individuals by assigning unfitness and sorting the population
        
        vector<double> UP = calculateUnfitnessParameters(Population, cities, which_loss);
        sortPopulationByUnfitness(Population, UP);
        UP.clear(); // dispose of the useless variables to free RAM

        //////////////////////////////// sort probabilities as well
        
        vector<double> sorted_UP = calculateUnfitnessParameters(Population, cities, which_loss);

        //////////////////////////////// remove individuals based on random numbers and unfitness parameters: poulation must be sorted first!
        
        killAndReplace(Population, sorted_UP, percentage_to_kill, rnd);

        //////////////////////////////// do some crossovers
        
        for (unsigned int i = 0; i < population_size-1; i+=2)
        {
            if (rnd.Rannyu()<P_crossover)
            { doCrossover(Population[i], Population[i+1], rnd); }
            
        }

        //////////////////////////////// sort by increasing fitness / decreasing unfitness
        
        UP.clear();
        UP = calculateUnfitnessParameters(Population, cities, which_loss);
        sortPopulationByUnfitness(Population, UP);
        UP.clear(); // delete then...
        if ( !checkPopulation(Population, number_of_cities) || population_size!=Population.size() ){
            cout << "BREEDING WARNING: something went wrong at step "<< gen << endl;
        }

        //#############################################################
        // MIGRATIONS: synchronize and exchange individuals
        //#############################################################
        
        if ( BOOLmigrations && (gen%migration_interval) == 0) { // every migration_interval generations

            //int migration_step = 1+gen/migration_interval;
            //cout << "MIGRATION NUMBER "<< migration_step << endl;
            
            //////////////////////////////// prepare data to send: select random individuals

            int num_individuals_to_exchange = population_size*percentage_that_migrates; // number of individuals to exchange
            vector<int> send_buffer(num_individuals_to_exchange * number_of_cities); // this is a 'flattened' collection of individuals, all stored in a unique auxiliary row vector
            for (int i = 0; i < num_individuals_to_exchange; i++) {
                int randomIndex = static_cast<int>(rnd.Rannyu(0, population_size)); // select individual at random (uniformly)
                copy(Population[randomIndex].begin(), Population[randomIndex].end(), send_buffer.begin() + i * number_of_cities);
            }

            //////////////////////////////// exchange individuals with other processes

            vector<int> recv_buffer(num_individuals_to_exchange * number_of_cities * size); // collect individuals from all processes
            MPI_Allgather(send_buffer.data(), num_individuals_to_exchange * number_of_cities, MPI_INT,
                          recv_buffer.data(), num_individuals_to_exchange * number_of_cities, MPI_INT,
                          MPI_COMM_WORLD);
            /*           Note: these are the rules to use MPI_Allgather 
            int MPI_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                    void *recvbuf, int recvcount, MPI_Datatype recvtype,
                              MPI_Comm comm)      
            'MPI_Allgather is similar to MPI_Gather, except that all processes receive the result, instead of just the root. 
            In other words, all processes contribute to the result, and all processes receive the result.'
            See also:
            https://www.mpich.org/static/docs/latest/www3/MPI_Allgather.html
            https://docs.open-mpi.org/en/v5.0.x/man-openmpi/man3/MPI_Allgather.3.html         */
            

            //////////////////////////////// incorporate received individuals

            for (int i = 0; i < size; i++) {
                if (i != rank) {
                    for (int j = 0; j < num_individuals_to_exchange; j++) {
                        vector<int> received_individual(recv_buffer.begin() + (i * num_individuals_to_exchange + j) * number_of_cities,
                                                        recv_buffer.begin() + (i * num_individuals_to_exchange + j + 1) * number_of_cities);
                        Population.push_back(received_individual);
                    }
                }
            }

            //////////////////////////////// sort again population and do checks

            UP = calculateUnfitnessParameters(Population, cities, which_loss);
            sortPopulationByUnfitness(Population, UP);
            UP.clear();
            if (Population.size() > population_size) { Population.erase(Population.begin(), Population.end() - population_size);  }
            if (!checkPopulation(Population, number_of_cities)){ cout << "MIGRATION WARNING: something went wrong!" << endl; exit(0); }

        }

        //#############################################################
        //  end of migrations
        //#############################################################

        //////////////////////////////// calculate the average of the best half of the population
        
        double averageBestHalf = 0;
        for (unsigned int i = 0; i<population_size/2; i++ ){
            averageBestHalf = averageBestHalf*static_cast<double>(i)/static_cast<double>(i+1) + Loss(Population[population_size/2+i], cities, which_loss)/static_cast<double>(i+1);
        }
        bestHalf.push_back(averageBestHalf);

        //////////////////////////////// show fittest individual on video
        
        bestPath.push_back(Loss(Population.back(), cities, which_loss));
        if (gen%50==0) { cout << "Shortest route found = "<< bestPath[gen-1] << endl; }
    }

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    // SORT POPULATION ONE LAST TIME AND PRINT RESULTS BOTH TO VIDEO AND TO FILE

    vector<double> UP = calculateUnfitnessParameters(Population, cities, which_loss);
    sortPopulationByUnfitness(Population, UP);

    //////////////////////////////// video
    
    /*cout << "Final population" << endl;
    printPopulation(Population);*/
    
    //////////////////////////////// file: best path
    
    string outputBESTPATH = "migr" + to_string(migrations_enabled) + "node" + to_string(rank)  + "capitals" + to_string(number_of_cities) + "loss" + to_string(which_loss)+".bestpath";
    writeFittestToFile(Population.back(), cities, outputBESTPATH);
    
    //////////////////////////////// show fittest individual
    
    cout << "Node: " << to_string(rank) << "\tShortest route found = "<< Loss(Population[population_size-1], cities, which_loss) << endl;
    
    //////////////////////////////// print results
    
    string outputLOSS = "migr" + to_string(migrations_enabled) + "node" + to_string(rank) + "capitals" + to_string(number_of_cities) + "loss" + to_string(which_loss) + ".decrease";
    ofstream OUT;
    OUT.open(outputLOSS, ios::out);
    for (unsigned int i = 1; i <= number_of_generations; i++)
    { OUT << i << "\t" << bestPath[i-1] << "\t" << bestHalf[i-1] << endl; }

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    // TERMINATES MPI PROCESSING AND CLOSES PROGRAM
    
    cout << "===================================" << endl;
    cout << "===========  END RANK " << to_string(rank) << " ===========" << endl;
    cout << "===================================" << endl;
    MPI_Finalize();

    return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////
// function to read run card
//////////////////////////////////////////////////////////////
void Input(void) {
    cout << "===================================" << endl;
    cout << "======== GENETIC ALGORITHM ========" << endl;
    cout << "========        TSP        ========" << endl;
    cout << "===================================" << endl;
    //////////////////////////////////////////////////////////////
    // read seed for random numbers
    int p1, p2;
    ifstream Primes("Primes");
    if (!Primes.is_open()) {
        cerr << "Error opening Primes! " << endl;
        exit(0);
    }
    Primes >> p1 >> p2;
    Primes.close();
    ifstream input("seed.in");
    if (!input.is_open()) {
        cerr << "Error opening seed.in! " << endl;
        exit(0);
    }
    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed, p1, p2);
    input.close();
    //////////////////////////////////////////////////////////////
    // read run card
    ifstream ReadRunCard;
    ReadRunCard.open("run_card.dat");
    if (!ReadRunCard.is_open()) {
        cerr << "Error opening run_card.dat! " << endl;
        exit(0);
    }
    // population size
    ReadRunCard >> population_size;
    cout << "Population size = " << population_size << endl;
    if (population_size<=0) {
        cout << "Select valid population size!" << endl;
        cout << "Exiting..." << endl;
        exit(0);
    }
    // loss function preferred
    ReadRunCard >> which_loss;
    if (which_loss == 1 || which_loss == 2) {
        cout << "Loss function in use = " << which_loss << endl;
    } else {
        cout << "Select a valid loss function! (1 or 2)" << endl;
        cout << "Exiting..." << endl;
        exit(0);
    }
    // probabilities of mutations and crossover
    bool badProbabilities = false;
    // pair perm
    ReadRunCard >> P_pair_permutation;
    cout << "Probability of pair permutation = " << P_pair_permutation << endl;
    if (P_pair_permutation<0||P_pair_permutation>1) { badProbabilities=true; }
    // shift
    ReadRunCard >> P_shift;
    cout << "Probability of shift = " << P_shift << endl;
    if (P_shift<0||P_shift>1) { badProbabilities=true; }
    // group perm
    ReadRunCard >> P_group_permutation;
    cout << "Probability of group permutation = " << P_group_permutation << endl;
    if (P_group_permutation<0||P_group_permutation>1) { badProbabilities=true; }
    // inversion
    ReadRunCard >> P_inversion;
    cout << "Probability of inversion = " << P_inversion << endl;
    if (P_inversion<0||P_inversion>1) { badProbabilities=true; }
    // crossover
    ReadRunCard >> P_crossover;
    cout << "Probability of crossover = " << P_crossover << endl;
    if (P_crossover<0||P_crossover>1) { badProbabilities=true; }
    // check probabilities
    if (badProbabilities) {
        cout << "Select valid probabilities! (between 0 and 1)" << endl;
        cout << "Exiting..." << endl;
        exit(0);
    }
    // number of generations
    ReadRunCard >> number_of_generations;
    cout << "Number of generations = " << number_of_generations << endl;
    if (number_of_generations<=0) {
        cout << "Select valid number of generations!" << endl;
        cout << "Exiting..." << endl;
        exit(0);
    }
    // migrations enabled or not
    ReadRunCard >> migrations_enabled;
    if (migrations_enabled==1) {
        cout << "Migrations are ENABLED (migr. percentage = " << percentage_that_migrates*100 << "%)" << endl;
        BOOLmigrations = true;
    } else if (migrations_enabled==0) {
        cout << "Migrations are NOT ENABLED" << endl;
    } else {
        cout << "Select valid input for migration activation! (1=yes, 0=no)" << endl;
        exit(0);
    }
    // how often migrations should occur
    ReadRunCard >> migration_interval;
    cout << "Interval between migrations = " << migration_interval << endl;
    if (migration_interval<=0) {
        cout << "Select valid migration interval!" << endl;
        cout << "Exiting..." << endl;
        exit(0);
    }
    // closing file reading stream
    ReadRunCard.close();
    cout << "===================================" << endl;
}

//////////////////////////////////////////////////////////////
// function to read cities from file
//////////////////////////////////////////////////////////////
vector<City> readCitiesFromFile(const string& filename) {
    vector<City> cities;
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening cities file: " << filename << endl;
        exit(0);
    }
    unsigned int id;
    double x, y;
    while (file >> id >> x >> y) {
        cities.emplace_back(id, x, y); // more efficient than move constructor!
    }
    file.close();
    return cities;
}

//////////////////////////////////////////////////////////////
// function to initialize the population
//////////////////////////////////////////////////////////////
vector<vector<int>> initializePopulation(int populationSize, int numberOfCities, Random& rnd) {
    vector<vector<int>> population;
    for (int i = 0; i < populationSize; ++i) {
        vector<int> permutation;
        permutation.push_back(1); // always start with 1 at the first position
        // fill the permutation with unique integers from 2 to numberOfCities
        for (int j = 2; j <= numberOfCities; ++j) {
            bool unique = true;
            int newRandomNumber;
            // Generate a random number until it's unique
            do {
                newRandomNumber = rnd.Rannyu(2, numberOfCities + 1); // Generate a random number from 2 to numberOfCities
                if (find(permutation.begin(), permutation.end(), newRandomNumber) != permutation.end()) { // the result is the last position of the iterator; if the iterator ends up at .end(), the number was not in the list
                    unique = false; // newRandomNumber already exists in the permutation
                } else {
                    unique = true; // newRandomNumber is unique
                }
            } while (!unique);
            permutation.push_back(newRandomNumber);
        }
        population.push_back(permutation);
    }
    return population;
}

//////////////////////////////////////////////////////////////
// function to check whether an individual permutation 
// satisfies the properties of the problem
//////////////////////////////////////////////////////////////
bool checkIndividual(const vector<int>& individual, int numberOfCities) {
    // check whether the size of the individual permutation is correct
    if (static_cast<int>(individual.size()) != numberOfCities) {
        return false;
    }
    // check whether the first element is 1
    if (individual[0] != 1) {
        return false;
    }
    // check whether all integers from 2 to numberOfCities appear exactly once in the permutation
    vector<int> sortedIndividual = individual; // copy the individual permutation to sort
    sort(sortedIndividual.begin(), sortedIndividual.end()); // sort it...
    for (int i = 2; i <= numberOfCities; ++i) { // ...and see if every number is there
        if (sortedIndividual[i - 1] != i) { // easy check if vector is sorted :P
            return false;
        }
    }
    // ok, if all conditions satisfied :)
    return true;
}

//////////////////////////////////////////////////////////////
// function to check the entire population
//////////////////////////////////////////////////////////////
bool checkPopulation(const vector<vector<int>>& population, int numberOfCities) {
    for (const auto& individual : population) {
        if (!checkIndividual(individual, numberOfCities)) {
            return false;
        }
    }
    // if all individuals are okay, so is the population :):):):):)
    return true;
}

//////////////////////////////////////////////////////////////
// function to calculate loss
//////////////////////////////////////////////////////////////
template <typename T> double Loss(const vector<int>& permutation, const vector<T>& locations, const unsigned int which_loss) {
    double totalDistance = 0.0;
    int N = permutation.size();
    // take the sum of distances between consecutive locations
    for (int i = 0; i < N - 1; ++i) {
        unsigned int loc1_id = permutation[i];     // once the location is identified 
        unsigned int loc2_id = permutation[i + 1]; // by the number in the permutation...
        const T& loc1 = locations[loc1_id - 1];    // ...and the locations 'book' is parsed...
        const T& loc2 = locations[loc2_id - 1];
        // ... the positions are readily recovered!
        double distance = (which_loss == 1) 
            ? sqrt(pow(loc2.getX() - loc1.getX(), 2.0) + pow(loc2.getY() - loc1.getY(), 2.0)) 
            : pow(loc2.getX() - loc1.getX(), 2.0) + pow(loc2.getY() - loc1.getY(), 2.0);
        totalDistance += distance;
    }
    // last term: distance from the last location back to the starting location
    unsigned int first_loc_id = permutation[0];
    unsigned int last_loc_id = permutation[N - 1];
    const T& first_loc = locations[first_loc_id - 1];
    const T& last_loc = locations[last_loc_id - 1];
    totalDistance += (which_loss == 1) 
        ? sqrt(pow(last_loc.getX() - first_loc.getX(), 2.0) + pow(last_loc.getY() - first_loc.getY(), 2.0)) 
        : pow(last_loc.getX() - first_loc.getX(), 2.0) + pow(last_loc.getY() - first_loc.getY(), 2.0);
    // all terms summed
    return totalDistance;
}

//////////////////////////////////////////////////////////////
// function to print a single individual
//////////////////////////////////////////////////////////////
void printIndividual(const vector<int>& individual) {
    cout << "\t| ";
    for (int city_id : individual) {
        cout << city_id << " ";
    }
    cout << "| " << endl;
}

//////////////////////////////////////////////////////////////
// function to print the entire population
//////////////////////////////////////////////////////////////
void printPopulation(const vector<vector<int>>& population) {
    unsigned int counter = 1;
    for (const auto& individual : population) {
        cout << counter++ << "\t- th path";
        printIndividual(individual);
    }
}

//////////////////////////////////////////////////////////////
// function to print details of a single city
//////////////////////////////////////////////////////////////
void printCity(const City& city) {
    cout << "\tCity ID: " << city.getId() << ",\tX: " << city.getX() << ",\tY: " << city.getY() << endl;
}

//////////////////////////////////////////////////////////////
// function to generate a pair permutation mutation
//////////////////////////////////////////////////////////////
void GenMut_PairPerm(vector<int>& individual, Random& rnd, int numberOfCities) {
    // note: NEVER swap position 0 (there must be always 1)
    int position = rnd.Rannyu(1, numberOfCities - 1); // result: integer between 1 and numberOfCities-2
    swap(individual[position], individual[position + 1]); // swap from STL
    if (!checkIndividual(individual, numberOfCities)) { // check if the mutated individual is valid
        cerr << "Mutation PairPerm resulted in an invalid individual!" << endl;
        // exit(0);
    }
}

//////////////////////////////////////////////////////////////
// function to generate a shift mutation
//////////////////////////////////////////////////////////////
void GenMut_Shift(vector<int>& individual, Random& rnd, int numberOfCities) {
    // ensure m is less than numberOfCities - 1
    int maxM = numberOfCities - 1 - 1; 
    int m = rnd.Rannyu(1, maxM + 1); // extract block size
    int maxN = numberOfCities - m - 1;
    int n = rnd.Rannyu(1, maxN + 1); // extract shift amount
    // define block start position
    int blockStart = rnd.Rannyu(1, numberOfCities - m - n); // start position of the block
    // copy the block to be shifted
    vector<int> block(individual.begin() + blockStart, individual.begin() + blockStart + m);
    // perform the shift by inserting block at the new position and removing the original block
    individual.erase(individual.begin() + blockStart, individual.begin() + blockStart + m);
    for (int i = 0; i < m; i++)
    {
        individual.insert(individual.begin() + blockStart + n + i, block[i] );
    }
    // check if the mutation produced a valid individual
    if (!checkIndividual(individual, numberOfCities)) {
        cout << "Mutation Shift resulted in an invalid individual!" << std::endl;
        exit(0);
    }
}

//////////////////////////////////////////////////////////////
// function to generate a group permutation
//////////////////////////////////////////////////////////////
void GenMut_GroupsPermutation(vector<int>& individual, Random& rnd, int numberOfCities) {
    int maxBlockSize = numberOfCities / 2 - 1; 
    int m = rnd.Rannyu(1, maxBlockSize + 1); // extract block size
    // define two non-overlapping blocks
    int startBlock1 = rnd.Rannyu(1, numberOfCities - 2 * m); // start of first block
    int startBlock2 = rnd.Rannyu(startBlock1 + m, numberOfCities - m); // start of second block
    // ensure blocks are non-overlapping
    if (startBlock1 + m > startBlock2) {
        startBlock2 = startBlock1 + m;
    }
    // copy the blocks
    vector<int> block1(individual.begin() + startBlock1, individual.begin() + startBlock1 + m);
    //cout << "block 1 = " ;
    //printIndividual(block1);
    vector<int> block2(individual.begin() + startBlock2, individual.begin() + startBlock2 + m);
    //cout << "block 2 = " ;
    //printIndividual(block2);
    // perform the swap
    copy(block2.begin(), block2.end(), individual.begin() + startBlock1);
    copy(block1.begin(), block1.end(), individual.begin() + startBlock2);
    // check if the mutation produced a valid individual
    if (!checkIndividual(individual, numberOfCities)) {
        cerr << "Mutation PermutationContiguousCities resulted in an invalid individual!" << std::endl;
        exit(0);
    }
}

//////////////////////////////////////////////////////////////
// function to generate a inversion mutation
//////////////////////////////////////////////////////////////
void GenMut_Inversion(vector<int>& individual, Random& rnd, int numberOfCities) {
    // Define block size
    int maxBlockSize = numberOfCities - 2; // excluding the first city (and the last, otherwise, there is no space for inversion...)
    int m = rnd.Rannyu(2, maxBlockSize + 1); // extract block size
    // Define block start position
    int blockStart = rnd.Rannyu(1, numberOfCities - m); // start position of the block (1 <= blockStart <= numberOfCities-m-1)
    // Perform inversion
    for (int i = 0; i < m / 2; ++i) { // iterate over half of the block size
        swap(individual[blockStart + i], individual[blockStart + m - 1 - i]); // swap elements from the beginning and end of the block
    }
    if (!checkIndividual(individual, numberOfCities)) {
        cout << "Mutation Inversion resulted in an invalid individual!" << endl;
    }
}

//////////////////////////////////////////////////////////////
// swap genetic material between parents
//////////////////////////////////////////////////////////////
void doCrossover(vector<int> &parent1, vector<int> &parent2, Random &rnd) {
  int numberOfCities = parent1.size();
  int numberOfCities_alias = parent2.size();
  if (numberOfCities!=numberOfCities_alias)
  {
    cout << "Uncompatible parents length!" << endl;
    exit(0);
  }  
  // choose a random crossover point
  int crossoverPoint = rnd.Rannyu(1, numberOfCities-1); // not (0,numberOfCities+1), to ensure both parts are non-empty
  // create temporary vectors to store the new offspring
  vector<int> offspring1(parent1.begin(), parent1.begin() + crossoverPoint);
  vector<int> offspring2(parent2.begin(), parent2.begin() + crossoverPoint);
  // fill in the missing cities from the second part of the opposite parent
  for (int i = 0; i < numberOfCities; ++i) {
    if (find(offspring1.begin(), offspring1.end(), parent2[i]) == offspring1.end()) {
      offspring1.push_back(parent2[i]);
    }
    if (find(offspring2.begin(), offspring2.end(), parent1[i]) == offspring2.end()) {
      offspring2.push_back(parent1[i]);
    }
  }
  if (!checkIndividual(offspring1, numberOfCities) || !checkIndividual(offspring2, numberOfCities)) {
    cout << "Crossover resulted in an invalid individual!" << endl;
    exit(0);
  }
  // overwrite the parents with the new offspring
  parent1 = offspring1;
  parent2 = offspring2;
}

//////////////////////////////////////////////////////////////
// function to calculate unfitness parameters for the population
//////////////////////////////////////////////////////////////
template <typename T> vector<double> calculateUnfitnessParameters(const vector<vector<int>>& population, const vector<T>& locations, unsigned int which_loss) {
    vector<double> unfitnessParameters(population.size());
    double totalLoss = 0.0;
    for (size_t i = 0; i < population.size(); ++i) {
        unfitnessParameters[i] = Loss(population[i], locations, which_loss); // individual loss
        totalLoss += unfitnessParameters[i]; // sum to get total loss (for normalization...)
    }
    for (size_t i = 0; i < population.size(); ++i) { // normalize individual losses to get death probabilities
        unfitnessParameters[i] /= totalLoss;
    }
    return unfitnessParameters;
}

//////////////////////////////////////////////////////////////
// function to print unfitness parameters for the population
//////////////////////////////////////////////////////////////
void printUP(const vector<double>& UnfitnessParameters) {
    cout << "UP:\t( ";
    for (auto x : UnfitnessParameters) {
        cout << x << " ";
    }
    cout << ")" << endl;
}

//////////////////////////////////////////////////////////////
// function to sort the population by decreasing unfitness parameters
//////////////////////////////////////////////////////////////
void sortPopulationByUnfitness(vector<vector<int>>& population, const vector<double>& unfitnessParameters) {
    // Create a vector of indices for the population
    vector<size_t> indices(population.size());
    for (size_t i = 0; i < population.size(); ++i) {
        indices[i] = i;
    }
    // Sort the indices based on unfitness parameters in descending order
    sort(indices.begin(), indices.end(), [&](size_t a, size_t b) { return unfitnessParameters[a] > unfitnessParameters[b]; });
    // Reorder the population based on sorted indices
    vector<vector<int>> sortedPopulation(population.size());
    for (size_t i = 0; i < population.size(); ++i) {
        sortedPopulation[i] = population[indices[i]];
    }
    // Update the original population with the sorted population
    population = sortedPopulation;
}

//////////////////////////////////////////////////////////////
// function to remove individuals from a SORTED population based on random numbers and unfitness parameters
//////////////////////////////////////////////////////////////
void killAndReplace(vector<vector<int>>& Population, const vector<double>& sorted_UP, double percentage_to_kill, Random& rnd) {
    //int originalSize = Population.size();
    double cumulativeProbability = 0.0;
    int killed = 0;
    for (int i = 0; killed<percentage_to_kill*population_size; ++i) {
        double randomNumber = rnd.Rannyu();
        cumulativeProbability = 0.0;
        for (unsigned int j = 0; j < sorted_UP.size(); ++j) {
            cumulativeProbability += sorted_UP[j];
            if (randomNumber < cumulativeProbability) {
                //cout << "pulled my trigger now he's dead"<< endl;
                Population.erase(Population.begin()+j);
                //cout << "i've killed , now there are "<< Population.size() << " individuals left" << endl;
                killed++;
                Population.insert(Population.begin()+j,Population.back());
                //cout << "i've reproduced , now there are "<< Population.size() << " individuals left" << endl;
                break;
            }
        }
    }
    //cout << "i've done it, i have replaced "<< killed << " individuals" << endl;
    //cout << "now there are "<< Population.size() << " individuals left" << endl;
    //exit(0);
}

//////////////////////////////////////////////////////////////
// function to write result to file
//////////////////////////////////////////////////////////////
template <typename T> void writeFittestToFile(const vector<int>& fittest, const vector<T>& locations, const string& filename) {
    ofstream outFile(filename);
    cout << "Writing fittest individual to file "<< filename << endl;
    if (!outFile) {
        cerr << "Error opening file " << filename << " for writing!" << endl;
        return;
    }
    for (int cityIndex : fittest) {
        outFile << cityIndex << "\t";
    }
    outFile.close();
}

//////////////////////////////////////////////////////////////
// function to read capitals from file
//////////////////////////////////////////////////////////////
vector<Capital> readCapitalsFromFile(const string& filename) {
    vector<Capital> capitals;
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        exit(1);
    }
    string line;
    getline(file, line); // Skip the header line
    unsigned int id = 1;
    while (getline(file, line)) {
        istringstream iss(line);
        string state, capitalName;
        double longitude, latitude;
        // Reading state and capital name can be tricky if they contain spaces, hence using getline with delimiters
        getline(iss >> ws, state, '\t');
        getline(iss >> ws, capitalName, '\t');
        iss >> longitude >> latitude;
        // Trim any potential leading or trailing spaces (tab or space)
        state.erase(remove_if(state.begin(), state.end(), ::isspace), state.end());
        capitalName.erase(remove_if(capitalName.begin(), capitalName.end(), ::isspace), capitalName.end());
        capitals.emplace_back(id++, longitude, latitude, state, capitalName);
    }
    file.close();
    return capitals;
}