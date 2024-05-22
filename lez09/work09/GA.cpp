#include "GA.h"

using namespace std;

int main(int argc, char* argv[]) { 

    //////////////////////////////////////////////////////////////
    // initialization
    Input();

    //////////////////////////////////////////////////////////////
    // get filename from command line arguments and do a check
    if (argc != 2) {
        cout << "Usage: " << argv[0] << " <filename>" << endl;
        return 1;
    }
    string filename = argv[1];
    vector<City> cities = readCitiesFromFile(filename);
    int number_of_cities = cities.size();

    //////////////////////////////////////////////////////////////
    // display cities
    cout << "The salesman has to visit the following "<< number_of_cities << " cities:" << endl;
    for (const City& city : cities) { // use iterators
        printCity(city);
    }

    //////////////////////////////////////////////////////////////
    // initialize population
    vector<vector<int>> Population = initializePopulation(population_size, number_of_cities, rnd);

    //////////////////////////////////////////////////////////////
    // consistency checks (uncomment if needed) /**/
    if(checkIndividual(Population[0],number_of_cities)) {
        cout << "Alles klar Herr Kommissar (Individuum)" << endl;
    } else {
        cout << "nein!" << endl;
    }
    if (checkPopulation(Population,number_of_cities))
    {
        cout << "Alles klar Herr Kommissar (Bevölkerung)" << endl;
    } else {
        cout << "nein!" << endl;
    }
    cout << "Am Anfang ist die Bevölkerung wie folgt:" << endl;
    printPopulation(Population);

    return 0;
}


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
    // invers
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
    if (individual.size() != numberOfCities) {
        return false;
    }
    // check whether the first element is 1
    if (individual[0] != 1) {
        return false;
    }
    // check whether all integers from 2 to numberOfCities appear exactly once in the permutation
    vector<int> sortedIndividual = individual; // copy the individual permutation to sort
    sort(sortedIndividual.begin(), sortedIndividual.end()); 
    for (int i = 2; i <= numberOfCities; ++i) {
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
double Loss(const vector<int>& permutation, const vector<City>& cities, const unsigned int which_loss) {
    double totalDistance = 0.0;
    int N = permutation.size();
    // take the sum of distances between consecutive cities
    for (int i = 0; i < N - 1; ++i) {
        unsigned int city1_id = permutation[i];     // once the city is identified 
        unsigned int city2_id = permutation[i + 1]; //by the number in the permutation...
        const City& city1 = cities[city1_id];       //...and the cities 'book' is parsed...
        const City& city2 = cities[city2_id];
        // ... the positions are readily recovered!
        double distance = (which_loss == 1) ? (sqrt(pow(city2.getX() - city1.getX(), 2.) + pow(city2.getY() - city1.getY(), 2.))) : (pow(city2.getX() - city1.getX(), 2.) + pow(city2.getY() - city1.getY(), 2.));
        
        totalDistance += distance;
    }
    // last term: distance from the last city back to the starting city
    unsigned int first_city_id = permutation[0];
    unsigned int last_city_id = permutation[N - 1];
    const City& first_city = cities[first_city_id];
    const City& last_city = cities[last_city_id];

    totalDistance += (which_loss == 1) ? (sqrt(pow(last_city.getX() - first_city.getX(), 2.) + pow(last_city.getY() - first_city.getY(), 2.))) : (pow(last_city.getX() - first_city.getX(), 2.) + pow(last_city.getY() - first_city.getY(), 2.));
    // all terms summed
    return totalDistance;
}

//////////////////////////////////////////////////////////////
// function to print a single individual
//////////////////////////////////////////////////////////////
void printIndividual(const vector<int>& individual) {
    cout << "\t| ";
    for (int city : individual) {
        cout << city << " ";
    }
    cout << "| " << endl;
}

//////////////////////////////////////////////////////////////
// function to print the entire population
//////////////////////////////////////////////////////////////
void printPopulation(const vector<vector<int>>& population) {
    for (const auto& individual : population) {
        printIndividual(individual);
    }
}

//////////////////////////////////////////////////////////////
// function to print details of a single city
//////////////////////////////////////////////////////////////
void printCity(const City& city) {
    cout << "\tCity ID: " << city.getId() << ",\tX: " << city.getX() << ",\tY: " << city.getY() << endl;
}