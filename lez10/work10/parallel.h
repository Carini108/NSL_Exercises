#ifndef __GA__
#define __GA__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <cmath>

using namespace std;

//////////////////////////////////////////////////////////////
// random numbers
#include "random.h"
int seed[4];
Random rnd;

//////////////////////////////////////////////////////////////
// genetic algorithm parameters
double P_pair_permutation, P_shift, P_group_permutation, P_inversion, P_crossover;
const double percentage_to_kill = 0.99;
const double percentage_that_migrates = 0.5;
unsigned int population_size, which_loss, number_of_generations;
int migrations_enabled, migration_interval;
bool BOOLmigrations = false;

//////////////////////////////////////////////////////////////
// class for cities
class City {
private:
    unsigned int id;
    double x;
    double y;
public:
    City(unsigned int id = 0, double x = 0.0, double y = 0.0) : id(id), x(x), y(y) {}
    unsigned int getId() const { return id; }
    void setId(unsigned int newId) { id = newId; }
    double getX() const { return x; }
    void setX(double newX) { x = newX; }
    double getY() const { return y; }
    void setY(double newY) { y = newY; }
};

class Capital : public City {
private:
    string state;
    string capitalName;

public:
    Capital(unsigned int id = 0, double x = 0.0, double y = 0.0, const string& state = "", const string& capitalName = "")
        : City(id, x, y), state(state), capitalName(capitalName) {}

    string getState() const { return state; }
    void setState(const string& newState) { state = newState; }

    string getCapitalName() const { return capitalName; }
    void setCapitalName(const string& newCapitalName) { capitalName = newCapitalName; }
};

//////////////////////////////////////////////////////////////
// functions
void Input(void);
vector<City> readCitiesFromFile(const string& filename);
vector<Capital> readCapitalsFromFile(const string& filename);
vector<vector<int>> initializePopulation(int populationSize, int numberOfCities, Random& rnd);
bool checkIndividual(const vector<int>& individual, int numberOfCities);
bool checkPopulation(const vector<vector<int>>& population, int numberOfCities);
template <typename T> double Loss(const vector<int>& permutation, const vector<T>& locations, const unsigned int which_loss);
void printIndividual(const vector<int>& individual);
void printPopulation(const vector<vector<int>>& population);
void printCity(const City& city);
// 
void GenMut_PairPerm(vector<int>& individual, Random& rnd, int numberOfCities);
void GenMut_Shift(vector<int>& individual, Random& rnd, int numberOfCities);
void GenMut_GroupsPermutation(vector<int>& individual, Random& rnd, int numberOfCities);
void GenMut_Inversion(vector<int>& individual, Random& rnd, int numberOfCities);
void doCrossover(vector<int>& parent_1, vector<int>& parent_2, Random& rnd);
// 
template <typename T> vector<double> calculateUnfitnessParameters(const vector<vector<int>>& population, const vector<T>& locations, unsigned int which_loss);
void printUP(const vector<double>& UnfitnessParameters);
void sortPopulationByUnfitness(vector<vector<int>>& population, const vector<double>& unfitnessParameters);
void killAndReplace(vector<vector<int>>& Population, const vector<double>& sorted_UP, double percentage_to_kill, Random& rnd);
template <typename T> void writeFittestToFile(const vector<int>& fittest, const vector<T>& locations, const string& filename);

#endif