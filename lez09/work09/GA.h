#ifndef __GA__
#define __GA__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <string>
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
unsigned int population_size, which_loss;

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

//////////////////////////////////////////////////////////////
// functions
void Input(void);
vector<City> readCitiesFromFile(const string& filename);
vector<vector<int>> initializePopulation(int populationSize, int numberOfCities, Random& rnd);
bool checkIndividual(const vector<int>& individual, int numberOfCities);
bool checkPopulation(const vector<vector<int>>& population, int numberOfCities);
double Loss(const vector<int>& permutation, const vector<City>& cities, const unsigned int which_loss);
void printIndividual(const vector<int>& individual);
void printPopulation(const vector<vector<int>>& population);
void printCity(const City& city);
//void pair_permutation();
//void shift();
//void group_permutation();
//void inversion();
//void crossover();

#endif