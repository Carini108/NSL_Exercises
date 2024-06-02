#ifndef __fluid__
#define __fluid__

#include <string>
//Random numbers
#include "random.h"
using namespace std;
int seed[4];
Random rnd;

//parameters, observables
const int m_props=1000;
int n_props, iv, ik, it, ie, iw, ip;
double vtail, ptail;
double walker[m_props];
double tail_v, tail_p; // tail corrections
int phase;
string state, filename;

struct RDFData {
    double g_r_avg;
    double g_r_err;
};

// radial distribution function
const int N_bins = 500;  // Number of bins for g(r)
double bin_size = 0.;
std::vector<double> rdf_hist(N_bins, 0.0);        // histogram for RDF
std::vector<std::vector<double>> rdf_blocks;      // RDF values for each block
std::vector<RDFData> rdf_results(N_bins);         // final RDF results with errors

// thermodynamical state
int npart;
double beta,temp,energy,vol,rho,box,rcut;

// averages
double blk_av[m_props], blk_norm, accepted, attempted;
double glob_av[m_props], glob_av2[m_props];
double stima_pot, stima_pres, stima_kin, stima_etot, stima_temp;
double err_pot, err_pres, err_kin, err_etot, err_temp;


//configuration 
const int m_part=108;
double x[m_part],    y[m_part],    z[m_part];
double xold[m_part], yold[m_part], zold[m_part];
double vx[m_part],  vy[m_part],   vz[m_part];

// simulation
int iNVET, nstep, nblk, restart;
double delta;

//pigreco
const double pi=3.1415927;

//functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
void NormalizeAndSaveBlock();
void CalculateFinalRDF();
void PrintRDFResults(const string&, const string&);
double Boltzmann(double, double, double, int);
double Pbc(double);
double Error(double,double,int);
double Force(int, int);

#endif