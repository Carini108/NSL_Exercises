#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{ 
  Input(); //Inizialization
  int n_of_plt_points = 100; // 100 points in the interval from T_max to T_min
  for (int n_th_point = n_of_plt_points; n_th_point >= 1; n_th_point--) // start from highest T
  {  
    for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
    {
      Reset(iblk);   //Reset block averages
      for(int istep=1; istep <= nstep; ++istep)
      {
        //cout << "\t step = " << istep << endl;
        Move(metro);
        Measure();
        Accumulate(); //Update block averages
      }
      Averages(iblk);   //Print results for current block
      if (iblk==nblk)
      {
        SavePlotPointAndUpdateTemp(n_th_point,n_of_plt_points); //Save only last block
      }
    }
    //ConfFinal(); //Uncomment to write final configuration
  }
  return 0;
}

void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature upper bound = " << temp << endl;

  ReadInput >> finaltemp;
  cout << "Temperature lower bound = " << finaltemp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> restart;

  if(metro==1) cout << "The program performs Metropolis moves" << endl;
  else cout << "The program performs Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

//Prepare arrays for measurements (4 observables)
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration
  if ( restart==0 ) { // if 0, do not restart and produce random config
    cout << "restart 0: starting from scratch (random config)" << endl;
    for (int i=0; i<nspin; ++i)
    {
      if(rnd.Rannyu() >= 0.5) s[i] = 1;
      else s[i] = -1;
    }
  } else if ( restart==1 ) { // if 1, restart from config.final
    cout << "restart 1: starting from last config (in 'config.final')" << endl;
    ifstream restart_config("config.final");
    for (int i=0; i<nspin; ++i)
    {
      restart_config >> s[i];
    }
  } else {
    cout << "ERROR: not valid restart option! Exiting..." << endl;
    exit(0);
  }
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}

void Move(int metro)
{
  int o;
  double p, energy_old, energy_new, sm;
  double energy_up, energy_down;

  for(int i=0; i<nspin; ++i) // try with all spins
  {
  //Select randomly a particle to flip (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1) //Metropolis
    { 
      if ( rnd.Rannyu()<min(1.,exp(-beta*( Boltzmann(-s[o],o)-Boltzmann(s[o], o) ))) )  
      {
        s[o]=-s[o];
        accepted+=1.; // sometimes rejected
      }
    }
    else //Gibbs sampling
    { 
      accepted+=1.; // always accepted
      if ( rnd.Rannyu()<pow( (1.+exp(-beta*( Boltzmann(+1,o)-Boltzmann(-1,o) ))) , -1. ) ) // slide 21
      {
        s[o]=-1;
      }
      else 
      {
        s[o]=+1;
      }
    }
    attempted+=1.;
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  int bin;
  double u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     m += s[i];
  }
  walker[iu] = u;
  walker[im] = m;
  walker[ix] = beta*m*m;
  walker[ic] = u*u;
}

void Reset(int iblk) //Reset block averages
{
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}

void Accumulate(void) //Update block averages
{
   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}

void Averages(int iblk) //Print results for current block
{
   ofstream Ene, Heat, Mag, Chi;
   const int wd=15;
   
   // Determine the suffixes for filenames
   string field_suffix = (h == 0) ? "fieldOFF" : "fieldON";
   string method_suffix = (metro == 1) ? "METRO" : "GIBBS";
   
   if (iblk%100==0) {
      cout << "Block number " << iblk << "/ " << nblk << endl;
      cout << "Acceptance rate " << accepted/attempted << endl << endl;
      cout << "----------------------------" << endl << endl;
   }
   
   // Construct filenames with the suffixes
   string ene_filename = "output.ene.0.isi." + field_suffix + "." + method_suffix;
   string mag_filename = "output.mag.0.isi." + field_suffix + "." + method_suffix;
   string heat_filename = "output.heat.0.isi." + field_suffix + "." + method_suffix;
   string chi_filename = "output.chi.0.isi." + field_suffix + "." + method_suffix;
   
   //------------------------------------
   Ene.open(ene_filename, ios::app); // append (add and do not delete already existing data on file)
   stima_u       = blk_av[iu]/blk_norm/(double)nspin; //Energy
   glob_av[iu]  += stima_u;
   glob_av2[iu] += stima_u*stima_u;
   err_u         = Error(glob_av[iu],glob_av2[iu],iblk);
   //                n blk                 instant                 blk progression ave                   err progress
   Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
   Ene.close();
   //------------------------------------
   Mag.open(mag_filename, ios::app);
   stima_m       = blk_av[im]/blk_norm/(double)nspin; //Magnetization
   glob_av[im]  += stima_m;
   glob_av2[im] += stima_m*stima_m;
   err_m         = Error(glob_av[im],glob_av2[im],iblk);
   //                n blk                 instant                 blk progression ave                   err progress
   Mag << setw(wd) << iblk << setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
   Mag.close();
   //------------------------------------
   Heat.open(heat_filename, ios::app);
   stima_c       = k_B*beta*beta*( blk_av[ic]/blk_norm - pow( blk_av[iu]/blk_norm, 2 ) ) / (double)nspin; //Specific heat
   glob_av[ic]  += stima_c;
   glob_av2[ic] += stima_c*stima_c;
   err_c         = Error(glob_av[ic],glob_av2[ic],iblk);
   //                n blk                 instant                 blk progression ave                   err progress
   Heat << setw(wd) << iblk << setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
   Heat.close();
   //------------------------------------
   Chi.open(chi_filename, ios::app);
   stima_x       = blk_av[ix]/blk_norm/(double)nspin; //Susceptibility
   glob_av[ix]  += stima_x;
   glob_av2[ix] += stima_x*stima_x;
   err_x         = Error(glob_av[ix],glob_av2[ix],iblk);
   //                n blk                 instant                 blk progression ave                   err progress
   Chi << setw(wd) << iblk << setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x<< endl;
   Chi.close();
   //------------------------------------
}

void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Printing final configuration to file config.final... " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

void SavePlotPointAndUpdateTemp(int n_th_point, int n_of_plt_points)
{
  ofstream Ene, Heat, Mag, Chi;
  const int wd=15;

  string field_suffix = (h == 0) ? "fieldOFF" : "fieldON";
  string method_suffix = (metro == 1) ? "METRO" : "GIBBS";

  string plot_errorbar_ene_filename = "plot_errorbar_ene.0.isi." + field_suffix + "." + method_suffix;
  string plot_errorbar_mag_filename = "plot_errorbar_mag.0.isi." + field_suffix + "." + method_suffix;
  string plot_errorbar_heat_filename = "plot_errorbar_heat.0.isi." + field_suffix + "." + method_suffix;
  string plot_errorbar_chi_filename = "plot_errorbar_chi.0.isi." + field_suffix + "." + method_suffix;

  double temperature_point = finaltemp+(double)n_th_point/(double)n_of_plt_points*(temp-finaltemp);
  double t  = finaltemp+(double)(n_th_point-1)/(double)(n_of_plt_points)*(temp-finaltemp); //
  beta = 1./t;

  //------------------------------------
  Ene.open(plot_errorbar_ene_filename, ios::app);
  err_u         = Error(glob_av[iu],glob_av2[iu],nblk);
  //                n th point of plot         last blk progression ave              err progress
  Ene << setw(wd) << temperature_point << setw(wd) << glob_av[iu]/(double)nblk << setw(wd) << err_u << endl;
  Ene.close();
  //------------------------------------
  Mag.open(plot_errorbar_mag_filename, ios::app);
  err_m         = Error(glob_av[im],glob_av2[im],nblk);
  //                n th point of plot         last blk progression ave              err progress
  Mag << setw(wd) << temperature_point << setw(wd) <<  glob_av[im]/(double)nblk << setw(wd) << err_m << endl;
  Mag.close();
  //------------------------------------
  Heat.open(plot_errorbar_heat_filename, ios::app);
  err_c         = Error(glob_av[ic],glob_av2[ic],nblk);
  //                n th point of plot         last blk progression ave              err progress
  Heat << setw(wd) << temperature_point << setw(wd) << glob_av[ic]/(double)nblk << setw(wd) << err_c << endl;
  Heat.close();
  //------------------------------------
  Chi.open(plot_errorbar_chi_filename, ios::app);
  err_x         = Error(glob_av[ix],glob_av2[ix],nblk);
  //                n th point of plot         last blk progression ave              err progress
  Chi << setw(wd) << temperature_point << setw(wd) << glob_av[ix]/(double)nblk << setw(wd) << err_x<< endl;
  Chi.close();
  //------------------------------------
  cout << "==============================" << endl;
  cout << "PLOT COMPLETION: "<< static_cast<double>(n_of_plt_points-n_th_point)/static_cast<double>(n_of_plt_points)*100 << "%" << endl;
  cout << "==============================" << endl;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/