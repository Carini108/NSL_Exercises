#include <iostream>
#include <fstream>
#include <ostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <math.h>
#include <numeric>
#include <thread> // For sleep_for
#include <chrono> // For chrono::seconds
#include "MD_MC.h"

using namespace std;

int main(int argc, char* argv[]) {
    if (argc != 2) { cerr << "Usage: " << argv[0] << " <solid|liquid|gas>" << endl; return 1; }
    state = argv[1];
    filename = "input_" + state + ".in";
    if (state != "solid" && state != "liquid" && state != "gas") { cerr << "Error: Invalid argument. Must be 'solid', 'liquid', or 'gas'." << endl; return 1; }
  Input(); //Inizialization
  int nconf = 1;
  for(int iblk=1; iblk <= nblk; iblk++) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; istep++)
    {
      Move();
      Measure();
      Accumulate(); //Update block averages
      //if(istep%10 == 0){
//        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
      //}
    }
    Averages(iblk);   //Print results for current block
    NormalizeAndSaveBlock(); //Normalize and save the current block's RDF
  }
  ConfFinal(); //Write final configuration
  CalculateFinalRDF();  // Calculate final RDF by averaging over all blocks and computing errors

  string filename_rdf_results = "RDF_RESULTS_"+state+".dat";
  string filename_rdf_block = "RDF_BLOCKS_"+state+".dat";
  PrintRDFResults(filename_rdf_results ,filename_rdf_block);
  
  return 0;
}


void Input(void)
{
  ifstream ReadInput, ReadConf, ReadVelocity, Primes, Seed;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "MD(NVE) / MC(NVT) simulation       " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

//Read seed for random numbers
  int p1, p2;
  Primes.open("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

//Read input informations
  ReadInput.open("input_" + state + ".in");
  if (!ReadInput.is_open()) {
      cerr << "Error: Cannot open input_" + state + ".in file" << endl;
      exit(1);
  }

  ReadInput >> iNVET;
  ReadInput >> restart;

  if(restart) Seed.open("seed.out");
  else Seed.open("seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  Seed.close();

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  box = pow(vol,1.0/3.0); // "box" is the side of the cubic box 
  bin_size = (box/2.0) / static_cast<double>(N_bins);
  cout << "Volume of the simulation box = " << vol << endl;
  cout << "Edge of the simulation box = " << box << endl;
  cout << "Bin size for g(r) histogram = " << bin_size << endl;

  ReadInput >> rcut;
  cout << "Cutoff of the interatomic potential = " << rcut << endl << endl;
    
  ReadInput >> delta;

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> phase;
  if (phase==1) { state = "solid"; }
  else if (phase==2) { state = "liquid"; }
  else if (phase==3) { state = "gas"; }
  else { cerr << "enter valid phase" << endl; exit(6); }
  

  // tail corrections
  tail_v = (8.*M_PI/3.)*(rho/pow(rcut,3.))*(1./(3.*pow(rcut,6.))-1.);
  tail_p = (32.*M_PI*static_cast<double>(npart))*(rho/pow(rcut,3.))*(1./(3.*pow(rcut,6.))-0.5);

  cout << "The program perform Metropolis moves with uniform translations" << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl;
  cout << "Tail corrections enabled"<< endl << endl;
  ReadInput.close();

//Prepare arrays for measurements
  iv = 0; //Potential energy
  it = 1; //Temperature
  ik = 2; //Kinetic energy
  ie = 3; //Total energy
  ip = 4; //Pressure
  n_props = 5; //Number of observables 

//Read initial configuration
  cout << "Read initial configuration" << endl << endl;
  if(restart)
  {
    cout << "Restarting from last configuration saved " << endl;
    ReadConf.open("config.out");
    ReadVelocity.open("velocity.out");
    for (int i=0; i<npart; ++i) ReadVelocity >> vx[i] >> vy[i] >> vz[i];
  }
  else 
  {
    ReadConf.open("config.in");
    cout << "Prepare velocities with center of mass velocity equal to zero " << endl;
    double sumv[3] = {0.0, 0.0, 0.0};
    for (int i=0; i<npart; ++i)
    {
      vx[i] = rnd.Gauss(0.,sqrt(temp));
      vy[i] = rnd.Gauss(0.,sqrt(temp));
      vz[i] = rnd.Gauss(0.,sqrt(temp));
      sumv[0] += vx[i];
      sumv[1] += vy[i];
      sumv[2] += vz[i];
    }
    for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
    double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i)
    {
      vx[i] = vx[i] - sumv[0];
      vy[i] = vy[i] - sumv[1];
      vz[i] = vz[i] - sumv[2];
      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }
    sumv2 /= (double)npart;
    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
    cout << "velocity scale factor: " << fs << endl << endl;
    for (int i=0; i<npart; ++i)
    {
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;
    }
  }

  for (int i=0; i<npart; ++i)
  {
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = Pbc( x[i] * box );
    y[i] = Pbc( y[i] * box );
    z[i] = Pbc( z[i] * box );
  }
  ReadConf.close();

  for (int i=0; i<npart; ++i)
  {
    if(iNVET)
    {
      xold[i] = x[i];
      yold[i] = y[i];
      zold[i] = z[i];
    }
    else
    {
      xold[i] = Pbc(x[i] - vx[i] * delta);
      yold[i] = Pbc(y[i] - vy[i] * delta);
      zold[i] = Pbc(z[i] - vz[i] * delta);
    }
  }
  
//Evaluate properties of the initial configuration
  Measure();
  
//Print initial values for measured properties
  cout << "Initial potential energy = " << walker[iv]/(double)npart << endl;
  cout << "Initial temperature      = " << walker[it] << endl;
  cout << "Initial kinetic energy   = " << walker[ik]/(double)npart << endl;
  cout << "Initial total energy     = " << walker[ie]/(double)npart << endl;
  cout << "Initial pressure         = " << walker[ip] << endl;
  cout << "Input was successful. Waiting for 5 seconds..." << endl;
  this_thread::sleep_for(chrono::seconds(5));
  return;
}


void Move()
{
  int o;
  double p, energy_old, energy_new;
  double xnew, ynew, znew;

  if(iNVET) // Monte Carlo (NVT) move
  {
    for(int i=0; i<npart; ++i)
    {
    //Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)
      o = (int)(rnd.Rannyu()*npart);

    //Old
      energy_old = Boltzmann(x[o],y[o],z[o],o);

    //New
      x[o] = Pbc( x[o] + delta*(rnd.Rannyu() - 0.5) );
      y[o] = Pbc( y[o] + delta*(rnd.Rannyu() - 0.5) );
      z[o] = Pbc( z[o] + delta*(rnd.Rannyu() - 0.5) );

      energy_new = Boltzmann(x[o],y[o],z[o],o);

    //Metropolis test
      p = exp(beta*(energy_old-energy_new));
      if(p >= rnd.Rannyu())  
      {
      //Update
        xold[o] = x[o];
        yold[o] = y[o];
        zold[o] = z[o];
        accepted = accepted + 1.0;
        //cout << "accepted" <<endl;
      } else {
        x[o] = xold[o];
        y[o] = yold[o];
        z[o] = zold[o];
        //cout << "rejected" <<endl;
      }
      attempted = attempted + 1.0;
    }
  } else // Molecular Dynamics (NVE) move
  {
    double fx[m_part], fy[m_part], fz[m_part];

    for(int i=0; i<npart; ++i){ //Force acting on particle i
      fx[i] = Force(i,0);
      fy[i] = Force(i,1);
      fz[i] = Force(i,2);
    }

    for(int i=0; i<npart; ++i){ //Verlet integration scheme

      xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
      ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
      znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

      vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
      vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
      vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

      xold[i] = x[i];
      yold[i] = y[i];
      zold[i] = z[i];

      x[i] = xnew;
      y[i] = ynew;
      z[i] = znew;

      accepted = accepted + 1.0;
      attempted = attempted + 1.0;
    }
  }
  return;
}

double Boltzmann(double xx, double yy, double zz, int ip)
{
  double ene=0.0;
  double dx, dy, dz, dr;

  for (int i=0; i<npart; ++i)
  {
    if(i != ip)
    {
// distance ip-i in pbc
      dx = Pbc(xx - x[i]);
      dy = Pbc(yy - y[i]);
      dz = Pbc(zz - z[i]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut)
      {
        ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      }
    }
  }

  return 4.0*ene;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
 
  return f;
}

void Measure() //Properties measurement
{
  double v = 0.0, kin=0.0, pressure=0.0;
  double vij;
  double dx, dy, dz, dr;
  std::fill(rdf_hist.begin(), rdf_hist.end(), 0.0);

// cycle over pairs of particles
  for (int i=0; i<npart-1; ++i)
  {
    for (int j=i+1; j<npart; ++j)
    {
// distance i-j in pbc
      dx = Pbc(x[i] - x[j]);
      dy = Pbc(y[i] - y[j]);
      dz = Pbc(z[i] - z[j]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut)
      {
        vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
        v += vij;
      }
      if (dr < 0.5 * box) // fill g(r) histo
      {
        int bin = int( dr / bin_size);
        rdf_hist[bin] += 2.0; // Increment by 2 for each pair (i, j) and (j, i)
      }
    }          
  }

  for (int i=0; i<npart; ++i) kin += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

  for (int i=0; i<npart-1; ++i) {
      for (int j=i+1; j<npart; ++j) {
        dx = Pbc(x[i] - x[j]);
        dy = Pbc(y[i] - y[j]);
        dz = Pbc(z[i] - z[j]);
        dr = dx*dx + dy*dy + dz*dz;
        dr = sqrt(dr);

        pressure = pressure + pow(dr,-12.)-0.5*pow(dr,-6.);
      }
  }


  walker[iv] = 4.0 * v; // Potential energy
  walker[ik] = kin; // Kinetic energy
  walker[it] = (2.0 / 3.0) * kin/(double)npart; // Temperature
  walker[ie] = 4.0 * v + kin;  // Total energy;
  walker[ip] = rho*temp + (16./vol)*pressure; // Pressure

  return;
}

void NormalizeAndSaveBlock() {
    std::vector<double> rdf_blk(N_bins, 0.0);
    for (int i = 0; i < N_bins; ++i) {
        double r_lower = i * bin_size;
        double r_upper = (i + 1) * bin_size;
        double shell_volume = (4.0 / 3.0) * M_PI * (std::pow(r_upper, 3) - std::pow(r_lower, 3)); // function of r!
        rdf_blk[i] = rdf_hist[i] / (npart * rho * shell_volume);
    }

    rdf_blocks.push_back(rdf_blk); // fill the vector of vectors
}

void CalculateFinalRDF() {
    for (int i = 0; i < N_bins; ++i) {

        std::vector<double> block_values;
        for (int block = 0; block < nblk; ++block) { block_values.push_back(rdf_blocks[block][i]); }

        // Calculate average
        double avg = std::accumulate(block_values.begin(), block_values.end(), 0.0) / nblk;

        // Calculate error
        double error = 0.0; // sum the squared differences
        for (const double& value : block_values) { error += (value - avg) * (value - avg); }
        error = std::sqrt(error / (nblk*(nblk - 1)));
        rdf_results[i] = {avg, error}; // Store final RDF data
    }
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i) // observables 
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
   // note: in principle, the blk scheme above could have been applied to g(r) as well,
   //       but it is not very convenient, since this would imply adding to n_props 
   //       other N_bins observables...
   
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   ofstream Epot, Ekin, Etot, Temp, Pres;
   const int wd=12;
    
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    Epot.open("output_epot_"+state+".dat",ios::app);
    Ekin.open("output_ekin_"+state+".dat",ios::app);
    Temp.open("output_temp_"+state+".dat",ios::app);
    Etot.open("output_etot_"+state+".dat",ios::app);
    Pres.open("output_pres_"+state+".dat",ios::app);
    
    stima_pot = blk_av[iv]/blk_norm/(double)npart; //Potential energy
    //stima_pot += tail_v/blk_norm;
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot*stima_pot;
    err_pot=Error(glob_av[iv],glob_av2[iv],iblk);
    
    stima_kin = blk_av[ik]/blk_norm/(double)npart; //Kinetic energy
    glob_av[ik] += stima_kin;
    glob_av2[ik] += stima_kin*stima_kin;
    err_kin=Error(glob_av[ik],glob_av2[ik],iblk);

    stima_etot = blk_av[ie]/blk_norm/(double)npart; //Total energy
    glob_av[ie] += stima_etot;
    glob_av2[ie] += stima_etot*stima_etot;
    err_etot=Error(glob_av[ie],glob_av2[ie],iblk);

    stima_temp = blk_av[it]/blk_norm; //Temperature
    glob_av[it] += stima_temp;
    glob_av2[it] += stima_temp*stima_temp;
    err_temp=Error(glob_av[it],glob_av2[it],iblk);

    stima_pres = blk_av[ip]/blk_norm; //Pressure
    //stima_pres += tail_p/blk_norm;
    glob_av[ip] += stima_pres;
    glob_av2[ip] += stima_pres*stima_pres;
    err_pres=Error(glob_av[ip],glob_av2[ip],iblk);

//Potential energy per particle
    Epot << setw(wd) << iblk <<  setw(wd) << stima_pot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_pot << endl;
//Kinetic energy
    Ekin << setw(wd) << iblk <<  setw(wd) << stima_kin << setw(wd) << glob_av[ik]/(double)iblk << setw(wd) << err_kin << endl;
//Total energy
    Etot << setw(wd) << iblk <<  setw(wd) << stima_etot << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_etot << endl;
//Temperature
    Temp << setw(wd) << iblk <<  setw(wd) << stima_temp << setw(wd) << glob_av[it]/(double)iblk << setw(wd) << err_temp << endl;
//Pressure
    Pres << setw(wd) << iblk <<  setw(wd) << stima_pres << setw(wd) << glob_av[ip]/(double)iblk << setw(wd) << err_pres << endl;

    
// line for clarity
    cout << "----------------------------" << endl << endl;

    Epot.close();
    Ekin.close();
    Etot.close();
    Temp.close();
    Pres.close();
}


void ConfFinal(void)
{
  ofstream WriteConf, WriteVelocity, WriteSeed;

  cout << "Print final configuration to file config.out" << endl << endl;
  WriteConf.open("config.out");
  WriteVelocity.open("velocity.out");
  for (int i=0; i<npart; ++i)
  {
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
    WriteVelocity << vx[i] << "   " <<  vy[i] << "   " << vz[i] << endl;
  }
  WriteConf.close();
  WriteVelocity.close();

  rnd.SaveSeed();
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r)  //Algorithm for periodic boundary conditions with side L=box
{
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt(fabs(sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

void PrintRDFResults(const string& filename_rdf_results, const string& filename_rdf_blocks) {
    const int wd=12;
    // Open file for rdf_results
    std::ofstream file_rdf_results(filename_rdf_results);
    if (!file_rdf_results.is_open()) {
        std::cerr << "Error opening file: " << filename_rdf_results << std::endl;
        exit(0);
    }
    // Print rdf_results to file
    for (unsigned int i = 0; i < rdf_results.size(); ++i) {
        //file_rdf_results << "Bin " << i << ": g(r) = " << rdf_results[i].g_r_avg << " +/- " << rdf_results[i].g_r_err << std::endl;
        file_rdf_results << i*bin_size << setw(wd) << rdf_results[i].g_r_avg << setw(wd) << rdf_results[i].g_r_err << endl; 
    }
    // Close file for rdf_results
    file_rdf_results.close();
    // Open file for rdf_blocks
    std::ofstream file_rdf_blocks(filename_rdf_blocks);
    if (!file_rdf_blocks.is_open()) {
        std::cerr << "Error opening file: " << filename_rdf_blocks << std::endl;
        exit(0);
    }
    // Print rdf_blocks to file
    for (unsigned int i = 0; i < rdf_blocks.size(); ++i) {
        for (unsigned int j = 0; j < rdf_blocks[i].size(); ++j) {
            //file_rdf_blocks << "Block " << i << ", Bin " << j << ": g(r) = " << rdf_blocks[i][j] << std::endl;
            file_rdf_blocks << i << setw(wd) << j*bin_size << setw(wd) << rdf_blocks[i][j] << endl;
        }
    }
    // Close file for rdf_blocks
    file_rdf_blocks.close();
}
