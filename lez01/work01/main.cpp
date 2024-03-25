/****************************************************************
*****************************************************************
LSN_Exercises_01
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "./prng/random.h"
#include "stats.h"

using namespace std;
 
int main (int argc, char *argv[]){

/****************************************************************************
setup
*****************************************************************************/

   //a generator is constructed
   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("./prng/Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();
   
   //seeds are taken in
   ifstream input("./prng/seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

/****************************************************************************
01.1.1 estimate <r> and its uncertainty, with blocks (plotted later in python) 
*****************************************************************************/
   int M = pow(10,6);  //number of throws
   int N = pow(10,3)/2;  //number of blocks
   int L = int(M/N);   //number of throws in each block
   vector<double> A_i; //each element is the result of a block

   auto id = [] (double x) { return x; }; //trivial function (id=identity)
   double inf = 0.; // integration boundaries
   double sup = 1.;

   // carry out N blocks of integrations
   cout << "01.1.1: estimating <r>..." << endl;
   for (int i=0; i<N; i++) {
      double integral = 0;                   
         for (int j=0; j<L; j++ ) { // which just means calculating an integral with L points
               double x = rnd.Rannyu(inf,sup);
               integral = integral*static_cast<double>(j)/static_cast<double>(j+1) + id(x)/static_cast<double>(j+1); //method of the average
         }
      integral = (sup-inf)*integral;
      A_i.push_back(integral);
      //cout << "Result A_"<<i+1<<": " << integral <<endl;
   }
   // let us build a progression of increasing blocks
   vector<double> averages_r;
   vector<double> sigma_r;
   for (int i=2; i<=N; i++) {
      vector<double> truncated_A_i(A_i.begin(), A_i.begin()+i); // the vector is sliced up
      //cout<<"Using "<<i<<" blocks"<<endl;                       // to account for the first i blocks only
      averages_r.push_back(Average(truncated_A_i));
      sigma_r.push_back(sqrt(Variance(truncated_A_i)/static_cast<double>(i-1)));
   }
   Print( averages_r , "averages.dat" );
   Print( sigma_r , "errors.dat" );

   averages_r.clear();
   sigma_r.clear();

/****************************************************************************
01.1.2 estimate the variance sigma^2 and its uncertainty (plotted later in python) 
*****************************************************************************/
   vector<double> B_i; //each element is the result of a block
   auto sq_err = [] (double x, double true_x) { return pow(x-true_x,2); }; //squared error

   // carry out N blocks of integrations (same inf & sup as above)
   cout << "01.1.2: estimating variance..." << endl;
   for (int i=0; i<N; i++) {
      double integral = 0;                   
         for (int j=0; j<L; j++ ) { // which just means calculating an integral with L points
               double x = rnd.Rannyu(inf,sup);
               integral = integral*static_cast<double>(j)/static_cast<double>(j+1) + sq_err(x,0.5)/static_cast<double>(j+1); //method of the average
         }
      integral = (sup-inf)*integral;
      B_i.push_back(integral);
   }
   // let us build a progression of increasing blocks
   vector<double> variance_averages_r;
   vector<double> variance_sigma_r;
   for (int i=2; i<=N; i++) {
      vector<double> truncated_B_i(B_i.begin(), B_i.begin()+i); // the vector is sliced up
      //cout<<"Using "<<i<<" blocks"<<endl;                       // to account for the first i blocks only
      variance_averages_r.push_back(Average(truncated_B_i));
      variance_sigma_r.push_back(sqrt(Variance(truncated_B_i)/static_cast<double>(i-1)));
   }
   Print( variance_averages_r , "variance_averages.dat" );
   Print( variance_sigma_r , "variance_errors.dat" );

   variance_averages_r.clear();
   variance_sigma_r.clear();

/****************************************************************************
01.1.3 plot chi^2 (plotted later in python) 
*****************************************************************************/
   int M_blocks=100;
   vector<double> chi_2(M_blocks);
   int N_generations=pow(10,4);

   //calculate the chi^2 for increasingly thin blocks (sub-intervals)
   cout << "01.1.3: calculating chi^2..." << endl;
   for (int number_of_blocks=1; number_of_blocks<=M_blocks; number_of_blocks++) {
      vector<double> frequencies(M_blocks);//bin populations
      for (int i=0; i<N_generations; i++){ //draw N_generations random numbers
         double x = rnd.Rannyu(0.,1.);
         for (int j=0; j<M_blocks; j++) { //this cycle checks whether the random number falls within [j/n,(j+1)/n)
            if ( x>=(static_cast<double>(j)/static_cast<double>(M_blocks)) && x<(static_cast<double>((j+1))/static_cast<double>(M_blocks)) ) {
               frequencies[j]+=1;//increases the population of the corresponding bin
               break;//stops the for cycle, since x can fall in one bin only
            }
         }
      }
      chi_2[number_of_blocks - 1] = 0; // Corrected index! vectors are filled from position 0
      for (int k = 0; k < M_blocks; k++)
        chi_2[number_of_blocks - 1] += pow(frequencies[k] - static_cast<double>(N_generations) / static_cast<double>(M_blocks), 2) / (static_cast<double>(N_generations) / static_cast<double>(M_blocks));
      // cout << chi_2[number_of_blocks - 1] << endl; // Corrected index!
   }
   Print( chi_2 , "chi_sq.dat" );
   chi_2.clear();

/****************************************************************************
01.2.1 create distributions 
*****************************************************************************/
   // exponential distr: parameters and vector for storage
   double lambda = 1.;
   vector<double> expo_distr;
   // lorentzian distr: parameters and vector for storage
   double gamma = 1.;
   double mu = 0.;
   vector<double> cau_lor_distr;
   //number of draws
   int K = pow(10,4);
   
   cout << "01.2.1: producting distrubutions..." << endl;
   for (int i=0; i<K; i++) {
      expo_distr.push_back(rnd.Expo(lambda));
      cau_lor_distr.push_back(rnd.Lorentz(mu,gamma));
   }

   Print( expo_distr , "exponential.dat" );
   Print( cau_lor_distr , "cauchy_lorentz.dat" );

   expo_distr.clear();
   cau_lor_distr.clear();

/****************************************************************************
01.2.2 play with dice 
*****************************************************************************/
   int W = 1; // then 2, 10, 100
   int n_realizations = pow(10,4); // number of dice extractions
   vector<double> S_N; 
   std::string N_title = std::to_string(W);
   std::string filename_str = "uniform_dice_" + N_title + ".dat";
   char* filename = &filename_str[0];

   cout << "01.2.2: casting dice..." << endl;
   while (W <= 100) {
   // cout<<"Dice extraction with N = "<<W<<endl;
   // set the title
   N_title = std::to_string(W);
   filename_str = "uniform_dice_" + N_title + ".dat";
   // uniform ===============================================
   for (int i=0; i<n_realizations; i++) {
      double S_N_i = 0;
      for (int l=0; l<W; l++) {
         S_N_i = S_N_i*static_cast<double>(l)/static_cast<double>(l+1) + rnd.Rannyu()/static_cast<double>(l+1);
      }
      S_N.push_back(S_N_i);
   }
   filename = &filename_str[0];
   Print(S_N, filename);
   S_N.clear();
   // exponential ===============================================
   // set the title
   N_title = std::to_string(W);
   filename_str = "exponential_dice_" + N_title + ".dat";
   //
   for (int i=0; i<n_realizations; i++) {
      double S_N_i = 0;
      for (int l=0; l<W; l++) {
         S_N_i = S_N_i*static_cast<double>(l)/static_cast<double>(l+1) + rnd.Expo(1.)/static_cast<double>(l+1);
      }
      S_N.push_back(S_N_i);
   }
   filename = &filename_str[0];
   Print(S_N, filename);
   S_N.clear();
   // lorenztian ===============================================
   // set the title
   N_title = std::to_string(W);
   filename_str = "lorentzian_dice_" + N_title + ".dat";
   //
   for (int i=0; i<n_realizations; i++) {
      double S_N_i = 0;
      for (int l=0; l<W; l++) {
         S_N_i = S_N_i*static_cast<double>(l)/static_cast<double>(l+1) + rnd.Lorentz(0.,1.)/static_cast<double>(l+1);
      }
      S_N.push_back(S_N_i);
   }
   filename = &filename_str[0];
   Print(S_N,filename);
   S_N.clear();
   // change W
   if (W == 1) {
      W = 2;
   } else if (W == 2) {
      W = 10;
   } else if (W == 10) {
      W = 100;
   } else if (W == 100) {
      W = 101;
   }
}
/****************************************************************************
01.3 simulate Buffon
*****************************************************************************/
   // parameters of the experiment
   double needle_L = 1.; // arbitrary
   double grating_d = 2.*needle_L; // spacing>needle_L, arbitrary 
   int N_thr_tot = pow(10,9); // Number of experiments (number of needles)
   int N_blk = pow(10,3); // Number of blocks
   int N_thr_blk = N_thr_tot/N_blk; // Number of throws in each experiment
   vector<double> pi_estimates; // vector for storage

   // carry out N blocks of experiments
   cout << "01.3: simulating Buffon..." << endl;
   for (int i=0; i<N_blk; i++) { cout<<i<<endl;
      int N_hit=0;                
         for (int j=0; j<N_thr_blk; j++ ) { // which just means doing the experiment with N_thr_blk blocks
            // == BUFFON'S EXPERIMENT ==
            // 1. extract the position of the centre of the needle
            double c = rnd.Rannyu(-grating_d*0.5,grating_d*0.5); 
            // 2. extract the needle orientation
            double x = 0.;
            double y = 0.;
            do {
               x = rnd.Rannyu();
               y = rnd.Rannyu();
            } while (x*x+y*y>1.); // only if inside of a circle sector
            double theta = atan(y/x); 
            // 3. check if the needle hits the grating
            if ( abs(0.5*needle_L*sin(theta))-abs(c)>0 ) {
               N_hit += 1;
            }
            // =========================
         }
      double approx_pi = (2.*needle_L/grating_d)*static_cast<double>(N_thr_blk)/static_cast<double>(N_hit);
      pi_estimates.push_back(approx_pi);
   } 
   // let us build a progression of increasing blocks
   vector<double> averages_PI;
   vector<double> sigma_PI;
   for (int i=2; i<=N_blk; i++) {
      vector<double> truncated_PI_i(pi_estimates.begin(), pi_estimates.begin()+i); // the vector is sliced up
      cout<<"Using "<<i<<" blocks"<<endl;                       // to account for the first i blocks only
      averages_PI.push_back(Average(truncated_PI_i));
      sigma_r.push_back(sqrt(Variance(truncated_PI_i)/static_cast<double>(i-1)));
   }

   Print( averages_PI , "averages_PI.dat" );
   Print( sigma_r , "errors_PI.dat" );
   averages_PI.clear();
   sigma_PI.clear();

/****************************************************************************
closing
*****************************************************************************/
   cout << "END" << endl;
   rnd.SaveSeed();
   return 0;
}
/****************************************************************************
*****************************************************************************/