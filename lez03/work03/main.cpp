/****************************************************************
*****************************************************************
LSN_Exercises_03
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "./prng/random.h"
#include "stats.h"

using namespace std;

int main(int argc, char *argv[]) {
    /****************************************************************************
    setup
    *****************************************************************************/

    //a generator is constructed
    Random rnd;
    int seed[4];
    int p1, p2;
    ifstream Primes("./prng/Primes");

    if (Primes.is_open()) {
        Primes >> p1 >> p2;
    } else {
        cerr << "PROBLEM: Unable to open Primes" << endl;
    }
    Primes.close();

    //seeds are taken in
    ifstream input("./prng/seed.in");
    string property;

    if (input.is_open()) {
        while (!input.eof()) {
            input >> property;
            if (property == "RANDOMSEED") {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed, p1, p2);
            }
        }
        input.close();
    } else {
        cerr << "PROBLEM: Unable to open seed.in" << endl;
    }

    /****************************************************************************
    SIMULATION PARAMETERS and S(t) to be evolved
    *****************************************************************************/
    double S_0 = 100., T=1., K=100., r=0.1, sigma=0.25;
    double S_t = S_0;
    int M = pow(10,5); // # of throws
    int N = pow(10,2); // # of blocks
    int L = M/N; // # of throws per block

    /****************************************************************************
    03.1.1 direct CALL
    *****************************************************************************/
    cout<<"03.1.1 direct CALL..."<<endl;
    // fix the time step
    double delta_t = T;
    vector<double> A_J(L); //auxiliary vector
    // do N blocks
    for (int j=0; j<N; j++) {
        // do one block
        vector<double> A_blk(L); //auxiliary vector
        for (int i=0; i<L; i++) {
            S_t = S_0; // initialize the starting point
            for (double time=0; time<T; time+=delta_t){ // take all steps to get to final time T (here, just one step)
                S_t = S_t*exp( (r-sigma*sigma*0.5)*delta_t + sigma*rnd.Gauss(0.,1.)*sqrt(delta_t) ); // take a step of GBM
            }
            A_blk[i] = max( 0. , exp(-r*T)*(S_t-K) );// save the final value (shows whether it's a good idea to use the option or not)
        }
        A_J[j] = Average(A_blk);
        A_blk.clear();
    }
    // let us build a progression of increasing blocks
    vector<double> ave_CALL_direct;
    vector<double> sigma_CALL_direct;
    for (int i=2; i<=N; i++) {
        vector<double> truncated_CALL_direct(A_J.begin(), A_J.begin()+i); // the vector is sliced up
        ave_CALL_direct.push_back(Average(truncated_CALL_direct));
        sigma_CALL_direct.push_back(sqrt(Variance(truncated_CALL_direct)/static_cast<double>(i-1)));
    }
    // save in a file
    Print( ave_CALL_direct , "ave_CALL_direct.dat" );
    Print( sigma_CALL_direct , "sigma_CALL_direct.dat" );
    ave_CALL_direct.clear();
    sigma_CALL_direct.clear();
    A_J.clear();
    
    /****************************************************************************
    03.1.1 direct PUT
    *****************************************************************************/
    cout<<"03.1.1 direct PUT..."<<endl;
    // fix the time step
    delta_t = T;
    // do N blocks
    for (int j=0; j<N; j++) {
        // do one block
        vector<double> A_blk(L); //auxiliary vector
        for (int i=0; i<L; i++) {
            S_t = S_0; // initialize the starting point
            for (double time=0; time<T; time+=delta_t){ // take all steps to get to final time T (here, just one step)
                S_t = S_t*exp( (r-sigma*sigma*0.5)*delta_t + sigma*rnd.Gauss(0.,1.)*sqrt(delta_t) ); // take a step of GBM
            }
            A_blk[i] = max( 0. , exp(-r*T)*(K-S_t) );// save the final value (shows whether it's a good idea to use the option or not)
        }
        A_J[j] = Average(A_blk);
        A_blk.clear();
    }
    // let us build a progression of increasing blocks
    vector<double> ave_PUT_direct;
    vector<double> sigma_PUT_direct;
    for (int i=2; i<=N; i++) {
        vector<double> truncated_PUT_direct(A_J.begin(), A_J.begin()+i); // the vector is sliced up
        ave_PUT_direct.push_back(Average(truncated_PUT_direct));
        sigma_PUT_direct.push_back(sqrt(Variance(truncated_PUT_direct)/static_cast<double>(i-1)));
    }
    // save in a file
    Print( ave_PUT_direct , "ave_PUT_direct.dat" );
    Print( sigma_PUT_direct , "sigma_PUT_direct.dat" );
    ave_PUT_direct.clear();
    sigma_PUT_direct.clear();
    A_J.clear();
    
    /****************************************************************************
    03.1.2 discrete CALL
    *****************************************************************************/
    cout<<"03.1.2 discrete CALL..."<<endl;
    // fix the time step
    delta_t = T/100.;
    // do N blocks
    for (int j=0; j<N; j++) {
        // do one block
        vector<double> A_blk(L); //auxiliary vector
        for (int i=0; i<L; i++) {
            S_t = S_0; // initialize the starting point
            for (double time=0; time<T; time+=delta_t){ // take all steps to get to final time T (here, just one step)
                S_t = S_t*exp( (r-sigma*sigma*0.5)*delta_t + sigma*rnd.Gauss(0.,1.)*sqrt(delta_t) ); // take a step of GBM
            }
            A_blk[i] = max( 0. , exp(-r*T)*(S_t-K) );// save the final value (shows whether it's a good idea to use the option or not)
        }
        A_J[j] = Average(A_blk);
        A_blk.clear();
    }
    // let us build a progression of increasing blocks
    vector<double> ave_CALL_discrete;
    vector<double> sigma_CALL_discrete;
    for (int i=2; i<=N; i++) {
        vector<double> truncated_CALL_direct(A_J.begin(), A_J.begin()+i); // the vector is sliced up
        ave_CALL_discrete.push_back(Average(truncated_CALL_direct));
        sigma_CALL_discrete.push_back(sqrt(Variance(truncated_CALL_direct)/static_cast<double>(i-1)));
    }
    // save in a file
    Print( ave_CALL_discrete , "ave_CALL_discrete.dat" );
    Print( sigma_CALL_discrete , "sigma_CALL_discrete.dat" );
    ave_CALL_discrete.clear();
    sigma_CALL_discrete.clear();
    A_J.clear();


    /****************************************************************************
    03.1.2 discrete PUT
    *****************************************************************************/
    cout<<"03.1.2 discrete PUT..."<<endl;
    // fix the time step
    delta_t = T/100.;
    // do N blocks
    for (int j=0; j<N; j++) {
        // do one block
        vector<double> A_blk(L); //auxiliary vector
        for (int i=0; i<L; i++) {
            S_t = S_0; // initialize the starting point
            for (double time=0; time<T; time+=delta_t){ // take all steps to get to final time T (here, just one step)
                S_t = S_t*exp( (r-sigma*sigma*0.5)*delta_t + sigma*rnd.Gauss(0.,1.)*sqrt(delta_t) ); // take a step of GBM
            }
            A_blk[i] = max( 0. , exp(-r*T)*(K-S_t) );// save the final value (shows whether it's a good idea to use the option or not)
        }
        A_J[j] = Average(A_blk);
        A_blk.clear();
    }
    // let us build a progression of increasing blocks
    vector<double> ave_PUT_discrete;
    vector<double> sigma_PUT_discrete;
    for (int i=2; i<=N; i++) {
        vector<double> truncated_PUT_direct(A_J.begin(), A_J.begin()+i); // the vector is sliced up
        ave_PUT_discrete.push_back(Average(truncated_PUT_direct));
        sigma_PUT_discrete.push_back(sqrt(Variance(truncated_PUT_direct)/static_cast<double>(i-1)));
    }
    // save in a file
    Print( ave_PUT_discrete , "ave_PUT_discrete.dat" );
    Print( sigma_PUT_discrete , "sigma_PUT_discrete.dat" );
    ave_PUT_discrete.clear();
    sigma_PUT_discrete.clear();
    A_J.clear();

    /****************************************************************************
    closing
    *****************************************************************************/
    cout << "END" << endl;
    rnd.SaveSeed();
    return 0;
}
/****************************************************************************
*****************************************************************************/