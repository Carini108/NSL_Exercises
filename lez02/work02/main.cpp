/****************************************************************
*****************************************************************
LSN_Exercises_02
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
    02.1.1 integral with uniform sampling
    *****************************************************************************/
    int M = pow(10, 5);  // number of throws
    int N = pow(10, 2);         // number of blocks
    int L = int(M / N);   // number of throws in each block

    vector<double> uniform_I_i; // each element is the result of a block

    // integration boundaries
    double inf = 0.;
    double sup = 1.;
    // integrand and parameters
    auto cosinusoid = [](double x, double A, double omega, double phi) { return A * cos(omega * x + phi); };
    double A = M_PI * 0.5;
    double omega = M_PI * 0.5;
    double phi = 0.;

    // carry out N blocks of integrations
    cout << "02.1.1: estimating <I> with uniform sampling..." << endl;
    for (int i = 0; i < N; i++) {
        double integral = 0;
        for (int j = 0; j < L; j++) {
            double x = rnd.Rannyu(inf, sup);
            integral = integral * static_cast<double>(j) / static_cast<double>(j + 1) +
                       cosinusoid(x, A, omega, phi) / static_cast<double>(j + 1); // method of the average
        }
        integral = (sup - inf) * integral;
        uniform_I_i.push_back(integral);
    }
    // let us build a progression of increasing blocks
    vector<double> uniform_I_ave;
    vector<double> uniform_I_sigma;
    for (int i = 2; i <= N; i++) {
        vector<double> truncated_uniform_I_i(uniform_I_i.begin(), uniform_I_i.begin() + i); // the vector is sliced up
        //cout<<"Using "<<i<<" blocks"<<endl;                       // to account for the first i blocks only
        uniform_I_ave.push_back(Average(truncated_uniform_I_i));
        uniform_I_sigma.push_back(sqrt(Variance(truncated_uniform_I_i) / static_cast<double>(i - 1)));
    }
    Print(uniform_I_ave, "I_uniform_averages.dat");
    Print(uniform_I_sigma, "I_uniform_errors.dat");

    uniform_I_ave.clear();
    uniform_I_sigma.clear();

    /****************************************************************************
    02.1.2 integral with importance sampling
    *****************************************************************************/
    vector<double> importance_I_i; // each element is the result of a block

    // carry out N blocks of integrations
    cout << "02.1.2: estimating <I> with importance sampling..." << endl;
    for (int i = 0; i < N; i++) {
        double integral = 0.;
        for (int j = 0; j < L; j++) {
            double y = rnd.Rannyu(inf, sup);
            double x = 1. - sqrt(1. - y); // importance sampling
            integral = integral * static_cast<double>(j) / static_cast<double>(j + 1) +
                       cosinusoid(y, A, omega, phi) / static_cast<double>(j + 1); // method of the average
        }
        integral = (sup - inf) * integral;
        importance_I_i.push_back(integral);
    }
    // let us build a progression of increasing blocks
    vector<double> importance_I_ave;
    vector<double> importance_I_sigma;
    for (int i = 2; i <= N; i++) {
        vector<double> truncated_importance_I_i(importance_I_i.begin(), importance_I_i.begin() + i); // the vector is sliced up
        //cout<<"Using "<<i<<" blocks"<<endl;                       // to account for the first i blocks only
        importance_I_ave.push_back(Average(truncated_importance_I_i));
        importance_I_sigma.push_back(sqrt(Variance(truncated_importance_I_i) / static_cast<double>(i - 1)));
    }
    Print(importance_I_ave, "I_importance_averages.dat");
    Print(importance_I_sigma, "I_importance_errors.dat");

    importance_I_ave.clear();
    importance_I_sigma.clear();

    /****************************************************************************
    02.2.1 RW discrete lattice
    *****************************************************************************/
    int a = 1;             // lattice spacing
    int m_throws = pow(10, 4);
    int n_blocks = pow(10, 2);
    int steps = 100;
    // **************************************************************************

    cout << "02.2.1: discrete RW..." << endl;
    // this matrix stores, for each step i, a number n_blocks of the positions
    vector<vector<double>> r2_DISC;
    r2_DISC.resize(n_blocks, vector<double>(steps, 0.0)); // r^2_i for different blocks
    // fill the 'reduced' matrix
    for (int k = 0; k < n_blocks; k++) {
        // this matrix helps perform the calculation for a block
        // after a block, it is 'reduced' in a column of r2_DISC
        vector<vector<double>> matrix_DISC;
        matrix_DISC.resize(m_throws / n_blocks, vector<double>(steps, 0.0));
        // the matrix to be reduced is now filled:
        for (int j = 0; j < m_throws / n_blocks; j++) {
            // a walk is now performed ==============================================
            // ======================================================================
            double X_DISC = 0, Y_DISC = 0, Z_DISC = 0; // start at zero
            for (int i = 0; i < steps; i++) {
                // take a step
                double t = rnd.Rannyu();
                if (t < 1. / 6.) { X_DISC += 1; }
                else if (1. / 6. <= t && t < 2. / 6.) { X_DISC -= 1; }
                else if (2. / 6. <= t && t < 3. / 6.) { Y_DISC += 1; }
                else if (3. / 6. <= t && t < 4. / 6.) { Y_DISC -= 1; }
                else if (4. / 6. <= t && t < 5. / 6.) { Z_DISC += 1; }
                else { Z_DISC -= 1; }
                // step taken ======================================================================
                // SAVE R^2_i for the j-th execution
                matrix_DISC[j][i] = X_DISC * X_DISC + Y_DISC * Y_DISC + Z_DISC * Z_DISC;
            }
            // ======================================================================
            // ======================================================================
        }
        // the k-th column of the matrix of the vectors of the reduced matrices is filled:
        for (int i = 0; i < steps; i++) { // each row is indexed by i
            for (int j = 0; j < m_throws / n_blocks; j++) { // average on the columns of the matrix to be reduced
                r2_DISC[k][i] += matrix_DISC[j][i];
            }
            r2_DISC[k][i] = r2_DISC[k][i] / static_cast<double>(m_throws / n_blocks);
        }
    }

    // now the data to be plotted are evaluated
    vector<double> DISC_r_2_AVG(steps, 0.0); // the last matrix is reduced into a column vector
    vector<double> DISC_r_2_ERR(steps, 0.0);

    for (int i = 0; i < steps; i++) {
        vector<double> auxiliary_vector(n_blocks, 0.0);
        for (int j = 1; j < n_blocks; j++) {
            auxiliary_vector[j] = r2_DISC[j][i];
        }
        DISC_r_2_AVG[i] = Average(auxiliary_vector);
        DISC_r_2_ERR[i] = sqrt(Variance(auxiliary_vector)) / static_cast<double>(n_blocks - 1);
    }

    // we are interested in the square roots
    vector<double> DISC_sqrt_of_r_2_AVG(steps, 0.0);
    for (int i = 0; i < steps; i++) {
        DISC_sqrt_of_r_2_AVG[i] = (sqrt(DISC_r_2_AVG[i]));
    }
    // propagate errors
    vector<double> DISC_sqrt_of_r_2_ERR(steps, 0.0);
    for (int i = 0; i < steps; i++) {
        DISC_sqrt_of_r_2_ERR[i] = (DISC_r_2_ERR[i] * 0.5 * (1. / DISC_sqrt_of_r_2_AVG[i]));
    }
    // save in a file
    Print(DISC_sqrt_of_r_2_AVG, "discrete_AVG.dat");
    Print(DISC_sqrt_of_r_2_ERR, "discrete_ERR.dat");

    /****************************************************************************
    02.2.2 RW continuum
    *****************************************************************************/
    cout << "02.2.2: RW in the continuum..." << endl;
    // this matrix stores, for each step i, a number n_blocks of the positions
    vector<vector<double>> r2_CONT;
    r2_CONT.resize(n_blocks, vector<double>(steps, 0.0)); // r^2_i for different blocks
    // fill the 'reduced' matrix
    for (int k = 0; k < n_blocks; k++) {
        // this matrix helps perform the calculation for a block
        // after a block, it is 'reduced' in a column of r2_CONT
        vector<vector<double>> matrix_CONT;
        matrix_CONT.resize(m_throws / n_blocks, vector<double>(steps, 0.0));
        // do many (m_throws/n_blocks) blocks
        for (int j = 0; j < m_throws / n_blocks; j++) {
            // a walk is now performed ==============================================
            // ======================================================================
            double X_CONT = 0., Y_CONT = 0., Z_CONT = 0.; // start the walk from the origin
            for (int i = 0; i < steps; i++) {
                // take a step: first extract the angle in 3d =====================================
                double x = 0., y = 0., z = 0.;
                do {
                    x = rnd.Rannyu(-1., 1.);
                    y = rnd.Rannyu(-1., 1.);
                    z = rnd.Rannyu(-1., 1.);
                } while (x * x + y * y + z * z > 1. || x * x + y * y + z * z == 0); // only if inside of the unit sphere
                double THETA = atan2(y, x);
                double PHI = acos(z / sqrt(x * x + y * y + z * z));
                // spherical to cartesian: displacement
                X_CONT += static_cast<double>(a) * sin(THETA) * cos(PHI);
                Y_CONT += static_cast<double>(a) * sin(THETA) * sin(PHI);
                Z_CONT += static_cast<double>(a) * cos(THETA);
                // step taken ======================================================================
                // SAVE R^2_i for the j-th execution
                matrix_CONT[j][i] = X_CONT * X_CONT + Y_CONT * Y_CONT + Z_CONT * Z_CONT;
            }
            // ======================================================================
            // ======================================================================
        }
        for (int i = 0; i < steps; i++) {
            for (int j = 0; j < m_throws / n_blocks; j++) {
                r2_CONT[k][i] += matrix_CONT[j][i];
            }
            r2_CONT[k][i] = r2_CONT[k][i] / static_cast<double>(m_throws / n_blocks);
        }
    }

    vector<double> r_2_AVG(steps, 0.0);
    vector<double> r_2_ERR(steps, 0.0);

    for (int i = 0; i < steps; i++) {
        vector<double> auxiliary_vector(n_blocks, 0.0);
        for (int j = 1; j < n_blocks; j++) {
            auxiliary_vector[j] = r2_CONT[j][i];
        }
        r_2_AVG[i] = Average(auxiliary_vector);
        r_2_ERR[i] = sqrt(Variance(auxiliary_vector)) / static_cast<double>(n_blocks - 1);
    }

    // we are interested in the square roots
    vector<double> sqrt_of_r_2_AVG(steps, 0.0);
    for (int i = 0; i < steps; i++) {
        sqrt_of_r_2_AVG[i] = (sqrt(r_2_AVG[i]));
    }
    // propagate errors
    vector<double> sqrt_of_r_2_ERR(steps, 0.0);
    for (int i = 0; i < steps; i++) {
        sqrt_of_r_2_ERR[i] = (r_2_ERR[i] * 0.5 * (1. / sqrt_of_r_2_AVG[i]));
    }
    //
    Print(sqrt_of_r_2_AVG, "continuum_AVG.dat");
    Print(sqrt_of_r_2_ERR, "continuum_ERR.dat");

    /****************************************************************************
    closing
    *****************************************************************************/
    cout << "END" << endl;
    rnd.SaveSeed();
    return 0;
}
/****************************************************************************
*****************************************************************************/