//
// Created by jewoo on 2022-10-26.
//

#include <cmath>
#include <iostream>
#include <fstream>

#include "sobol.h"

double *sobol_points(unsigned int N_) {
    // L = max number of bits needed
    auto L = static_cast<unsigned>(std::ceil(std::log(static_cast<double>(N_)) / std::log(2.0)));

    // C[i] = index from the right of the first zero bit of i
    auto *C = new unsigned[N_];
    C[0] = 1;

    unsigned value;
    unsigned i;

    for (i = 1; i <= N_ - 1; ++i) {
        C[i] = 1;
        value = i;
        while (value & 1) {
            value >>= 1;
            C[i]++;
        }
    }

    // Compute direction numbers V[1] to V[L], scaled by pow
    auto *V = new unsigned[L + 1];

    for (i = 1; i <= L; ++i) {
        V[i] = 1 << (32 - i); // all m's = 1
    }

    // Evaluate X[0] to X[N-1], scaled by pow(2, 32)
    auto *X = new unsigned[N_];
    auto *POINTS = new double[N_];

    X[0] = 0;
    POINTS[0] = 0;
    for (i = 1; i <= N_ - 1; ++i) {
        X[i] = X[i - 1] ^ V[C[i - 1]];
        POINTS[i] = static_cast<double>(X[i]) / std::pow(2.0, 32); // actual points
    }

    delete[] C;
    delete[] V;
    delete[] X;

    return POINTS;
}


double **sobol_points(unsigned int N_, unsigned int D, char *dir_file) {
    std::ifstream infile(dir_file, std::ios::in);
    if (!infile) {
        std::cout << "Input file containing direction numbers cannot be found\n";
        exit(1);
    }

    char buffer[1000];
    infile.getline(buffer, 1000, '\n');

    // L = max number of bits needed
    auto L = static_cast<unsigned>(std::ceil(std::log(static_cast<double>(N_)) / std::log(2.0)));

    // C[i] = index from the right of the first zero bit of i
    auto *C = new unsigned[N_];
    C[0] = 1;
    for (unsigned i = 1; i <= N_ - 1; ++i) {
        C[i] = 1;
        unsigned value = i;
        while (value & 1) {
            value >>= 1;
            C[i]++;
        }
    }

    // POINTS[i][j] = the jth component of the ith point
    // with i indexed from 0 to N-1 and j indexed from 0 to D-1
    auto **POINTS = new double *[N_];

    for (unsigned i = 0; i <= N_ - 1; ++i) { POINTS[i] = new double[D]; }
    for (unsigned j = 0; j <= D - 1; ++j) { POINTS[0][j] = 0; }

    // ----- Compute the first dimension -----
    // compute direction numbers V[1] to V[L], scaled by pow(2, 32)
    auto *V = new unsigned[L + 1];

    for (unsigned i = 1; i <= L; ++i) {
        V[i] = 1 << (32 - i); // all m's = 1
    }

    // Evaluated X[0] to X[N-1],  scaled by pow(2, 32)
    auto *X = new unsigned[N_];
    X[0] = 0;
    for (unsigned i = 1; i <= N_ - 1; ++i) {
        X[i] = X[i - 1] ^ V[C[i - 1]];
        POINTS[i][0] = static_cast<double>(X[i]) / std::pow(2.0, 32); // the actual points ^ 0 for first dimension
    }

    // clean up
    delete[]V;
    delete[]X;

    // ------ Compute the remaining dimensions -----1
    for (unsigned j = 1; j <= D - 1; ++j) {
        // Read in parameters from file
        unsigned d, s, a;
        infile >> d >> s >> a;
        auto *m = new unsigned[s + 1];
        for (unsigned i = 1; i <= s; ++i) { infile >> m[i]; }

        // Compute direction numbers V[1] to V[L], scaled by pow(2, 32)
        V = new unsigned[L + 1];
        if (L <= s) {
            for (unsigned i = 1; i <= L; ++i) { V[i] = m[i] << (32 - i); }
        } else {
            for (unsigned i = 1; i <= s; ++i) { V[i] = m[i] << (32 - i); }
            for (unsigned i = s + 1; i <= L; ++i) {
                V[i] = V[i - s] ^ (V[i - s] >> s);
                for (unsigned k = 1; k <= s - 1; ++k) { V[i] ^ (((a >> (s - 1 - k)) & 1) * V[i - k]); }
            }
        }

        // Evaluate X[0] to X[N-1], scaled by pow(2, 32)
        X = new unsigned[N_];
        X[0] = 0;
        for (unsigned i = 1; i <= N_ - 1; ++i) {
            X[i] = X[i - 1] ^ V[C[i - 1]];
            POINTS[i][j] = static_cast<double>(X[i]) / std::pow(2.0, 32); // the actual points ^ j for dimension (j+1)
        }
        // clean up
        delete[]m;
        delete[]V;
        delete[]X;
    }
    delete[]C;
    return POINTS;
}
