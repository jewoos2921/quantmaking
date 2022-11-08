//
// Created by jewoo on 2022-11-06.
//

# include <python3.12/Python.h>

#include <cstdio>
#include <cstdlib>
#include <complex>

#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>
#include <numpy/npy_3kcompat.h>

typedef std::complex<double> cmplxd;
typedef std::complex<float> cmplxf;
static const int nbits = 22;
static cmplxd *psi;


void apply_single(cmplxd *psi, cmplxd gate[2][2],
                  int nbits, int qubit) {

    int q2 = 1 << qubit;
    for (int g = 0; g < 1 << nbits; g += (1 << (qubit + 1))) {
        for (int i = g; i < g + q2; ++i) {
            cmplxd t1 = gate[0][0] * psi[i] + gate[0][1] * psi[i + q2];
            cmplxd t2 = gate[1][0] * psi[i] + gate[1][1] * psi[i + q2];
            psi[i] = t1;
            psi[i + q2] = t2;
        }
    }
}

// optimized Gates.
void apply_single_opt(cmplxd *psi,
                      int nbits, int qubit) {

    int q2 = 1 << qubit;
    for (int g = 0; g < 1 << nbits; g += (1 << (qubit + 1))) {
        for (int i = g; i < g + q2; ++i) {
            cmplxd t1 = psi[i + q2];
            cmplxd t2 = psi[i];
            psi[i] = t1;
            psi[i + q2] = t2;
        }
    }
}