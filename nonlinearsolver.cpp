//
// Created by jewoo on 2022-10-15.
//

#include <iostream>
#include <cmath>
#include "nonlinearsolver.h"
#include "mathlib.h"

void gauss_newton_parameter_solver(int n_ref, double *ref_x, double *ref_y, int n_parameter, double *parameter,
                                   double (*object_func)(int, double *, double)) {

    // 비선형 파라미터 추정
    int iter, i, k, j;
    double *guess, *unknown, *known, *perr;
    double **s_matrix;
    double error, old_error, tmp_error, tolerance, h;
    double fx_value;

    iter = 0;
    tolerance = 1.0e-15;
    h = 0.00001;

    perr = new double[n_ref];
    guess = new double[n_parameter];
    unknown = new double[n_parameter];
    known = new double[n_parameter];
    s_matrix = new double *[n_parameter];

    for (i = 0; i < n_parameter; ++i) {
        s_matrix[i] = new double[n_parameter];
    }

    for (i = 0; i < n_parameter; ++i) {
        guess[i] = parameter[i];
    }

    error = 0.0;
    for (i = 0; i < n_ref; ++i) {
        fx_value = object_func(n_parameter, guess, ref_x[i]);
        perr[i] = fx_value - ref_y[i];
        error += perr[i] * perr[i];
    }

    old_error = error;

    double df1, df2, faph, famh, fbph, fbmh;

    iter = 0;
    tmp_error = 1.0e15;

    do {
        for (i = 0; i < n_parameter; ++i) {
            for (j = 0; j < n_parameter; ++j) {
                s_matrix[i][j] = 0.0;
            }
            known[i] = unknown[i] = 0.0;
        }

        for (i = 0; i < n_parameter; ++i) {
            for (j = 1; j < n_parameter; ++j) {
                for (k = 0; k < n_ref; ++k) {
                    guess[i] += h;
                    faph = object_func(n_parameter, guess, ref_x[k]);
                    guess[i] -= 2.0 * h;
                    famh = object_func(n_parameter, guess, ref_x[k]);
                    guess[i] += h;

                    guess[j] += h;
                    fbph = object_func(n_parameter, guess, ref_x[k]);
                    guess[j] -= 2.0 * h;
                    fbmh = object_func(n_parameter, guess, ref_x[k]);
                    guess[j] += h;

                    df1 = (faph - famh) / (2.0 * h);
                    df2 = (fbph - fbmh) / (2.0 * h);
                    s_matrix[i][j] += df1 * df2;
                }
            }
        }

        for (i = 1; i < n_parameter; ++i) {
            for (j = 0; j < i; ++j) {
                s_matrix[i][j] = s_matrix[j][i];
            }
        }

        for (i = 0; i < n_parameter; ++i) {
            for (k = 0; k < n_ref; ++k) {
                guess[i] += h;
                faph = object_func(n_parameter, guess, ref_x[k]);
                guess[i] -= 2.0 * h;
                famh = object_func(n_parameter, guess, ref_x[k]);
                guess[i] += h;

                df1 = (faph - famh) / (2.0 * h);
                known[i] -= df1 * perr[k];
            }
        }

        gaussian_elimination(s_matrix, known, unknown, n_parameter);

        for (i = 0; i < n_parameter; ++i) {
            guess[i] += unknown[i];
        }

        tmp_error = old_error;
        error = 0.0;

        for (i = 0; i < n_ref; ++i) {
            fx_value = object_func(n_parameter, guess, ref_x[i]);
            perr[i] = fx_value - ref_y[i];
            error += perr[i] * perr[i];
        }

        old_error = error;
        std::cout << iter << "번쨰" << std::endl;
        for (i = 0; i < n_parameter; ++i) {
            std::cout << "p[" << i << "] : " << guess[i] << std::endl;
        }
        iter++;
    } while (error > tolerance && std::fabs(tmp_error - error) > tolerance && iter < 1000);

    for (i = 0; i < n_parameter; ++i) {
        parameter[i] = guess[i];
    }
    std::cout << iter << "번쨰" << std::endl;
    for (i = 0; i < n_parameter; ++i) {
        delete s_matrix[i];
    }

    delete[] unknown;
    delete[]known;
    delete[] guess;
    delete[]s_matrix;

}

double fx_function(int n_parameter, double *parameter, double x) {
    double *p;
    int i;
    p = new double[n_parameter];
    for (i = 0; i < n_parameter; ++i) {
        p[i] = parameter[i];
    }

    return p[0] * x + p[1] * x * x + p[2] * std::exp(p[3] * x);
}
