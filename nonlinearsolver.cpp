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

void levenberg_marquardt_parameter_solver(double S, double cc, int n_x,
                                          double *v_x, int n_y, double *v_y,
                                          double **m_ref, int n_parameter,
                                          double *parameter,
                                          double object_func(int n_parameter, double *parameter, double x, double t)) {

    // 빈선형 파라미터 추정
    int i, j, n, m, iter;
    double *t_guess, *unknown, *known, *min_guess, *guess;
    double **s_matrix, **perr, **t_perr;
    double error, old_error, tmp_error, tolerance, h, min_error;
    double fx_value, x_value;

    iter = 0;
    tolerance = 1.0e-15;
    h = 0.00001;

    perr = new double *[n_x];
    t_perr = new double *[n_x];

    for (i = 0; i < n_x; ++i) {
        perr[i] = new double[n_y];
        t_perr[i] = new double[n_y];
    }

    guess = new double[n_parameter];
    t_guess = new double[n_parameter];
    min_guess = new double[n_parameter];
    unknown = new double[n_parameter];
    known = new double[n_parameter];
    s_matrix = new double *[n_parameter];

    for (i = 0; i < n_parameter; ++i) {
        s_matrix[i] = new double[n_parameter];
    }

    for (i = 0; i < n_parameter; ++i) {
        guess[i] = parameter[i];
    }

    double df1, df2, faph, famh, fbph, fbmh, lambda;
    lambda = 0.0001;

    error = 0.0;
    for (i = 0; i < n_x; ++i) {
        for (j = 0; j < n_y; ++j) {
            x_value = std::log(S * std::exp(cc * v_y[j]) / v_x[i]);
            fx_value = object_func(n_parameter, guess, x_value, v_y[j]);
            perr[i][j] = (fx_value - m_ref[i][j]);
            error += perr[i][j] * perr[i][j];
        }
    }
    old_error = error;
    iter = 0;
    tmp_error = 0.0;
    min_error = 1.0e10;

    do {
        for (i = 0; i < n_parameter; ++i) {
            for (j = 0; j < n_parameter; ++j) {
                s_matrix[i][j] = 0.0;
            }
            known[i] = unknown[i] = 0.0;
        }

        for (i = 0; i < n_parameter; ++i) {
            for (j = 0; j < n_parameter; ++j) {
                for (m = 0; m < n_x; ++m) {
                    for (n = 0; n < n_y; ++n) {
                        x_value = std::log(S * std::exp(cc * v_y[n]) / v_x[m]);
                        guess[i] += h;
                        faph = object_func(n_parameter, guess, x_value, v_y[n]);
                        guess[i] -= 2.0 * h;
                        famh = object_func(n_parameter, guess, x_value, v_y[n]);
                        guess[i] += h;

                        guess[j] += h;
                        fbph = object_func(n_parameter, guess, x_value, v_y[n]);
                        guess[j] -= 2.0 * h;
                        fbmh = object_func(n_parameter, guess, x_value, v_y[n]);
                        guess[j] += h;

                        df1 = (faph - famh) / (2.0 * h);
                        df2 = (fbph - fbmh) / (2.0 * h);
                        s_matrix[i][j] += df1 * df2;
                    }
                }
            }
        }

        for (i = 0; i < n_parameter; ++i) {
            for (j = 0; j < i; ++j) {
                s_matrix[i][j] = s_matrix[j][i];
            }
        }

        for (i = 0; i < n_parameter; ++i) {
            s_matrix[i][i] = s_matrix[i][i] * (1.0 + lambda);
            for (m = 0; m < n_x; ++m) {
                for (n = 0; n < n_y; ++n) {
                    x_value = std::log(S * std::exp(cc * v_y[n]) / v_x[n]);
                    guess[i] += h;
                    faph = object_func(n_parameter, guess, x_value, v_y[n]);
                    guess[i] -= 2.0 * h;
                    famh = object_func(n_parameter, guess, x_value, v_y[n]);
                    guess[i] += h;

                    df1 = (faph - famh) / (2.0 * h);
                    known[i] -= df1 * perr[m][n];
                }
            }
        }

        gaussian_elimination(s_matrix, known, unknown, n_parameter);

        for (i = 0; i < n_parameter; ++i) {
            t_guess[i] = guess[i] + unknown[i];
        }
        error = 0.0;
        for (i = 0; i < n_x; ++i) {
            for (j = 0; j < n_y; ++j) {
                x_value = std::log(S * std::exp(cc * v_y[j]) / v_x[i]);
                fx_value = object_func(n_parameter, t_guess, x_value, v_y[j]);
                t_perr[i][j] = (fx_value - m_ref[i][j]);
                error += t_perr[i][j] * t_perr[i][j];
            }
        }

        if (iter != 0 && error > old_error) {
            lambda *= 10.0;
        } else {
            lambda /= 10.0;
            for (i = 0; i < n_parameter; ++i) {
                guess[i] = t_guess[i];
            }
            for (i = 0; i < n_x; ++i) {
                for (j = 0; j < n_y; ++j) {
                    perr[i][j] = t_perr[i][j];
                }
            }
            tmp_error = old_error;
            old_error = error;
        }
        if (error < min_error) {
            for (i = 0; i < n_parameter; ++i) {
                min_guess[i] = guess[i];
            }
            min_error = error;
        }

        std::cout << iter << "번째" << std::endl;
        for (i = 0; i < n_parameter; ++i) {
            std::cout << "p[" << i << "] : " << guess[i] << std::endl;
        }

        std::cout << "error: " << min_error << std::endl;
        iter++;
    } while (error > tolerance && std::fabs(tmp_error - error) > tolerance && iter < 1000);
    for (i = 0; i < n_parameter; ++i) {
        parameter[i] = guess[i];
    }

    std::cout << iter << "번 실행 후 최적해 구함 " << std::endl;

    for (i = 0; i < n_x; ++i) {
        delete perr[i];
        delete t_perr[i];
    }

    delete[] perr;
    delete[] t_perr;

    for (i = 0; i < n_parameter; ++i) {
        delete s_matrix[i];
    }

    delete[] unknown;
    delete[] known;
    delete[] guess;
    delete[] t_guess;
    delete[] min_guess;
    delete[] s_matrix;
}

double implied_volatility_function(int n_parameter, double *parameter, double x, double t) {

    int i;
    double *p;
    double imvol;
    p = new double[n_parameter];
    for (i = 0; i < n_parameter; ++i) { p[i] = parameter[i]; }

    imvol = p[0] + (p[2] * std::exp(p[1] * t))
            + (p[3] * x) + (p[4] * x * x) + (p[5] * x * x * x) +
            (p[6] * x * x * x * x);
    delete[] p;
    return imvol;
}

void implied_volatility(double S, double cc, int n_x, double *v_x,
                        int n_y, double *v_y, double **mat, int n_parameter, double *parameter) {
    int i, j;
    double x_value;
    for (i = 0; i < n_x; ++i) {
        for (j = 0; j < n_y; ++j) {
            x_value = std::log(S * std::exp(cc * v_y[j]) / v_x[i]);
            mat[i][j] = implied_volatility_function(n_parameter,
                                                    parameter,
                                                    x_value, v_y[j]);
        }
    }
}

void nonlinear_parameter_solver(double S, double cc, int n_x, double *v_x, int n_y, double *v_y, double **mat,
                                int n_parameter, double *parameter,
                                double object_func(int n_parameter, double *parameter, double x, double t)) {

    // 비선형 파라미터 추정
    int i, j, m, n, iter;
    double *unknown, *known, *guess, *tguess, *minguess;
    double **smatrix, **perr, **tperr;
    double error, olderror, minerror, tmperror, tolerance, h;
    double fx_value, xvalue;

    iter = 0;
    tolerance = 1.0e-15;
    h = 0.001;

    perr = new double *[n_x];
    tperr = new double *[n_x];

    for (i = 0; i < n_x; ++i) {
        perr[i] = new double[n_y];
        tperr[i] = new double[n_y];
    }

    guess = new double[n_parameter];
    known = new double[n_parameter];
    unknown = new double[n_parameter];
    tguess = new double[n_parameter];
    minguess = new double[n_parameter];
    smatrix = new double *[n_parameter];
    for (i = 0; i < n_parameter; ++i) { smatrix[i] = new double[n_parameter]; }
    for (i = 0; i < n_parameter; ++i) { guess[i] = parameter[i]; }

    double df1, df2;
    double faph, famh, fbph, fbmh, lamda;
    lamda = 0.0001;
    error = 0.0;

    for (i = 0; i < n_x; ++i) {
        for (j = 0; j < n_y; ++j) {
            xvalue = std::log(S / v_x[i]) + (cc * v_y[j]); // std::log(S * std::exp(cc * v_y[j]) / v_x[i]);
            fx_value = object_func(n_parameter, guess, xvalue, v_y[j]);
            perr[i][j] = (fx_value - mat[i][j]) / mat[i][j];
            error += perr[i][j] * perr[i][j];
        }
    }

    olderror = error;
    iter = 0;
    tmperror = 0.0;
    minerror = 1.0e10;

    do {
        for (i = 0; i < n_parameter; ++i) {
            for (j = 0; j < n_parameter; ++j) { smatrix[i][j] = 0.0; }
            known[i] = unknown[i] = 0.0;
        }
        for (i = 0; i < n_parameter; ++i) {
            for (j = i; j < n_parameter; ++j) {
                for (m = 0; m < n_x; ++m) {
                    for (n = 0; n < n_y; ++n) {
                        xvalue = std::log(S / v_x[m]) + (cc * v_y[n]);
                        guess[i] += h;
                        faph = object_func(n_parameter, guess, xvalue, v_y[n]);
                        guess[i] -= 2.0 * h;
                        famh = object_func(n_parameter, guess, xvalue, v_y[n]);
                        guess[i] += h;

                        guess[i] += h;
                        fbph = object_func(n_parameter, guess, xvalue, v_y[n]);
                        guess[i] -= 2.0 * h;
                        fbmh = object_func(n_parameter, guess, xvalue, v_y[n]);
                        guess[i] += h;

                        df1 = (faph - famh) / (2.0 * h);
                        df2 = (fbph - fbmh) / (2.0 * h);
                        smatrix[i][j] += df1 * df2;
                    }
                }
            }
        }

        for (i = 1; i < n_parameter; ++i) {
            for (j = 0; j < i; ++j) { smatrix[i][j] = smatrix[j][i]; }
        }

        for (i = 0; i < n_parameter; ++i) {
            smatrix[i][i] = smatrix[i][i] * (1.0 + lamda);
            for (m = 0; m < n_x; ++m) {
                for (n = 0; n < n_y; ++n) {
                    xvalue = std::log(S / v_x[m]) + (cc * v_y[n]);
                    guess[i] += h;
                    faph = object_func(n_parameter, guess, xvalue, v_y[n]);
                    guess[i] -= 2.0 * h;
                    famh = object_func(n_parameter, guess, xvalue, v_y[n]);
                    guess[i] += h;

                    df1 = (faph - famh) / (2.0 * h);
                    known[i] -= df1 * perr[m][n];
                }
            }
        }

        gaussian_elimination(smatrix, known, unknown, n_parameter);
        for (i = 0; i < n_parameter; ++i) { tguess[i] = guess[i] + unknown[i]; }

        error = 0.0;
        for (i = 0; i < n_x; ++i) {
            for (j = 0; j < n_y; ++j) {
                xvalue = std::log(S / v_x[i]) + (cc * v_y[j]);
                fx_value = object_func(n_parameter, tguess, xvalue, v_y[j]);
                tperr[i][j] = (fx_value - mat[i][j]) / mat[i][j];
                error += tperr[i][j] * tperr[i][j];
            }
        }
        if (error > olderror) { lamda *= 10.0; } // iter != 0
        else {
            lamda /= 10.0;
            for (i = 0; i < n_parameter; ++i) { guess[i] = tguess[i]; }

            for (i = 0; i < n_x; ++i) {
                for (j = 0; j < n_y; ++j) { perr[i][j] = tperr[i][j]; }
            }

            tmperror = olderror;
            olderror = error;
        }

        if (error < minerror) {
            for (i = 0; i < n_parameter; ++i) { minguess[i] = guess[i]; }
            minerror = error;
        }

        std::cout << iter << " 번째" << std::endl;
        for (i = 0; i < n_parameter; ++i) { std::cout << "p[" << i << "]: " << guess[i] << std::endl; }
        std::cout << "error: " << minerror << std::endl;
        iter++;

    } while (error > tolerance && std::fabs(tmperror - error) > tolerance && iter < 10000);

    for (i = 0; i < n_parameter; ++i) { parameter[i] = minguess[i]; }
    for (i = 0; i < n_x; ++i) {
        delete perr[i];
        delete tperr[i];
    }

    delete[] perr;
    delete[] tperr;

    for (i = 0; i < n_parameter; ++i) { delete smatrix[i]; }

    delete[] unknown;
    delete[] known;
    delete[] guess;
    delete[] tguess;
    delete[] minguess;
    delete[] smatrix;

}

double local_volatility(double S, double cc, int n_x, double *v_x, int n_y, double *v_y, double **mat, int n_parameter,
                        double *parameter) {
    int i, j;
    double d1, iv, k, x, t, sigmat, sigmak, sigmak2;
    double *p;
    p = new double[n_parameter];
    for (i = 0; i < n_parameter; ++i) { p[i] = parameter[i]; }

    for (i = 0; i < n_x; ++i) {
        k = v_x[i];
        for (j = 0; j < n_y; ++j) {
            t = v_y[j];
            x = std::log(S / k) + (cc * t); // log(s*exp(cc*t)/k)
            iv = implied_volatility_function(n_parameter, parameter, x, t);
            d1 = (std::log(S / k) + (cc + 0.5 * iv * iv) * t) / (iv * std::sqrt(t));

            sigmat = p[2] * p[1] * std::exp(p[1] * t) +
                     cc * (p[3] + 2.0 * p[4] * x + 3.0 * p[5] * x * x + 4.0 * p[6] * x * x * x);

            sigmak = -(p[3] + 2.0 * p[4] * x + 3.0 * p[5] * x * x +
                       4.0 * p[6] * x * x * x) / k;
            sigmak2 = (p[3] + 2.0 * p[4] * (x + 1.0) + 3.0 * p[5] * x * (x + 2.0)
                       + 4.0 * p[6] * x * x * (x + 3.0)) / (k * k);

            mat[i][j] = (iv * iv + 2.0 * iv * t * sigmat + 2.0 * cc * iv * k * t * sigmak);
            mat[i][j] = mat[i][j] / (std::pow((1.0 + k * d1 * std::sqrt(t) * sigmak), 2.0)
                                     + iv * k * k * t * (sigmak2 - d1 * sigmak * sigmak * std::sqrt(t)));
            mat[i][j] = std::sqrt(std::fabs(mat[i][j]));
        }
    }

    delete[] p;
}

void
implied_volatility_2(double S, double cc, int n_x, double *v_x, int n_y, double *v_y, double **mat, int n_parameter,
                     double *parameter) {
    int i, j;
    double x_value;

    for (i = 0; i < n_x; ++i) {
        for (j = 0; j < n_y; ++j) {
            x_value = std::log(S / v_x[i]) + (cc * v_y[j]);
            mat[i][j] = implied_volatility_function(n_parameter,
                                                    parameter,
                                                    x_value, v_y[j]);
        }
    }
}
