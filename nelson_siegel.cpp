//
// Created by jewoo on 2022-10-29.
//

#include <iostream>
#include <cmath>

#include "nelson_siegel.h"
#include "mathlib.h"

void nelson_siegel_zerorate(int cal_type, int nztime, double *ztime, double *zrate, double *nsrate) {
    int i, nparameter;
    double *parameter;

    nparameter = 4;
    parameter = new double[nparameter];

    // caltype 1: bond price, 0:zerorate; // 비선형 파라미터 추정
    if (cal_type == 1) {
        get_ns_parameter_with_bondprice(nztime, ztime, zrate, nparameter, parameter);
    } else {
        get_ns_parameter_with_zerorate(nztime, ztime, zrate, nparameter, parameter);
    }

    for (i = 0; i < nztime; ++i) { nsrate[i] = ns_zerorate_function(parameter, ztime[i]); }

    delete[] parameter;
}

void get_ns_parameter_with_zerorate(int nztime, double *ztime, double *zrate, int nguess, double *guess) {
    std::cout << "parameter calibration with zero rate" << std::endl;

    int i, k, j, iter;
    double *unknown, *known;
    double **hessian;
    double *tguess, *minguess;
    double error, tolerance, h;
    double *rfrate, *err, *terr, *perr, *tperr;

    rfrate = new double[nztime];
    err = new double[nztime];
    terr = new double[nztime];
    perr = new double[nztime];
    tperr = new double[nztime];

    unknown = new double[nguess];
    known = new double[nguess];
    tguess = new double[nguess];
    minguess = new double[nguess];

    guess[0] = zrate[nztime - 1] * 1.1;
    guess[1] = zrate[0] * 0.9 - zrate[nztime - 1];
    guess[2] = 0.001;
    guess[3] = 0.5;

    // 파라미터 계산 시작
    for (i = 0; i < nztime; ++i) { rfrate[i] = zrate[i]; } // 절대값 계산

    hessian = new double *[nguess];
    for (i = 0; i < nguess; ++i) { hessian[i] = new double[nguess]; }
    double lambda = 0.001;
    double df1, df2, faph, famh, fbph, fbmh, olderror, tmperror, minerror, nsrate;

    olderror = 0.0;
    for (i = 0; i < nztime; ++i) {
        nsrate = ns_zerorate_function(guess, ztime[i]);
        perr[i] = rfrate[i] - nsrate; // 절대오차 계산, Y - Fx
        err[i] = perr[i] / rfrate[i]; // 상대오차 계산 (Y- Fx)/Y
        olderror += err[i] * err[i]; // 최초 누적 에러 계산
    }

    tmperror = 0.0;
    minerror = 1.0e10;
    tolerance = 1.0e-15;
    h = 0.0001;
    lambda = 0.001;
    iter = 0;

    do {

        for (i = 0; i < nguess; ++i) {
            for (j = 0; j < nguess; ++j) { hessian[i][j] = 0.0; }
            known[i] = 0.0;
        }

        for (i = 0; i < nguess; ++i) {
            for (j = i; j < nguess; ++j) {
                for (k = 0; k < nztime; ++k) {
                    faph = ns_zerorate_function(guess, ztime[k], i, h);
                    famh = ns_zerorate_function(guess, ztime[k], i, -h);
                    fbph = ns_zerorate_function(guess, ztime[k], j, h);
                    fbmh = ns_zerorate_function(guess, ztime[k], j, -h);
                    df1 = (faph - famh) / (2.0 * h); // round Fi / round xi
                    df2 = (fbph - fbmh) / (2.0 * h); // round Fj / round xj
                    hessian[i][j] += df1 * df2;
                }
            }
        }

        for (i = 1; i < nguess; ++i) {
            for (j = 0; j < i; ++j) { hessian[i][j] = hessian[j][i]; }
        }

        for (i = 0; i < nguess; ++i) {
            hessian[i][i] = hessian[i][i] * (1.0 + lambda);
            for (k = 0; k < nztime - 1; ++k) {
                faph = ns_zerorate_function(guess, ztime[k], i, h);
                famh = ns_zerorate_function(guess, ztime[k], i, -h);
                df1 = (faph - famh) / (2.0 * h);
                known[i] += df1 * perr[k]; // (round Fi/ round xi) * 절대오차
            }
        }
        // 가우스 소거법에 의해 unknown 값 계산
        gaussian_elimination(hessian, known, unknown, nguess);

        // Target Guess 는 guess + unknwon 값
        for (i = 0; i < nguess; ++i) { tguess[i] = guess[i] + unknown[i]; }

        for (i = 0; i < nguess; ++i) {
            // Target Guess값이 모수의 법위를 초과할 땐 제약 조건 부여
            if (i == 0 || i == 3) {
                if (tguess[i] >= 1.0) {
                    tguess[i] = 0.999;
                    unknown[i] = 0.999 - guess[i];
                } else if (tguess[i] <= 0.0) {
                    tguess[i] = 0.001;
                    unknown[i] = 0.001 - guess[i];
                }
            } else if (i == 1) {
                if (tguess[i] >= 1.0) {
                    tguess[i] = 0.999;
                    unknown[i] = 0.999 - guess[i];
                } else if (tguess[i - 1] + tguess[i] <= 0.0) {
                    tguess[i] = 0.001 - tguess[i - 1];
                    unknown[i] = 0.001 - tguess[i - 1] - guess[i];
                }
            } else if (i == 2) {
                if (tguess[i] >= 1.0) {
                    tguess[i] = 0.999;
                    unknown[i] = 0.999 - guess[i];
                } else if (tguess[i] <= -1.0) {
                    tguess[i] = -0.999;
                    unknown[i] = -0.999 - guess[i];
                }
            }
        }

        error = 0.0;
        for (i = 0; i < nztime; ++i) {
            nsrate = ns_zerorate_function(tguess, ztime[i]);
            // Target Guess에 의한 가격
            tperr[i] = rfrate[i] - nsrate; // target Geuss에 의한 절대오차 계산
            terr[i] = tperr[i] / rfrate[i]; // target Geuss에 의한 상대오차 계산
            error += terr[i] * terr[i]; // target Geuss에 의한 누적에러 계산
        }

        if (error >= olderror) { lambda *= 10.0; } // Target Guess에 의한 누적 에러가 이전에 계산된 에러보다 크다면 lambda 값 10배 확대 조정
        else {
            lambda /= 10.0; // Target Guess에 의한 누적 에러가 이전에 계산된 에러보다 작다면 lambda 값 10배 축소 조정
            for (i = 0; i < nguess; ++i) { guess[i] = tguess[i]; } // Guess는 Target Guess로 수정
            for (i = 0; i < nztime - 1; ++i) { perr[i] = tperr[i]; }

            tmperror = olderror;
            olderror = error; // Guess에 의한 에러를 Target Guess에 의한 에러로 수정
        }

        if (error < minerror) {
            // 에러값이 가장 작다면 최소에러의 Guess값을  Target Guess값으로 수정
            for (i = 0; i < nguess; ++i) { minguess[i] = guess[i]; }
            minerror = error;
        }

        iter++;
    } while (error > tolerance && std::fabs(error - tmperror) > 1.0e-15 && iter < 1000);

    for (i = 0; i < nguess; ++i) { guess[i] = minguess[i]; }
    std::cout << "Nelson-sigel 파라미터 계산" << std::endl;
    for (i = 0; i < nguess; ++i) { std::cout << guess[i] << std::endl; }

    for (i = 0; i < nguess; ++i) { delete hessian[i]; }

    delete[] hessian;
    delete[] unknown;
    delete[] known;
    delete[] tguess;
    delete[] minguess;
    delete[] rfrate;
    delete[] err;
    delete[] terr;
    delete[] tperr;
    delete[] perr;

}

void get_ns_parameter_with_bondprice(int nztime, double *ztime, double *zrate, int nguess, double *guess) {
    std::cout << "parameter calibration with bond price" << std::endl;

    int i, k, j, iter;
    double *unknown, *known;
    double **hessian;
    double *tguess, *minguess;
    double error, tolerance, h;
    double *rfprice, *err, *terr, *perr, *tperr;

    rfprice = new double[nztime];
    err = new double[nztime];
    terr = new double[nztime];
    perr = new double[nztime];
    tperr = new double[nztime];

    unknown = new double[nguess];
    known = new double[nguess];
    tguess = new double[nguess];
    minguess = new double[nguess];

    guess[0] = zrate[nztime - 1] * 1.1;
    guess[1] = zrate[0] * 0.9 - zrate[nztime - 1];
    guess[2] = 0.001;
    guess[3] = 0.5;

    // 파라미터 계산 시작
    for (i = 0; i < nztime; ++i) { rfprice[i] = std::exp(-ztime[i] * zrate[i]); } // 절대값 계산

    hessian = new double *[nguess];
    for (i = 0; i < nguess; ++i) { hessian[i] = new double[nguess]; }

    double lambda;
    double df1, df2, faph, famh, fbph, fbmh;
    double olderror, tmperror, minerror, nsprice;

    olderror = 0.0;

    for (i = 0; i < nztime; ++i) {
        nsprice = price_ns_discountfactor(guess, ztime[i]);
        perr[i] = rfprice[i] - nsprice; // 절대오차 계산, Y - Fx
        err[i] = perr[i] / rfprice[i]; // 상대오차 계산 (Y- Fx)/Y
        olderror += err[i] * err[i]; // 최초 누적 에러 계산
    }

    tmperror = 0.0;
    minerror = 1.0e10;
    tolerance = 1.0e-15;
    h = 0.0001;
    lambda = 0.001;
    iter = 0;

    do {

        for (i = 0; i < nguess; ++i) {
            for (j = 0; j < nguess; ++j) { hessian[i][j] = 0.0; }
            known[i] = 0.0;
        }

        for (i = 0; i < nguess; ++i) {
            for (j = i; j < nguess; ++j) {
                for (k = 0; k < nztime; ++k) {
                    faph = price_ns_discountfactor(guess, ztime[k], i, h);
                    famh = price_ns_discountfactor(guess, ztime[k], i, -h);
                    fbph = price_ns_discountfactor(guess, ztime[k], j, h);
                    fbmh = price_ns_discountfactor(guess, ztime[k], j, -h);
                    df1 = (faph - famh) / (2.0 * h); // round Fi / round xi
                    df2 = (fbph - fbmh) / (2.0 * h); // round Fj / round xj
                    hessian[i][j] += df1 * df2;
                }
            }
        }

        for (i = 1; i < nguess; ++i) {
            for (j = 0; j < i; ++j) { hessian[i][j] = hessian[j][i]; }
        }

        for (i = 0; i < nguess; ++i) {
            hessian[i][i] = hessian[i][i] * (1.0 + lambda);
            for (k = 0; k < nztime - 1; ++k) {
                faph = price_ns_discountfactor(guess, ztime[k], i, h);
                famh = price_ns_discountfactor(guess, ztime[k], i, -h);
                df1 = (faph - famh) / (2.0 * h);
                known[i] += df1 * perr[k]; // (round Fi/ round xi) * 절대오차
            }
        }

        // 가우스 소거법에 의해 unknown 값 계산
        gaussian_elimination(hessian, known, unknown, nguess);

        // Target Guess 는 guess + unknwon 값
        for (i = 0; i < nguess; ++i) { tguess[i] = guess[i] + unknown[i]; }

        for (i = 0; i < nguess; ++i) {
            // Target Guess값이 모수의 법위를 초과할 땐 제약 조건 부여
            if (i == 0 || i == 3) {
                if (tguess[i] >= 1.0) {
                    tguess[i] = 0.999;
                    unknown[i] = 0.999 - guess[i];
                } else if (tguess[i] <= 0.0) {
                    tguess[i] = 0.001;
                    unknown[i] = 0.001 - guess[i];
                }
            } else if (i == 1) {
                if (tguess[i] >= 1.0) {
                    tguess[i] = 0.999;
                    unknown[i] = 0.999 - guess[i];
                } else if (tguess[i - 1] + tguess[i] <= 0.0) {
                    tguess[i] = 0.001 - tguess[i - 1];
                    unknown[i] = 0.001 - tguess[i - 1] - guess[i];
                }
            } else if (i == 2) {
                if (tguess[i] >= 1.0) {
                    tguess[i] = 0.999;
                    unknown[i] = 0.999 - guess[i];
                } else if (tguess[i] <= -1.0) {
                    tguess[i] = -0.999;
                    unknown[i] = -0.999 - guess[i];
                }
            }
        }

        error = 0.0;
        for (i = 0; i < nztime; ++i) {
            nsprice = price_ns_discountfactor(tguess, ztime[i]);
            // Target Guess에 의한 가격
            tperr[i] = rfprice[i] - nsprice; // target Geuss에 의한 절대오차 계산
            terr[i] = tperr[i] / rfprice[i]; // target Geuss에 의한 상대오차 계산
            error += terr[i] * terr[i]; // target Geuss에 의한 누적에러 계산
        }

        if (error >= olderror) { lambda *= 10.0; } // Target Guess에 의한 누적 에러가 이전에 계산된 에러보다 크다면 lambda 값 10배 확대 조정
        else {
            lambda /= 10.0; // Target Guess에 의한 누적 에러가 이전에 계산된 에러보다 작다면 lambda 값 10배 축소 조정
            for (i = 0; i < nguess; ++i) { guess[i] = tguess[i]; } // Guess는 Target Guess로 수정
            for (i = 0; i < nztime - 1; ++i) { perr[i] = tperr[i]; }

            tmperror = olderror;
            olderror = error; // Guess에 의한 에러를 Target Guess에 의한 에러로 수정
        }

        if (error < minerror) {
            // 에러값이 가장 작다면 최소에러의 Guess값을  Target Guess값으로 수정
            for (i = 0; i < nguess; ++i) { minguess[i] = guess[i]; }
            minerror = error;
        }

        iter++;
    } while (error > tolerance && std::fabs(error - tmperror) > 1.0e-15 && iter < 1000);

    for (i = 0; i < nguess; ++i) { guess[i] = minguess[i]; }
    std::cout << "Nelson-sigel 파라미터 계산" << std::endl;
    for (i = 0; i < nguess; ++i) { std::cout << guess[i] << std::endl; }

    for (i = 0; i < nguess; ++i) { delete hessian[i]; }

    delete[] hessian;
    delete[] unknown;
    delete[] known;
    delete[] tguess;
    delete[] minguess;
    delete[] rfprice;
    delete[] err;
    delete[] terr;
    delete[] tperr;
    delete[] perr;
}

double price_ns_discountfactor(double *guess, double time) {
    double nsspot = ns_zerorate_function(guess, time);
    return std::exp(-time * nsspot);
}

double price_ns_discountfactor(double *guess, double time, int j, double h) {
    int i;
    double nsspot, p[4];
    for (i = 0; i < 4; ++i) { p[i] = guess[i]; }

    p[j] = p[j] + h;
    nsspot = ns_zerorate_function(p, time);
    return std::exp(-time * nsspot);
}

double ns_zerorate_function(double *parameter, double time) {
    double nsspot = parameter[0] + parameter[1] * (1 - std::exp(-parameter[3] * time)) / (parameter[3] * time)
                    + parameter[2] *
                      ((1 - std::exp(-parameter[3] * time)) / (parameter[3] * time) -
                       std::exp(-parameter[3] * time));
    return nsspot;
}

double ns_zerorate_function(double *parameter, double time, int j, double h) {
    int i;
    double nsspot, p[4];
    for (i = 0; i < 4; ++i) { p[i] = parameter[i]; }
    p[j] = p[j] + h;

    nsspot = p[0] + p[1] * (1 - std::exp(-p[3] * time)) / (p[3] * time)
             + p[2] * ((1 - std::exp(-p[3] * time)) / (p[3] * time) - std::exp(-p[3] * time));

    return nsspot;
}
