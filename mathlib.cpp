//
// Created by jewoo on 2022-10-13.
//

#include "mathlib.h"
#include <iostream>
#include <ctgmath>
#include <random>
#include <chrono>

void gaussian_elimination(double **Smatrix, double *Known, double *Unknown, int n_eqns) {

    double ratio, sum, maxs, tmp;
    int i, k, j, maxi;

    for (k = 0; k < n_eqns - 1; ++k) {
        maxs = std::fabs(Smatrix[k][k]);
        maxi = k;
        for (i = k + 1; i < n_eqns; ++i) {
            if (maxs < std::fabs(Smatrix[i][k]))
                maxi = i;
        }
        if (maxi != k) {
            for (i = 0; i < n_eqns; ++i) {
                tmp = Smatrix[k][i];
                Smatrix[k][i] = Smatrix[maxi][i];
                Smatrix[maxi][i] = tmp;
            }
            tmp = Known[k];
            Known[k] = Known[maxi];
            Known[maxi] = tmp;
        }
        for (i = k + 1; i < n_eqns; ++i) {
            if (Smatrix[i][k] == 0.) continue;
            ratio = -Smatrix[i][k] / Smatrix[k][k];
            for (j = k + 1; j < n_eqns; ++j) {
                if (Smatrix[k][j] == 0.) continue;
                Smatrix[i][j] += ratio * Smatrix[k][j];
            }
            Known[i] += ratio * Known[k];
        }
    }

    Unknown[n_eqns - 1] = Known[n_eqns - 1] / Smatrix[n_eqns - 1][n_eqns - 1];
    for (i = n_eqns - 2; i >= 0; --i) {
        sum = 0.;
        for (j = i + 1; j < n_eqns; ++j) sum += Smatrix[i][j] * Unknown[j];
        Unknown[i] = (Known[i] - sum) / Smatrix[i][i];
    }
}

void tridiagonal_elimination(double **Smatrix, double *Known, double *Unknown, int n_eqns) {
    double ratio, sum;
    int i, k, j;

    // 첫번째 행부터 마지막 직전행까지는 동일한 로직
    for (k = 0; k < n_eqns - 1; ++k) {
        if (Smatrix[k][k] == 0.) {
            std::cout << "Matrix is singular in tridiagonal_elimination()" << std::endl;
            exit(0);
        }
        ratio = -Smatrix[k + 1][k] / Smatrix[k][k];
        if (k == n_eqns - 2) {
            for (j = k; j <= k + 1; j++) { // 2 열에 적용; 마지막 행에서는 2열만 적용하므로
                if (Smatrix[k][j] == 0.) continue;
                Smatrix[k + 1][j] += ratio * Smatrix[k][j];
            }
        } else { // k < n_eqns -2
            for (j = k; j <= k + 2; ++j) { // 3 열에 적용;
                if (Smatrix[k][j] == 0.) continue;
                Smatrix[k + 1][j] += ratio * Smatrix[k][j];
            }
        }
        Known[k + 1] += ratio * Known[k];
    }
    Unknown[n_eqns - 1] = Known[n_eqns - 1] / Smatrix[n_eqns - 1][n_eqns - 1];
    // 마지막행에서 바로 미지수 해 찾음
    for (i = n_eqns - 2; i >= 0; --i) {
        // 1열에 적용, 모든행에서 미지수 1개, 기지수 1개이므로
        sum = Smatrix[i][i + 1] * Unknown[i + 1];
        Unknown[i] = (Known[i] - sum) / Smatrix[i][i];
    }
}

void cholesky_decomposition(double **L, double **A, int n_col) {
    // A : in, L : out, n_col: 행렬 크기
    int i, j, k;
    double sum;
    for (i = 0; i < n_col; ++i) {
        for (j = 0; j < n_col; ++j) {
            L[i][j] = 0.0;
        }
    }
    L[0][0] = std::sqrt(A[0][0]);
    for (i = 1; i < n_col; ++i) {
        L[i][0] = A[i][0] / L[0][0];
        for (j = 0; j <= i; ++j) {
            sum = 0.0;
            if (j != i) {
                for (k = 0; k < j; ++k) {
                    sum += L[i][k] * L[j][k];
                }
                L[i][j] = (A[i][j] - sum) / L[j][j];
            } else {
                for (k = 0; k < i; ++k) {
                    sum += std::pow(L[i][k], 2);
                }
                L[i][j] = std::sqrt(A[i][j] - sum);
            }
        }
    }
}

void cholesky_elimination(double **L, double *UV, double *KV, int n_col) {
    double sum;
    for (int i = 0; i < n_col; ++i) {
        UV[i] = 0.0;
    }
    for (int i = 0; i < n_col; ++i) {
        sum = KV[i];
        for (int j = 0; j < i; ++j) {
            sum -= L[i][j] * UV[j];
        }
        UV[i] = sum / L[i][i];
    }
    for (int i = n_col - 1; i >= 0; --i) {
        sum = UV[i];
        for (int j = i + 1; j < n_col; ++j) {
            sum -= L[j][i] * UV[j]; // L[j][i] = LT[i][j]
        }
        UV[i] = sum / L[i][i];
    }
}

void cholesky_solver(double **A, double *UV, double *KV, int n_col) {
    double **L; // 하삼각 행렬

    L = new double *[n_col];
    for (int i = 0; i < n_col; ++i) {
        L[i] = new double[n_col];
    }
    cholesky_decomposition(L, A, n_col);

    std::cout << "하삼각 행렬" << std::endl;
    for (int i = 0; i < n_col; ++i) {
        for (int j = 0; j < n_col; ++j) {
            std::cout << L[i][j] << "  ";
        }
        std::cout << std::endl;
    }

    std::cout << "상삼각 행렬" << std::endl;
    for (int i = 0; i < n_col; ++i) {
        for (int j = 0; j < n_col; ++j) {
            std::cout << L[j][i] << "  ";
        }
        std::cout << std::endl;
    }

    cholesky_elimination(L, UV, KV, n_col);
    std::cout << "숄레스키분해를 이용한 연립방정식의 해는" << std::endl;
    for (int i = 0; i < n_col; ++i) {
        std::cout << "x_" << i + 1 << " : " << UV[i] << std::endl;
    }
    for (int i = 0; i < n_col; ++i) {
        delete[]L[i];
    }
    delete[] L;
}

double nonlinear_function(double X) {
    return std::exp(0.5 * X) - 5.0;
}

double bisection_method() {
    double xa, xb, xc, fxa, fxc;
    double error_tolerance, error;
    int iter;

    iter = 0;
    error_tolerance = 1.0e-15;

    xa = -10.0; // 해외 구간 최소값 -10
    xb = 10.0; // 해외 구간 최대값 10
    do {

        xc = (xa + xb) / 2.0;
        fxa = nonlinear_function(xa);
        fxc = nonlinear_function(xc);
        if (fxa * fxc < 0) {
            xb = xc;
        } else if (fxa * fxc > 0) {
            xa = xc;
        } else {
            error = 0;
            break;
        }
        error = std::fabs(xb - xa);
        iter++;
    } while (error > error_tolerance && iter < 100);
    if (iter == 100) {
        std::cout << " 허용오차범위는 불만족하나, 가장 근사한 값은 " << xc << " 입니다." << std::endl;
    } else {
        std::cout << "bisection method 계산횟수 : " << iter << std::endl;
    }
    return xc;
}

double newtonraphson_method() {
    double fx, fxprime;
    double x, dx, error, error_tolerance;
    int iter;
    iter = 0;
    error_tolerance = 1.0e-15;
    x = 5.0;
    dx = 0.001;
    do {

        fx = nonlinear_function(x);
        fxprime = (nonlinear_function(x + dx) - nonlinear_function(x - dx)) / (2.0 * dx);
        x = x - fx / fxprime;
        error = std::fabs(fx);
        iter++;
    } while (error > error_tolerance && iter < 100);
    if (iter == 100) {
        std::cout << " 허용오차범위는 불만족하나, 가장 근사한 값은 " << x << " 입니다." << std::endl;
    } else {
        std::cout << "newtonraphson method 계산횟수 : " << iter << std::endl;
    }
    return x;

}

double linear_interpolation(int n_data, int *x_data, double *y_data, int x) {
    int i;
    double y = 0;
    for (i = 0; i < n_data; ++i) {
        if (x_data[i] > x) {
            y = y_data[i - 1] + (y_data[i] - y_data[i - 1]) / (x_data[i] - x_data[i - 1]) * (x - x_data[i - 1]);
            break;
        } else
            continue;
    }
    return y;
}

void set_system_matrix(int n_data, double **Smatrix, double *known, int *x_data, double *y_data) {
    int i, j;
    double *hi;
    hi = new double[n_data - 1];

    for (i = 0; i < n_data - 2; ++i) {
        for (j = 0; j < n_data - 2; ++j) {
            Smatrix[i][j] = 0.0;
        }
        hi[i] = x_data[i + 1] - x_data[i];
    }

    hi[n_data - 2] = x_data[n_data - 1] - x_data[n_data - 2];

    // 0~n_data-3는 n_data-2개
    Smatrix[0][0] = 2 * (hi[0] + hi[1]);
    Smatrix[0][1] = hi[1];
    known[0] = 3 * ((y_data[2] - y_data[1]) / hi[1] - (y_data[1] - y_data[0]) / hi[0]);

    for (i = 1; i < n_data - 3; ++i) {
        Smatrix[i][i - 1] = hi[i];
        Smatrix[i][i] = 2 * (hi[i] + hi[i + 1]);
        Smatrix[i][i + 1] = hi[i + 1];
        known[i] = 3 * ((y_data[i + 2] - y_data[i + 1]) / hi[i + 1] - (y_data[i + 1] - y_data[i]) / hi[i]);
    }

    Smatrix[n_data - 3][n_data - 4] = hi[n_data - 3];
    Smatrix[n_data - 3][n_data - 3] = 2 * (hi[n_data - 3] + hi[n_data - 2]);

    known[n_data - 3] = 3 * ((y_data[n_data - 1] - y_data[n_data - 2]) / hi[n_data - 2] -
                             (y_data[n_data - 2] - y_data[n_data - 3]) / hi[n_data - 3]);

    delete[]hi;
}

double cubicspline_interpolation(int n_data, int *x_data, double *y_data, int x) {
    int i, n_eq;
    double **s_matrix; // 행렬선언
    double *known_value; // 기지값
    double *unknown_value; // 미지값

    n_eq = n_data - 2;
    s_matrix = new double *[n_eq];

    for (i = 0; i < n_eq; ++i) {
        s_matrix[i] = new double[n_eq];
    }

    known_value = new double[n_eq];
    unknown_value = new double[n_eq];

    set_system_matrix(n_data, s_matrix, known_value, x_data, y_data);
    tridiagonal_elimination(s_matrix, known_value, unknown_value, n_eq);

    double *ci;
    ci = new double[n_data];
    for (i = 0; i < n_eq; ++i) {
        ci[i + 1] = unknown_value[i];
    }
    ci[0] = 0.0; // 자연 경계 조건
    ci[n_data - 1] = 0.0;// 자연 경계 조건
    double a, b, d, y, hi, xp;

    if (x_data[0] > x) {
        hi = x_data[1] - x_data[0];
        xp = x - x_data[0];
        a = y_data[0];
        b = (y_data[1] - y_data[0]) / hi - hi * (2 * ci[0] + ci[1]) / 3;
        d = (ci[1] - ci[0]) / (3 * hi);
        y = a + b * xp + ci[0] * xp * xp + d * xp * xp * xp;
    } else if (x_data[n_data - 1] < x) {
        hi = x_data[n_data - 1] - x_data[n_data - 2];
        xp = x - x_data[n_data - 2];
        a = y_data[n_data - 2];
        b = (y_data[n_data - 1] - y_data[n_data - 2]) / hi - hi * (2 * ci[n_data - 2] + ci[n_data - 1]) / 3;
        d = (ci[n_data - 1] - ci[n_data - 2]) / (3 * hi);
        y = a + b * xp + ci[n_data - 2] * xp * xp + d * xp * xp * xp;
    } else {
        for (i = 1; i < n_data; ++i) {
            if (x_data[i] > x) {
                hi = x_data[i] - x_data[i - 1];
                xp = x - x_data[i - 1];
                a = y_data[i - 1];
                b = (y_data[i] - y_data[i - 1]) / hi - hi * (2 * ci[i - 1] + ci[i]) / 3;
                d = (ci[i] - ci[i - 1]) / (3 * hi);
                y = a + b * xp + ci[i - 1] * xp * xp + d * xp * xp * xp;
                break;
            } else continue;
        }
    }


    for (i = 0; i < n_eq; ++i) {
        delete s_matrix[i];
    }

    delete[] s_matrix, known_value, unknown_value;

    return y;

}

double normdistrand() {
    double rnd = 0;
    /* Initialise. Do this once (not for every
     random number). */
    std::random_device rd;
    std::mt19937_64 gen(rd());

    /* This is where you define the number generator for double: */
    std::uniform_real_distribution<double> dis(0, RAND_MAX);

    for (int i = 0; i < 12; ++i) {
        rnd += dis(gen) / static_cast<double>(RAND_MAX);
    }
    rnd -= 6.0;
    return rnd;
}

void normal_distribution_goodness_fit_test(unsigned long nrand, double *randnum) {
    unsigned long i;
    double mean, stddev;

    mean = 0.0;
    for (i = 0; i < nrand; ++i) {
        mean += randnum[i];
    }

    mean /= nrand;

    stddev = 0.0;
    for (i = 0; i < nrand; ++i) {
        stddev += std::pow((randnum[i] - mean), 2.0);
    }

    stddev /= (nrand - 1);
    stddev = std::sqrt(stddev);

    std::cout << "평균 : " << mean << std::endl;
    std::cout << "표준편차 : " << stddev << std::endl;

    double mu3, mu4, skewness, kurtosis;
    mu3 = mu4 = 0.0;

    for (i = 0; i < nrand; ++i) {
        mu3 += std::pow((randnum[i] - mean), 3.0);
        mu4 += std::pow((randnum[i] - mean), 4.0);
    }

    skewness = (mu3 / nrand) / std::pow(stddev, 3.0);
    kurtosis = (mu4 / nrand) / std::pow(stddev, 4.0);

    std::cout << "왜도 : " << skewness << std::endl;
    std::cout << "첨도 : " << kurtosis << std::endl;

    double jarque_bera;
    jarque_bera = nrand * (skewness * skewness + std::pow((kurtosis - 3.0), 2.0) / 4.0) / 6.0;

    std::cout << "Jarque-Bera Test: " << jarque_bera << std::endl;

}

double inverse_normal_cumulative_distribution_function(double p) {
    static const double a[] = {
            -3.969683028665376e+01, 2.209460984245205e+02,
            -2.759285104469687e+02, 1.383577518672690e+02,
            -3.066479806614716e+01, 2.506628277459239e+00
    };

    static const double b[] = {
            -5.447609879822406e+01, 1.615858368580409e+02,
            -1.556989798598866e+02, 6.680131188771872e+01,
            -1.328068155288572e+01
    };

    static const double c[] = {
            -7.784894002430293e-03, -3.223964580411365e-01,
            -2.400758277161838e+00, -2.549732539343734e+00,
            4.374664141464968e+00, 2.938163982698783e+00

    };

    static const double d[] = {
            7.784695709041462e+03, 3.224671290700398e+01,
            2.445134137142996e+00, 3.754408661907416e+00,
    };


    static double LOW{0.02425};
    static double HIGH{0.97575};

    double q, r;
    errno = 0;

    if (p < LOW) {
        if (p < EPS) {
            p = EPS;
        }
        // Rational approximation for lower region
        q = std::sqrt(-2 * std::log(p));
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
    } else if (p > HIGH) {
        // Rational approximation for upper region
        if (p > 1.0 - EPS) {
            p = 1.0 - EPS;
        }
        q = std::sqrt(-2 * std::log(1 - p));
        return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);

    } else {
        // Rational approximation for central region
        q = p - 0.5;
        r = q * q;
        return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
               (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + +b[4]) * r + 1);
    }
}

//double normdistrand_BoxMuller() {
//    double u1, u2, z1, z2;
//    static int iset{0};
//    static double gset;
//
//    std::random_device rd;
//    std::mt19937_64 gen(rd());
//    std::uniform_real_distribution<double> dis;
//
//    if (iset == 0) {
//        do {
//            u1 = dis(gen) / static_cast<double>(RAND_MAX);
//            u2 = dis(gen) / static_cast<double>(RAND_MAX);
//            z1 = std::sqrt(-2 * std::log(u1)) * std::cos(2 * PI * u2);
//            z2 = std::sqrt(-2 * std::log(u1)) * std::sin(2 * PI * u2);
//        } while (u1 == 0.0);
//        gset = z2; // 다음에 사용하기 위해서
//        iset = 1;
//        return z1;
//    } else {
//        iset = 0;
//        return gset; // 저장된 값
//    }
//}

double normdistrand_BoxMuller() {
    double u1, u2, z1, z2;
    static int iset{0};
    static double gset;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
//    std::random_device rd;
//    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0, RAND_MAX);

    if (iset == 0) {
        do {
            u1 = dis(generator) / static_cast<double>(RAND_MAX);
            u2 = dis(generator) / static_cast<double>(RAND_MAX);
            z1 = std::sqrt(-2 * std::log(u1)) * std::cos(2 * PI * u2);
            z2 = std::sqrt(-2 * std::log(u1)) * std::sin(2 * PI * u2);
        } while (u1 == 0.0);
        gset = z2; // 다음에 사용하기 위해서
        iset = 1;
        return z1;
    } else {
        iset = 0;
        return gset; // 저장된 값
    }
}

double N(double z) {
    double a1, a2, a3, a4, a5, r, c;

    a1 = 0.31938153;
    a2 = -0.356563782;
    a3 = 1.781477937;
    a4 = -1.821255978;
    a5 = 1.330274429;
    r = 0.2316419;
    c = 1 / std::sqrt(2 * PI);

    if (z > 100.0) { return 1.0; }
    else if (z < -100.0) { return 0.0; }


    double x = std::fabs(z);
    double k = 1.0 / (1.0 + r * x);
    double b = c * std::exp((-z * z) / 2.0);
    double nv = ((((a5 * k + a4) * k + a3) * k + a2) * k + a1) * k;
    nv = 1.0 - b * nv;
    if (z < 0.0) { nv = 1.0 - nv; }
    return nv;
}

void mean_stddev_error(unsigned long nrand, double *value) {
    unsigned long i;
    double mean, stddev;
    mean = 0.0;
    for (i = 0; i < nrand; ++i) {
        mean += value[i];
    }
    mean /= nrand;

    stddev = 0.0;
    for (i = 0; i < nrand; ++i) {
        stddev += std::pow((value[i] - mean), 2.0);
    }
    stddev /= (nrand - 1);
    stddev = std::sqrt(stddev);

    std::cout << "평균: " << mean << std::endl;
    std::cout << "표준 편차: " << stddev << std::endl;
    std::cout << "표준 오차: " << stddev / std::sqrt(static_cast<double>(nrand)) << std::endl;
}

double *halton_sequence(unsigned int num, unsigned int prime_number) {
    // 홀튼 수열 생성을 위하여 Brandimarte(2002), P236을 C++ 코드로 변환
    int dnum, digitnumber;
    unsigned i;
    double invprime;
    double haltonnumber;
    double *POINTS;
    POINTS = new double[num];

    for (i = 0; i < num; ++i) {
        invprime = 1.0 / static_cast<double>(prime_number);
        dnum = i;
        haltonnumber = 0.0;
        do {
            digitnumber = dnum % prime_number;
            haltonnumber += digitnumber * invprime;
            dnum = (dnum - digitnumber) / prime_number;
            invprime /= static_cast<double>(prime_number);
        } while (dnum > 0);

        POINTS[i] = haltonnumber;
    }

    return POINTS;
}

double halton_sequence_(int num, int prime_number) {
    // 홀튼 수열 생성을 위하여 Brandimarte(2002), P236을 C++ 코드로 변환
    int dnum, digitnumber;
    double invprime, haltonnumber;

    invprime = 1.0 / static_cast<double>(prime_number);
    dnum = num;
    haltonnumber = 0.0;

    do {
        digitnumber = dnum % prime_number;
        haltonnumber += digitnumber * invprime;
        dnum = (dnum - digitnumber) / prime_number;
        invprime /= static_cast<double>(prime_number);
    } while (dnum > 0);

    return haltonnumber;
}

double n(double z) { // The normal distribution function
    return 1.0 / std::sqrt(2.0 * PI) * std::exp(-0.5 * z * z);
}

double implied_volatility_bisection(double spot, double strike, double riskfree, double dividend, double maturity,
                                    double option_price) {
    double xa, xb, xc, fxa, fxc;
    double vol, error_tolerance, error;
    int iter = 0;
    error_tolerance = 1.0e-15;

    xa = 0.0; // 변동성 최저값 0
    xb = 10.0; // 변동성 최대값 10

    do {
        xc = (xa + xb) / 2.0;
        fxa = european_calloption_price(spot, strike, riskfree, dividend, xa, maturity);
        fxc = european_calloption_price(spot, strike, riskfree, dividend, xc, maturity);
        if (fxa * fxc < 0) {
            xb = xc;
        } else if (fxa * fxc > 0) {
            xa = xc;
        } else {
            error = 0.0;
            break;
        }
        error = std::fabs(xb - xa);
        iter++;
    } while (error > error_tolerance && iter < 100);

    if (iter == 100) {
        std::cout << "허용오차범위는 불만족하나, 가장 근사한 값은 " << xc << "입니다." << std::endl;
    } else {
        std::cout << "bisection method 계산횟수: " << iter << "입니다." << std::endl;
    }

    vol = xc;
    return vol;
}

double
european_calloption_price(double spot, double strike, double riskfree, double dividend, double vol, double maturity) {
    double d1, d2;
    d1 = (std::log(spot / strike) + (riskfree - dividend + (vol * vol) / 2.0) * maturity);
    d2 = d1 - vol * std::sqrt(maturity);

    return spot * std::exp(-dividend * maturity) * N(d1) -
           strike * std::exp(-riskfree * maturity) * N(d2); // BS call option price
}

double european_calloption_vega(double spot, double strike, double riskfree, double dividend, double volatility,
                                double maturity) {
    double d1 = (std::log(spot / strike) +
                 (riskfree - dividend + (volatility * volatility) / 2.0) * maturity)
                / (volatility * std::sqrt(maturity));

    return spot * std::exp(-dividend * maturity) * n(d1) * std::sqrt(maturity); // BS call option vega
}

double implied_volatility_newtonraphson(double spot, double strike, double riskfree, double dividend, double maturity,
                                        double option_price) {
    double fx, fxprime;
    double vol, error_tolerance, error;
    int iter;

    iter = 0;

    error_tolerance = 1.0e-15;
    vol = 0.5;

    do {
        fx = european_calloption_price(spot, strike, riskfree, dividend, vol, maturity);
        fxprime = european_calloption_vega(spot, strike, riskfree, dividend, vol, maturity);

        vol = vol - (fx - option_price) / fxprime;
        error = std::fabs(fx - option_price);
        iter++;
    } while (error > error_tolerance && iter < 100);

    if (iter == 100) {
        std::cout << "허용오차범위는 불만족하나, 가장 근사한 값은 " << vol << "입니다." << std::endl;
    } else {
        std::cout << "newton raphson method 계산횟수: " << iter << "입니다." << std::endl;
    }

    return vol;
}

void
par_to_zero_bootstrapping(int nytime, double *ytime, double *yrate, int coupon_frequency, int nztime, double *ztime,
                          double *zrate) {

    int i, j, si, ei, iter;
    double dt, sumdf, sumdfp, Error, Epsilon, estimated;
    double *df, *dfp;
    df = new double[nztime];
    dfp = new double[nztime];
    dt = 1.0 / coupon_frequency;
    Epsilon = 1.0e-15;
    zrate[0] = yrate[0];

    df[0] = std::pow((1.0 + zrate[0] * dt), -ztime[0] * coupon_frequency);
    dfp[0] = -ztime[0] * coupon_frequency * std::pow((1.0 + zrate[0] * dt), -ztime[0] * coupon_frequency - 1.0);
    si = 1;

    for (i = 1; i < nytime; ++i) {
        iter = 0;
        Error = 0.0;
        ei = static_cast<int>(ytime[i] * coupon_frequency);
        zrate[ei - 1] = yrate[i];
        do {
            for (j = si; j < ei - 1; ++j) {
                zrate[j] =
                        (zrate[ei - 1] - zrate[si - 1]) / (ztime[ei - 1] - ztime[si - 1]) *
                        (ztime[j] - ztime[si - 1]) + zrate[si - 1];
            }
            for (j = si; j < ei; ++j) {
                df[j] = std::pow((1.0 + zrate[j] * dt), -ztime[j] * coupon_frequency);
                dfp[j] = -ztime[j] * coupon_frequency *
                         std::pow((1.0 + zrate[j] * dt), -ztime[j] * coupon_frequency - 1.0);
            }

            sumdf = sumdfp = 0.0;

            for (j = 0; j < ei; ++j) {
                sumdf += yrate[i] * df[j] * dt;
                sumdfp += ytime[i] * dfp[j] * dt;
            }

            sumdf += df[ei - 1];
            sumdfp += dfp[ei - 1];
            estimated = zrate[ei - 1] - (sumdf - 1.0) / sumdfp;
            Error = std::fabs(estimated - zrate[ei - 1]);
            zrate[ei - 1] = estimated;

            iter++;
        } while (Error > Epsilon && iter < 10000);

        for (j = si; j < ei - 1; ++j) {
            zrate[j] =
                    (zrate[ei - 1] - zrate[si - 1]) / (ztime[ei - 1] - ztime[si - 1]) *
                    (ztime[j] - ztime[si - 1]) + zrate[si - 1];
        }
        for (j = si; j < ei; ++j) {
            df[j] = std::pow((1.0 + zrate[j] * dt), -ztime[j] * coupon_frequency);
            dfp[j] = -ztime[j] * coupon_frequency *
                     std::pow((1.0 + zrate[j] * dt), -ztime[j] * coupon_frequency - 1.0);
        }
        si = ei;
    }

    delete[]df;
    delete[]dfp;
}

