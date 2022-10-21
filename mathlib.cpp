//
// Created by jewoo on 2022-10-13.
//

#include "mathlib.h"
#include <iostream>
#include <ctgmath>

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
