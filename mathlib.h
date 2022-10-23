//
// Created by jewoo on 2022-10-13.
//

#pragma once
#define _CRT_SECURE_NO_DEPRECATE

void gaussian_elimination(double **Smatrix,
                          double *Known,
                          double *Unknown,
                          int n_eqns);

void tridiagonal_elimination(double **Smatrix,
                             double *Known,
                             double *Unknown,
                             int n_eqns);

void cholesky_decomposition(double **L, double **A, int n_col);

void cholesky_elimination(double **L, double *UV, double *KV, int n_col);

void cholesky_solver(double **A, double *UV, double *KV, int n_col);

double nonlinear_function(double X);

double bisection_method();

double newtonraphson_method();

double linear_interpolation(int n_data, int *x_data, double *y_data, int x);

void set_system_matrix(int n_data, double **Smatrix, double *known, int *x_data,
                       double *y_data);

double cubicspline_interpolation(int n_data,
                                 int *x_data,
                                 double *y_data, int x);


double normdistrand(); // 정규분포 난수 생성

void normal_distribution_goodness_fit_test(unsigned long nrand, double *randnum); // 정규뷴포 적합성 검증


#define EPS 12e-7

double inverse_normal_cumulative_distribution_function(double p); // 정규분포 난수 생성


#define PI 3.14159265368979323846264338327950288419716939937519582097

#define MAX(a, b) (((a)>(b))?(a):(b))

double normdistrand_BoxMuller();  // 정규분포 난수 생성

double N(double z);


