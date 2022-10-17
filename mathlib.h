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