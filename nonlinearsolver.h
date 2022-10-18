//
// Created by jewoo on 2022-10-15.
//

#pragma once
#define _CRT_SECURE_NO_DEPRECATE

void gauss_newton_parameter_solver(int n_ref,
                                   double *ref_x, double *ref_y,
                                   int n_parameter, double *parameter,
                                   double object_func(int n_parameter, double *parameter, double x));

double fx_function(int n_parameter, double *parameter, double x);

// =====================================================================================================================
void levenberg_marquardt_parameter_solver(double S, double cc, int n_x,
                                          double *v_x, int n_y, double *v_y,
                                          double **m_ref, int n_parameter,
                                          double *parameter,
                                          double object_func(int n_parameter, double *parameter, double x, double t));

double implied_volatility_function(int n_parameter, double *parameter, double x, double t);

void implied_volatility(double S, double cc, int n_x, double *v_x,
                        int n_y, double *v_y, double **mat, int n_parameter, double *parameter);