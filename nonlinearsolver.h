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