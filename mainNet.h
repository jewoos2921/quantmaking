//
// Created by jewoo on 2022-10-18.
//

#pragma once
#define _CRT_SECURE_NO_DEPRECATE

#include <iostream>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <random>
#include "index.h"
#include "yield.h"
#include "option.h"
#include "mathlib.h"
#include "nonlinearsolver.h"

void make_gaussian_elimination();

void make_tridiagonal_elimination();

void choleski();

void make_bisection_method_and_newton();

void make_gauss_newton();

void make_levenberg_marquardt();

void make_linear_interpolation();

void make_cubicspline_interpolation();

void make_normal_distribution_goodness_fit_test();

void make_inverse_normal_cumulative_distribution_function();