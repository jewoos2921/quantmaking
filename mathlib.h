//
// Created by jewoo on 2022-10-13.
//

#pragma once
#define CRT_SECURE_NO_DEPRECATE

void gaussian_elimination(double **Smatrix,
                          double *Known,
                          double *Unknown,
                          int n_eqns);

void tridiagonal_elimination(double **Smatrix,
                             double *Known,
                             double *Unknown,
                             int n_eqns);

void cholesky_decomposition(double **L, double **A, int n_col); // A: in, L: out, col: 배열 크기

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

// 분산 검증
void mean_stddev_error(unsigned long nrand, double *value);

// 준난수 생성
double *halton_sequence(unsigned num, unsigned prime_number);


double halton_sequence_(int num, int prime_number);

double n(double z);

// 이분법을 이용한 옵션의 내재변동성을 찾기
double implied_volatility_bisection(double spot, double strike, double riskfree,
                                    double dividend, double maturity, double option_price);

double european_calloption_price(double spot, double strike, double riskfree,
                                 double dividend, double vol, double maturity);


double european_calloption_vega(double spot, double strike, double riskfree,
                                double dividend, double volatility, double maturity);

double implied_volatility_newtonraphson(double spot, double strike, double riskfree,
                                        double dividend, double maturity, double option_price);

// 이자율 변환
void par_to_zero_bootstrapping(int nytime, double *ytime, double *yrate,
                               int coupon_frequency, int nztime, double *ztime, double *zrate);

double linear_interpolation(int n_time, double *time, double *rate, double t);


double caplet_price(int option_type, double forward, double strike, double tau, double maturity, double vol,
                    double discount_factor);

double cap_floor_price(int option_type, double forwardrate, double strike, double vol, int ntime, double *zerorate,
                       double *markettime, double basis);