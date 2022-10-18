#include "mainNet.h"

void make_linear_interpolation();

int main() {
//    make_gaussian_elimination();
//    printf("\n");
//    make_tridiagonal_elimination();
//    printf("\n");
//    choleski();
//    make_bisection_method_and_newton();
//    printf("\n");
//    make_gauss_newton();

//    make_levenberg_marquardt();
    make_linear_interpolation();
}

void make_linear_interpolation() {
    int i, n_days;
    int *maturity_day; // 이자율 잔존만기일
    double *zero_rate;
    n_days = 19;
    maturity_day = new int[n_days];
    zero_rate = new double[n_days];
    i = 0;

    maturity_day[i] = 1;
    zero_rate[i] = 2.542;

    maturity_day[++i] = 91;
    zero_rate[i] = 2.827;
    maturity_day[++i] = 182;
    zero_rate[i] = 2.989;
    maturity_day[++i] = 274;
    zero_rate[i] = 3.141;
    maturity_day[++i] = 366;
    zero_rate[i] = 3.275;
    maturity_day[++i] = 548;
    zero_rate[i] = 3.460;
    maturity_day[++i] = 732;
    zero_rate[i] = 3.615;
    maturity_day[++i] = 1099;
    zero_rate[i] = 3.808;

    maturity_day[++i] = 1463;
    zero_rate[i] = 3.963;
    maturity_day[++i] = 1827;
    zero_rate[i] = 4.097;
    maturity_day[++i] = 2193;
    zero_rate[i] = 4.191;
    maturity_day[++i] = 2558;
    zero_rate[i] = 4.263;
    maturity_day[++i] = 2923;
    zero_rate[i] = 4.331;
    maturity_day[++i] = 3290;
    zero_rate[i] = 4.395;
    maturity_day[++i] = 3654;
    zero_rate[i] = 4.460;
    maturity_day[++i] = 4019;
    zero_rate[i] = 4.494;
    maturity_day[++i] = 4384;
    zero_rate[i] = 4.523;
    maturity_day[++i] = 5481;
    zero_rate[i] = 4.587;
    maturity_day[++i] = 7308;
    zero_rate[i] = 4.630;

    int maturity;
    double zeros;
    maturity = 915;
    zeros = linear_interpolation(n_days, maturity_day, zero_rate, maturity);
    std::cout << "zero rate의 선형 보간 결과물: " << zeros << std::endl;
    delete[] maturity_day;
    delete[] zero_rate;

}

