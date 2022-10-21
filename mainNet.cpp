//
// Created by jewoo on 2022-10-18.
//

#include "mainNet.h"
#include "data.h"


void make_european_calloption() {
    CIndex index(110., 0.25, 0.03); // 기초 자산, 변동성, 배당률
    CYield yield(0.055); // 무위험이자율
    COption option(105.0, 1.0); // 행사가격, 잔존 만기


    option.european_calloption_price(index, yield);

    std::cout << "European Call option 가격: " << option.m_optionprice << std::endl;
}

void make_gaussian_elimination() {
    // 4원 1차 연립 방정식 해 찾기
    int i, neq;
    double **smatrix;
    double *known_value;
    double *unknown_value;

    neq = 4;
    smatrix = new double *[neq];
    for (i = 0; i < neq; ++i) {
        smatrix[i] = new double[neq];
    }
    known_value = new double[neq];
    unknown_value = new double[neq];

    smatrix[0][0] = 1;
    smatrix[0][1] = 1;
    smatrix[0][2] = 0;
    smatrix[0][3] = 3;

    smatrix[1][0] = 2;
    smatrix[1][1] = 1;
    smatrix[1][2] = -1;
    smatrix[1][3] = 1;

    smatrix[2][0] = 3;
    smatrix[2][1] = -1;
    smatrix[2][2] = -1;
    smatrix[2][3] = 2;

    smatrix[3][0] = -1;
    smatrix[3][1] = 2;
    smatrix[3][2] = 3;
    smatrix[3][3] = 1;

    known_value[0] = 4;
    known_value[1] = 1;
    known_value[2] = -3;
    known_value[3] = 4;

    gaussian_elimination(smatrix, known_value, unknown_value, neq);

    for (i = 0; i < neq; ++i) {
        std::cout << "x_" << i + 1 << " : " << unknown_value[i] << std::endl;
    }

    for (i = 0; i < neq; ++i) {
        delete smatrix[i];
    }
    delete[] smatrix, known_value, unknown_value;
}

void make_tridiagonal_elimination() {
    // 4원 1차 연립 방정식 해 찾기
    int i, neq;
    double **smatrix;
    double *known_value;
    double *unknown_value;

    neq = 4;

    smatrix = new double *[neq];
    for (i = 0; i < neq; ++i) {
        smatrix[i] = new double[neq];
    }
    known_value = new double[neq];
    unknown_value = new double[neq];

    smatrix[0][0] = 2;
    smatrix[0][1] = -1;
    smatrix[0][2] = 0;
    smatrix[0][3] = 0;

    smatrix[1][0] = -1;
    smatrix[1][1] = 2;
    smatrix[1][2] = -1;
    smatrix[1][3] = 0;

    smatrix[2][0] = 0;
    smatrix[2][1] = -1;
    smatrix[2][2] = 2;
    smatrix[2][3] = -1;

    smatrix[3][0] = 0;
    smatrix[3][1] = 0;
    smatrix[3][2] = -1;
    smatrix[3][3] = 2;

    known_value[0] = 1;
    known_value[1] = 0;
    known_value[2] = 0;
    known_value[3] = 1;

    tridiagonal_elimination(smatrix, known_value, unknown_value, neq);

    for (i = 0; i < neq; ++i) {
        std::cout << "x_" << i + 1 << " : " << unknown_value[i] << std::endl;
    }

    for (i = 0; i < neq; ++i) {
        delete smatrix[i];
    }
    delete[] smatrix, known_value, unknown_value;
}

void choleski() {

    int ncol;
    double **A; // 행렬 선언
    double *UV; // 미지수, A x = b 에서 x
    double *KV; // 기지수, A x = b 에서 b

    ncol = 3;
    A = new double *[ncol];
    for (int i = 0; i < ncol; ++i) {
        new double[ncol];
    }
    UV = new double[ncol];
    KV = new double[ncol];

    A[0][0] = 4.0;
    A[0][1] = -4.0;
    A[0][2] = 8.0;

    A[1][0] = -4.0;
    A[1][1] = 5.0;
    A[1][2] = -11.0;

    A[2][0] = 1.0;
    A[2][1] = -11.0;
    A[2][2] = 32.0;

    KV[0] = 1.0;
    KV[1] = 1.5;
    KV[2] = 2;

    cholesky_solver(A, UV, KV, ncol);
    gaussian_elimination(A, KV, UV, ncol);
    std::cout << "가우스 소거법을 이용한 연립방정식의 해는 " << std::endl;
    for (int i = 0; i < ncol; ++i) {
        std::cout << "x_" << i + 1 << " : " << UV[i] << std::endl;
    }
    for (int i = 0; i < ncol; ++i) {
        delete[]A[i];
    }

    delete[] A;
    delete[] KV;
    delete[] UV;
}

void make_bisection_method_and_newton() {
    double solution;
    std::cout << std::setprecision(30); // 결과 값의 소수점 단위를 확대 표시
    solution = bisection_method();
    std::cout << "bisection method 계산된 해 : " << solution << std::endl;
    solution = newtonraphson_method();
    std::cout << "newtonraphson method 계산된 해 : " << solution << std::endl;
}

void make_gauss_newton() {
    int i, n_data;
    double *mx, *my; // 데이터

    int n_parameter; // 미지 파라미터
    double *parameter; // 파라미터

    std::cout << std::setprecision(5);
    n_data = 10;
    mx = new double[n_data];
    my = new double[n_data];

    for (i = 0; i < n_data; ++i) {
        mx[i] = 1.0 + 0.2 * i;
        my[i] = 1.561 * mx[i] - 0.725 * mx[i] * mx[i] - 0.532 * std::exp(-0.288 * mx[i]);
    }

    n_parameter = 4;
    parameter = new double[n_parameter];
    for (i = 0; i < n_parameter; ++i) {
        parameter[i] = 1.0; // 파라미터 초기화
    }

    gauss_newton_parameter_solver(n_parameter, mx, my, n_parameter, parameter, fx_function);

    // 비선형 파라미터 추정
    for (i = 0; i < n_parameter; ++i) {
        std::cout << "p[" << i << "] : " << parameter[i] << std::endl;
    }

    delete[] parameter;
    delete[] mx;
    delete[] my;
}

void make_levenberg_marquardt() {

    // 비선형 파라미터 추정
    int i;
    int n_strike, n_maturity; // 행사 가격 수량, 잔존 만기 수량
    double risk_free_rate, dividend, S, carry_cost; // 무위험 이자율, 배당율, 기초자산, 캐리 비용

    double *v_strike, *v_maturity; // 행사가격, 잔존만기
    double **m_vol; // 마켓 볼륨

    int n_parameter; // 미지 파라미터
    double *parameter; // 파라미터

    // display points
    std::cout << std::setprecision(20);

    n_strike = 5;
    n_maturity = 6;

    v_strike = new double[n_strike];
    v_maturity = new double[n_maturity];
    m_vol = new double *[n_strike];

    for (i = 0; i < n_strike; ++i) {
        m_vol[i] = new double[n_maturity];
    }

    read_market_vol("marketvol.txt", n_strike, v_strike, n_maturity, v_maturity, m_vol);

    S = v_strike[2];
    risk_free_rate = 0.03;
    dividend = 0.0;
    carry_cost = risk_free_rate - dividend;

    n_parameter = 7;
    parameter = new double[n_parameter];

    for (i = 0; i < n_parameter; ++i) {
        parameter[i] = 0.1; // 파라미터 초기화
    }

    levenberg_marquardt_parameter_solver(S, carry_cost, n_strike, v_strike, n_maturity,
                                         v_maturity, m_vol, n_parameter, parameter,
                                         implied_volatility_function);
    // 비선형 파라미터 추정

    for (i = 0; i < n_parameter; ++i) {
        std::cout << "p[" << i << "] : " << parameter[i] << std::endl;
    }

    double **imvol;
    imvol = new double *[n_strike]; // implied vol function에 의한 내재변동성
    for (i = 0; i < n_strike; ++i) {
        imvol[i] = new double[n_maturity];
    }

    implied_volatility(S, carry_cost, n_strike, v_strike, n_maturity, v_maturity,
                       imvol,
                       n_parameter, parameter);

    save_vol("impliedvol.txt", n_strike, v_strike, n_maturity, v_maturity, imvol);

    for (i = 0; i < n_strike; ++i) {
        delete[]m_vol[i];
        delete[]imvol[i];
    }

    delete[]parameter;
    delete[]v_strike;
    delete[] v_maturity;
    delete[] imvol;
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

void make_cubicspline_interpolation() {
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
    zeros = cubicspline_interpolation(n_days, maturity_day, zero_rate, maturity);
    std::cout << "zero rate의 스팔라인 보간 결과물: " << zeros << std::endl;
    delete[] maturity_day;
    delete[] zero_rate;
}


