//
// Created by jewoo on 2022-10-23.
//
#include <iostream>
#include <ctime>
#include <iomanip>
#include <fstream>
#include "mathlib.h"
#include "index.h"
#include "yield.h"
#include "product.h"
#include "price.h"
#include "makeOption.h"
#include "sobol.h"
#include "priceAlpha.h"
#include "nonlinearsolver.h"
#include "data.h"


void make_option() {
    double S, X, r, q, sigma, T;

    int n_step, n_spot;
    unsigned long n_sim;


    S = 100.0; // 기초 자산 가격
    X = 110.0; // 행사 가격
    r = 0.05; // risk free rate
    q = 0.01; // dividend rate
    sigma = 0.3; // 변동성
    T = 1.0; // 잔존만기

    n_step = 200; // 잔존만기에 대해 n_step으로 계산, tree와 fdm에서
    n_spot = 400; // 기초자산의 가격을 n_stpo 수로 나누어 계산
    n_sim = 30000; // 시뮬레이션 횟수

    CIndex index(S, sigma, q);
    CYield yield(r);
    CProduct eoption("call", X, T);
    CPrice price;

    price.black_scholes_option_price(index, yield, eoption);
    std::cout << "해석해에 의한 옵션가격: " << price.m_price << std::endl;

    price.simulation_european_option_price(index, yield, eoption, n_sim);
    std::cout << "시뮬레이션에 의한 옵션가격: " << price.m_price << std::endl;

    price.binomial_tree_european_option_price(index, yield, eoption, n_step);
    std::cout << "Binomial Tree에 의한 옵션가격: " << price.m_price << std::endl;

    price.lognormal_binomial_tree_european_option_price(index, yield, eoption, n_step);
    std::cout << "로그정규분포 적용 Binomial Tree에 의한 옵션가격: " << price.m_price << std::endl;

    price.trinomial_tree_european_option_price(index, yield, eoption, n_step);
    std::cout << "Trinomial Tree에 의한 옵션가격: " << price.m_price << std::endl;

    price.implicit_fdm_european_option_price(index, yield, eoption, n_step, n_spot);
    std::cout << "Implicit FDM에 의한 옵션가격: " << price.m_price << std::endl;

}


void make_option2() {
    double S1, S2, X, r, q1, q2, sigma1, sigma2, T, rho;

    int n_step, n_spot;
    unsigned long n_sim;


    S1 = 110.0; // 기초 자산1 가격
    S2 = 100.0; // 기초 자산2 가격
    X = 110.0; // 행사 가격
    r = 0.05; // risk free rate
    q1 = 0.01; // 기초 자산1 dividend rate
    q2 = 0.015; // 기초 자산2 dividend rate

    sigma1 = 0.2; // 기초 자산1 변동성
    sigma2 = 0.3; // 기초 자산2 변동성
    rho = 0.5; // 기초 자산1과 기초 자산2의 상관계수
    T = 1.0; // 잔존만기

    n_step = 200; // 잔존만기에 대해 n_step으로 계산, simulation, tree와 fdm에서
    n_sim = 30000; // 시뮬레이션 횟수
    n_spot = 200; // 기초자산의 가격을 n_spot 수로 나누어 계산

    CIndex index1(S1, sigma1, q1);
    CIndex index2(S2, sigma2, q2);

    CYield yield(r);

    CProduct soption("call", X, T);

    soption.m_correlation = rho;

    CPrice price;

    price.spread_option_price(index1, index2, yield, soption);
    std::cout << "해석해에 의한 옵션가격: " << price.m_price << std::endl;

    price.simulation_spread_option_price(index1, index2, yield, soption, n_sim);
    std::cout << "시뮬레이션에 의한 옵션가격: " << price.m_price << std::endl;

    price.binomial_tree_spread_option_price(index1, index2, yield, soption, n_step);
    std::cout << "Binomial Tree에 의한 옵션가격: " << price.m_price << std::endl;

    price.trinomial_tree_spread_option_price(index1, index2, yield, soption, n_step);
    std::cout << "Trinomial Tree에 의한 옵션가격: " << price.m_price << std::endl;

    price.implicit_fdm_spread_option_price(index1, index2, yield, soption, n_step, n_spot);
    std::cout << "Implicit FDM에 의한 옵션가격: " << price.m_price << std::endl;

}


void make_option3() {
    double S, X, r, q, sigma, T;

    unsigned long n_sim;
    double *rn;


    S = 100.0; // 기초 자산 가격
    X = 110.0; // 행사 가격
    r = 0.05; // risk free rate
    q = 0.01; // dividend rate
    sigma = 0.3; // 변동성
    T = 1.0; // 잔존만기
    n_sim = 1000000; // 시뮬레이션 횟수

    rn = new double[n_sim]; // random number


    CIndex index(S, sigma, q);
    CYield yield(r);
    CProduct eoption("call", X, T);
    CPrice price;

    printf("\n");
    price.black_scholes_option_price(index, yield, eoption);
    std::cout << "해석해에 의한 옵션가격: " << price.m_price << std::endl;

    std::srand(static_cast<unsigned >(std::time(nullptr))); // Random seed 초기화

    printf("\n");
    price.simulation_european_option_price(index, yield, eoption, n_sim, rn);
    std::cout << "시뮬레이션에 의한 옵션가격: " << price.m_price << std::endl;
    printf("\n");
    price.antithetic_variation_simulation_european_option_space(index, yield, eoption, n_sim, rn);
    std::cout << "대조변수기법이 적용된 시물레이션에 의한 옵션가격: " << price.m_price << std::endl;

    printf("\n");
    price.control_variation_simulation_european_option_space(index, yield, eoption, n_sim, rn);
    std::cout << "통제변수기법이 적용된 시물레이션에 의한 옵션가격: " << price.m_price << std::endl;


    printf("\n");
    price.importance_sampling_simulation_european_option_space(index, yield, eoption, n_sim, rn);
    std::cout << "중요표본추출법이 적용된 시물레이션에 의한 옵션가격: " << price.m_price << std::endl;
    delete[] rn;

}

void make_option4() {
    unsigned i, ns;
    double *rn;
    ns = 1000000;
    rn = new double[ns];

    double *haltons = halton_sequence(ns, 2); // halton sequence를 이용한 난수 생성

    // 평균0, 분산1 인 정규분포형 난수 생성
    for (i = 0; i < ns; ++i) {
        rn[i] = inverse_normal_cumulative_distribution_function(haltons[i]);
    }

    normal_distribution_goodness_fit_test(ns, rn); // 정규 분포 적합성 검증

    delete[] rn;
    delete[] haltons;
}

void make_option5() {
    unsigned i, ns;
    double *rn;
    ns = 1000000;
    rn = new double[ns];
    double *sobol = sobol_points(ns);
    for (i = 0; i < ns; ++i) {
        rn[i] = inverse_normal_cumulative_distribution_function(sobol[i]);
    }
    normal_distribution_goodness_fit_test(ns, rn);
    delete[]sobol;
    delete[]rn;
}


void make_option6() {
    double S, X, r, q, sigma, T;

    unsigned long n_sim;


    S = 100.0; // 기초 자산 가격
    X = 110.0; // 행사 가격
    r = 0.05; // risk free rate
    q = 0.01; // dividend rate
    sigma = 0.3; // 변동성
    T = 1.0; // 잔존만기
    n_sim = 1000000; // 시뮬레이션 횟수




    CIndex index(S, sigma, q);
    CYield yield(r);
    CProduct eoption("call", X, T);
    CPrice price;

    printf("\n");
    price.black_scholes_option_price(index, yield, eoption);
    std::cout << "해석해에 의한 옵션가격: " << price.m_price << std::endl;


    printf("\n");
    price.halton_sequence_simulation_european_option_space(index, yield, eoption, n_sim);
    std::cout << "시뮬레이션에 의한 옵션가격: " << price.m_price << std::endl;
    printf("\n");

//    price.sobol_sequence_simulation_european_option_space(index, yield, eoption, n_sim);
//    std::cout << "대조변수기법이 적용된 시물레이션에 의한 옵션가격: " << price.m_price << std::endl;


}

void make_option7() {
    double S, X, r, q, sigma, T, alpha;

    int n_step, n_spot;


    S = 100.0; // 기초 자산 가격
    X = 110.0; // 행사 가격
    r = 0.05; // risk free rate
    q = 0.01; // dividend rate
    sigma = 0.3; // 변동성
    T = 1.0; // 잔존만기

    n_step = 200; // 잔존만기에 대해 n_step으로 계산, tree와 fdm에서
    n_spot = 400; // 기초자산의 가격을 n_stpo 수로 나누어 계산
    alpha = 5.0; // confidence Coefficient
    // alpha 1은 = 68.26%
    // alpha 1.65은 = 90.1%
    // alpha 2은 = 95.44%
    // alpha 2.58은 = 99%
    // alpha 3은 = 99.75%
    // alpha 4은 = 99.995%
    // alpha 5은 = 99.99995%
    // 의 신뢰수준

    CIndex index(S, sigma, q);
    CYield yield(r);
    CProduct eoption("call", X, T);
    CPrice price;

    price.black_scholes_option_price(index, yield, eoption);
    std::cout << "해석해에 의한 옵션가격: " << price.m_price << std::endl;

    price.ci_binomial_tree_european_option_space(index, yield, eoption, n_step, alpha);
    std::cout << "Binomial Tree에 의한 옵션가격: " << price.m_price << std::endl;

    price.ci_trinomial_tree_european_option_space(index, yield, eoption, n_step, alpha);
    std::cout << "Trinomial Tree에 의한 옵션가격: " << price.m_price << std::endl;

    price.ci_implicit_fdm_european_option_space(index, yield, eoption, n_step, n_spot, alpha);
    std::cout << "Implicit FDM에 의한 옵션가격: " << price.m_price << std::endl;

}

void make_option_greeks() {
    double S, X, r, q, sigma, T, alpha;

    int n_step, n_spot;
    unsigned long n_sim;


    S = 100.0; // 기초 자산 가격
    X = 110.0; // 행사 가격
    r = 0.05; // risk free rate
    q = 0.01; // dividend rate
    sigma = 0.3; // 변동성
    T = 1.0; // 잔존만기

    n_step = 200; // 잔존만기에 대해 n_step으로 계산, tree와 fdm에서
    n_sim = 30000; // 시뮬레이션 횟수
    n_spot = 400; // 기초자산의 가격을 n_stpo 수로 나누어 계산
    alpha = 5.0; // confidence Coefficient
    // alpha 1은 = 68.26%
    // alpha 1.65은 = 90.1%
    // alpha 2은 = 95.44%
    // alpha 2.58은 = 99%
    // alpha 3은 = 99.75%
    // alpha 4은 = 99.995%
    // alpha 5은 = 99.99995%
    // 의 신뢰수준

    CIndex index(S, sigma, q);
    CYield yield(r);
    CProduct eoption("call", X, T);
    CPriceAlpha price;

    price.m_alpha = alpha;

    price.black_scholes_option_price_greeks(index, yield, eoption);
    std::cout << "해석해에 의한 옵션가격: " << price.m_price << std::endl;

    price.simulation_european_option_price_greeks(index, yield, eoption, n_sim);
    std::cout << "시뮬레이션에 의한 옵션가격: " << price.m_price << std::endl;

    price.binomial_tree_european_option_price_greeks(index, yield, eoption, n_step);
    std::cout << "Binomial Tree에 의한 옵션가격: " << price.m_price << std::endl;

    price.trinomial_tree_european_option_price_greeks(index, yield, eoption, n_step);
    std::cout << "Trinomial Tree에 의한 옵션가격: " << price.m_price << std::endl;

    price.implicit_fdm_european_option_price_greeks(index, yield, eoption, n_step, n_spot);
    std::cout << "Implicit FDM에 의한 옵션가격: " << price.m_price << std::endl;

}


void make_volatility() {
    double spot, strike, ristfree, dividend, volaitirity, maturity;

    spot = 100.0;
    strike = 105.0;
    ristfree = 0.055;
    dividend = 0.02;
    volaitirity = 0.25;
    maturity = 1.0;

    double option_price = european_calloption_price(spot, strike, ristfree, dividend,
                                                    volaitirity, maturity);
    std::cout << "European Call option 가격 : " << option_price << std::endl;
    std::cout << "기초자산의 입력 변동성 : " << volaitirity << std::endl;

    double implied_vol = implied_volatility_newtonraphson(spot, strike, ristfree, dividend,
                                                          volaitirity, maturity);

    std::cout << "newton raphson method 계산된 내재변동성: " << implied_vol << std::endl;
}

void make_implied_volatility() {

    int i;
    int n_strike, n_maturity; // 행사가격, 잔존만기 수량
    double riskfreerate, dividend, S, carrycost; // 무위험이자율, 배당률, 기초자산, 캐리비용

    double *v_strike, *v_maturity; // 행사가격, 잔존만기
    double **m_vol; // market vol

    int n_parameter;  // 미지 파라미터
    double *parameter; // 파라미터

    // display points
    std::cout << std::setprecision(20);

    n_strike = 5;
    n_maturity = 6;

    v_strike = new double[n_strike];
    v_maturity = new double[n_maturity];
    m_vol = new double *[n_strike];

    for (i = 0; i < n_strike; ++i) { m_vol[i] = new double[n_maturity]; }

    read_market_vol("marketvol.txt", n_strike, v_strike, n_maturity, v_maturity, m_vol);

    S = v_strike[(n_strike - 1) / 2];
    riskfreerate = 0.03;
    dividend = 0.00;
    carrycost = riskfreerate - dividend;

    n_parameter = 7;
    parameter = new double[n_parameter];
    for (i = 0; i < n_parameter; ++i) { parameter[i] = 0.1; } // 파라미터 초기화

    nonlinear_parameter_solver(S, carrycost, n_strike, v_strike, n_maturity,
                               v_maturity, m_vol, n_parameter, parameter,
                               implied_volatility_function); // 비선형 파라미터 추정

    for (i = 0; i < n_parameter; ++i) { std::cout << "p[" << i << "]: " << parameter[i] << std::endl; }

    double **imvol, **localvol;
    imvol = new double *[n_strike]; // implied vol function에 의한 내재변동성
    localvol = new double *[n_strike]; // local vol

    for (i = 0; i < n_strike; ++i) {
        imvol[i] = new double[n_maturity];
        localvol[i] = new double[n_maturity];
    }

    implied_volatility_2(S, carrycost, n_strike, v_strike, n_maturity, v_maturity, imvol,
                         n_parameter, parameter);
    save_vol("impliedvol.txt", n_strike, v_strike, n_maturity, v_maturity, imvol);
    local_volatility(S, carrycost, n_strike, v_strike, n_maturity, v_maturity, localvol,
                     n_parameter, parameter);
    save_vol("localvol.txt", n_strike, v_strike, n_maturity, v_maturity, localvol);

    for (i = 0; i < n_strike; ++i) {
        delete[] m_vol[i];
        delete[] imvol[i];
        delete[] localvol[i];
    }

    delete[] parameter;
    delete[] v_strike;
    delete[] v_maturity;
    delete[] imvol;
    delete[] localvol;
}

void make_bootstrapping() {
    std::cout << std::setprecision(20);

    int i, nytime, coupon_frequency;
    double *ytime;
    double *yrate;
    nytime = 15;
    coupon_frequency = 4;

    ytime = new double[nytime];
    yrate = new double[nytime];


    ytime[0] = 0.25;
    yrate[0] = 0.05360;
    ytime[1] = 1.0;
    yrate[1] = 0.05030;
    ytime[2] = 2.0;
    yrate[2] = 0.04890;
    ytime[3] = 3.0;
    yrate[3] = 0.04780;
    ytime[4] = 4.0;
    yrate[4] = 0.04738;
    ytime[5] = 5.0;
    yrate[5] = 0.04758;
    ytime[6] = 6.0;
    yrate[6] = 0.04778;
    ytime[7] = 7.0;
    yrate[7] = 0.04785;
    ytime[8] = 8.0;
    yrate[8] = 0.04810;
    ytime[9] = 9.0;
    yrate[9] = 0.04830;
    ytime[10] = 10.0;
    yrate[10] = 0.04855;
    ytime[11] = 11.0;
    yrate[11] = 0.04880;
    ytime[12] = 12.0;
    yrate[12] = 0.04905;
    ytime[13] = 15.0;
    yrate[13] = 0.04930;
    ytime[14] = 20.0;
    yrate[14] = 0.04955;


    int nztime;
    double *ztime;
    double *zrate;

    nztime = static_cast<int>(ytime[nytime - 1] * coupon_frequency);
    ztime = new double[nztime];
    zrate = new double[nztime];

    for (i = 0; i < nztime; ++i) { ztime[i] = static_cast<double>(i + 1) / coupon_frequency; }

    par_to_zero_bootstrapping(nytime, ytime, yrate,
                              coupon_frequency, nztime, ztime,
                              zrate);

    std::ofstream outf("result.txt", std::ios::out | std::ios::trunc);
    // ios::trunc 사용시 기존 데이터 삭제후 추가
    // ios::app 사용시 데이터 추가 기능

    outf << std::setprecision(20);
    outf << "=입력된 액면가수익룰=" << std::endl;
    outf << "만기 액면가수익률" << std::endl;
    for (i = 0; i < nytime; ++i) {
        outf << ytime[i] << " " << yrate[i] << std::endl;
        outf << std::endl;
    }

    outf << "=계산된 현물이자율=" << std::endl;
    outf << "만기 현물이자율" << std::endl;
    for (i = 0; i < nztime; ++i) {
        outf << ztime[i] << " " << zrate[i] << std::endl;
        outf << std::endl;
    }

    delete[]ytime;
    delete[]yrate;
    delete[]ztime;
    delete[]zrate;
}