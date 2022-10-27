//
// Created by jewoo on 2022-10-23.
//
#include <iostream>
#include <ctime>
#include "mathlib.h"
#include "index.h"
#include "yield.h"
#include "product.h"
#include "price.h"
#include "makeOption.h"
#include "sobol.h"

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
    n_sim = 30000; // 시뮬레이션 횟수
    n_spot = 400; // 기초자산의 가격을 n_stpo 수로 나누어 계산

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

