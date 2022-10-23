//
// Created by jewoo on 2022-10-23.
//
#include <iostream>
#include "mathlib.h"
#include "index.h"
#include "yield.h"
#include "product.h"
#include "price.h"
#include "makeOption.h"

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
