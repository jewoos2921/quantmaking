//
// Created by jewoo on 2022-10-28.
//

#pragma once

class CIndex;

class CYield;

class CProduct;

class CPriceAlpha {
public:
    CPriceAlpha();

    ~CPriceAlpha();

    void black_scholes_option_price_greeks(CIndex &index,
                                           CYield &yield,
                                           CProduct &product);

    void simulation_european_option_price_greeks(CIndex &index,
                                                 CYield &yield,
                                                 CProduct &product, int n_sim);

    void binomial_tree_european_option_price_greeks(CIndex &index,
                                                    CYield &yield,
                                                    CProduct &product, int n_step);


    void trinomial_tree_european_option_price_greeks(CIndex &index,
                                                     CYield &yield,
                                                     CProduct &product, int n_step);

    void implicit_fdm_european_option_price_greeks(CIndex &index,
                                                   CYield &yield,
                                                   CProduct &product, int n_step, int n_spot);


    double black_scholes_option_price(int iop, double S, double X, double r, double q,
                                      double sigma, double T);

    void black_scholes_option_greeks(int iop, double S, double X, double r, double q,
                                     double sigma, double T);


    double ci_binomial_tree_european_option_space(int iop, double S, double X, double r, double q,
                                                  double sigma, double T, int n_step, double alpha);

    double ci_trinomial_tree_european_option_space(int iop, double S, double X, double r, double q,
                                                   double sigma, double T, int n_step, double alpha);

    double ci_implicit_fdm_european_option_space(int iop, double S, double X, double r, double q,
                                                 double sigma, double T, int n_step, int n_spot, double alpha);

public:
    double m_price; // 평가 값
    double m_alpha; // confidence coefficient

};

