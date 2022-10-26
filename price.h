//
// Created by jewoo on 2022-10-23.
//

#pragma once

class CIndex;

class CProduct;

class CYield;

class CPrice {
public:
    CPrice();

    ~CPrice();

    void black_scholes_option_price(CIndex &index,
                                    CYield &yield,
                                    CProduct &product);

    void simulation_european_option_price(CIndex &index,
                                          CYield &yield,
                                          CProduct &product, int n_sim);

    void binomial_tree_european_option_price(CIndex &index,
                                             CYield &yield,
                                             CProduct &product, int n_step);

    void lognormal_binomial_tree_european_option_price(CIndex &index,
                                                       CYield &yield,
                                                       CProduct &product, int n_step);

    void trinomial_tree_european_option_price(CIndex &index,
                                              CYield &yield,
                                              CProduct &product, int n_step);

    void implicit_fdm_european_option_price(CIndex &index,
                                            CYield &yield,
                                            CProduct &product, int n_step, int n_spot);

    // 스프레드 옵션
    void spread_option_price(CIndex &index1, CIndex &index2, CYield &yield,
                             CProduct &product);

    void simulation_spread_option_price(CIndex &index1, CIndex &index2, CYield &yield,
                                        CProduct &product, int n_sim);

    void binomial_tree_spread_option_price(CIndex &index1, CIndex &index2, CYield &yield,
                                           CProduct &product, int n_step);

    void trinomial_tree_spread_option_price(CIndex &index1, CIndex &index2, CYield &yield,
                                            CProduct &product, int n_step);

    void implicit_fdm_spread_option_price(CIndex &index1, CIndex &index2, CYield &yield,
                                          CProduct &product, int n_step, int n_spot);


    // 분산기법
    void simulation_european_option_price(CIndex &index,
                                          CYield &yield,
                                          CProduct &product, unsigned long n_sim,
                                          double *rn);

    void antithetic_variation_simulation_european_option_space(CIndex &index,
                                                               CYield &yield,
                                                               CProduct &product,
                                                               unsigned long n_sim,
                                                               double *rn);

public:
    double m_price;


};