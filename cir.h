//
// Created by jewoo on 2022-11-02.
//

#pragma once

#pragma once

class CCIR {
public:
    CCIR();

    ~CCIR();

    CCIR(double alpha, double theta, double sigma);

public:
    double m_alpha; // 평균회귀 계수
    double m_theta; // 이자율 장기 균형점
    double m_sigma; // 변동성
    double m_ini_shortrate{}; // 초기 단기 이자율

public:
    double zero_bond_value(double shortrate, double start_t, double end_t);

    double spot_rate(double shortrate, double start_t, double end_t);

    double zero_bond_option_value(int iop, double shortrate, double start_t,
                                  double option_maturity, double bond_maturity, double strike);

    double zero_bond_option_value(int iop, double notional, double shortrate,
                                  double start_t,
                                  double option_maturity, double bond_maturity, double strike);

    double get_rate_k(double notional, double coupon, double coupon_period,
                      int ncoupon, double shortrate, double option_maturity, double strike);

    double fixed_coupon_bond_option_value(int iop, double notional, double coupon_rate,
                                          int coupon_frequency,
                                          double shortrate, double start_t,
                                          double option_maturity, double bond_maturity, double strike);
};