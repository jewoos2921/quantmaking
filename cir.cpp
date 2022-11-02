//
// Created by jewoo on 2022-11-02.
//

#include <iostream>
#include <cmath>
#include "mathlib.h"
#include "cir.h"

CCIR::CCIR() {
    m_alpha = 0.1;
    m_theta = 0.1;
    m_sigma = 0.01;
    m_ini_shortrate = 0.05;
}

CCIR::~CCIR() = default;

CCIR::CCIR(double alpha, double theta, double sigma) {
    m_alpha = alpha;
    m_theta = theta;
    m_sigma = sigma;
}


double CCIR::zero_bond_value(double shortrate, double start_t, double end_t) {
    double tau, sigma2, a2;
    double AtT, BtT;
    double gamma, denominator;

    tau = end_t - start_t;
    sigma2 = m_sigma * m_sigma;
    a2 = m_alpha * m_alpha;

    gamma = std::sqrt(a2 + 2 * sigma2);
    denominator = (gamma + m_alpha) * (std::exp(gamma * tau) - 1)
                  + 2.0 * gamma;

    BtT = 2.0 * (std::exp(gamma * tau) - 1.0) / denominator;
    AtT = std::pow(((2.0 * gamma * std::exp((gamma + m_alpha) * tau / 2.0)) / denominator),
                   (2.0 * m_alpha * m_theta / sigma2));

    return AtT * std::exp(-BtT * shortrate);
}

double CCIR::spot_rate(double shortrate, double start_t, double end_t) {
    double PtT = this->zero_bond_value(shortrate, start_t, end_t);
    return -std::log(PtT) / (end_t - start_t);
}

double CCIR::zero_bond_option_value(int iop, double shortrate, double start_t, double option_maturity,
                                    double bond_maturity, double strike) {
    return 0;
}

double
CCIR::zero_bond_option_value(int iop, double notional, double shortrate, double start_t, double option_maturity,
                             double bond_maturity, double strike) {


    double gamma, phi, psi, rstar;
    double sigma2, t1, t2;
    double PtT, PtS;
    double denominator;
    double ATS, BTS;
    double nu, z1, kappa1, z2, kappa2;

    sigma2 = m_sigma * m_sigma;
    gamma = std::sqrt(m_alpha * m_alpha + 2 * sigma2);

    t1 = option_maturity - start_t;
    t2 = bond_maturity - option_maturity;

    phi = 2.0 * gamma / (sigma2 * (std::exp(gamma * t1) - 1.0));
    psi = (m_alpha + gamma) / sigma2;

    denominator = (gamma * m_alpha) * (std::exp(gamma * t2) - 1) + 2.0 * gamma;

    BTS = 2.0 * (std::exp(gamma * t2) - 1.0) / denominator;
    ATS = std::pow(((2.0 * gamma * std::exp((gamma + m_alpha) * t2 / 2.0)) / denominator),
                   (2.0 * m_alpha * m_theta / sigma2));

    rstar = std::log(notional * ATS / strike) / BTS;

    nu = 4.0 * m_alpha * m_theta / sigma2;
    z1 = 2.0 * rstar * (phi + psi + BTS);
    z2 = 2.0 * rstar * (psi + phi);

    kappa1 = 2.0 * phi * phi * shortrate * std::exp(gamma * t1) / (psi + phi + BTS);
    kappa2 = 2.0 * phi * phi * shortrate * std::exp(gamma * t1) / (psi + phi);


    PtT = this->zero_bond_value(shortrate, start_t, option_maturity);
    PtS = this->zero_bond_value(shortrate, start_t, bond_maturity);

    if (iop == 1) {
        return notional * PtS * QS(z1, nu, kappa1) - strike * PtT * QS(z2,
                                                                       nu, kappa2);
    } else {
        return strike * PtT * (1.0 - QS(z2, nu, kappa2))
               - notional * PtS * (1.0 - QS(z1, nu, kappa1));
    }
}

double CCIR::get_rate_k(double notional, double coupon, double coupon_period, int ncoupon, double shortrate,
                        double option_maturity, double strike) {
    int i, iter;
    double r_k, fx, fxp, fxm, fxprime, payoff;
    double spread, error, error_tolerance;

    r_k = shortrate;
    iter = 0;
    spread = 0.00010;
    error_tolerance = 1.0e-15;

    do {
        fx = fxp = fxm = 0.0;
        for (i = 0; i < ncoupon; ++i) {
            if (i < ncoupon - 1) {
                payoff = coupon;
            } else {
                payoff = coupon + notional;
            }
            fx += this->zero_bond_value(r_k, option_maturity, option_maturity + coupon_period * (i + 1)) * payoff;
            fxp += this->zero_bond_value(r_k + spread, option_maturity,
                                         option_maturity + coupon_period * (i + 1)) * payoff;
            fxm += this->zero_bond_value(r_k - spread, option_maturity, option_maturity + coupon_period * (i + 1)) *
                   payoff;
        }

        fxprime = (fxp - fxm) / (2.0 * spread);
        r_k = r_k - (fx - strike) / fxprime;
        error = std::fabs(fx - strike);
        iter++;
    } while (error > error_tolerance && iter < 100);

    if (iter == 100) {
        std::cout << "허용오차범위는 불만족하나, 가장 근사한 값은 " << r_k << "입니다." << std::endl;
    } else {
        std::cout << "newton raphson method 계산횟수: " << iter << "입니다." << std::endl;
        std::cout << "rate_k의 값: " << r_k << std::endl;
    }

    return r_k;
}

double CCIR::fixed_coupon_bond_option_value(int iop, double notional, double coupon_rate, int coupon_frequency,
                                            double shortrate, double start_t, double option_maturity,
                                            double bond_maturity, double strike) {

    int i, ncoupon;
    double rate_k, option_value;
    double coupon, coupon_period, payoff, strike_i;

    coupon_period = 1.0 / coupon_frequency;
    coupon = notional * coupon_rate * coupon_period;
    ncoupon = static_cast<int>((bond_maturity - option_maturity) * coupon_frequency);
    rate_k = this->get_rate_k(notional, coupon, coupon_period, ncoupon, shortrate,
                              option_maturity, strike);

    option_value = 0.0;
    for (i = 0; i < ncoupon; ++i) {
        if (i < ncoupon - 1) {
            payoff = coupon;
        } else {
            payoff = coupon + notional;
        }
        strike_i = this->zero_bond_value(rate_k, option_maturity,
                                         option_maturity + coupon_period * (i + 1)) * payoff;
        option_value += this->zero_bond_option_value(iop, payoff, shortrate, start_t,
                                                     option_maturity, option_maturity
                                                                      + coupon_period
                                                                        * (i + 1), strike_i);
    }

    return option_value;
}

