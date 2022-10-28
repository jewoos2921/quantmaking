//
// Created by jewoo on 2022-10-28.
//
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>

#include "index.h"
#include "yield.h"
#include "product.h"
#include "priceAlpha.h"
#include "mathlib.h"


CPriceAlpha::CPriceAlpha() {
    m_price = 0.0;
    m_alpha = 5.0;
}

CPriceAlpha::~CPriceAlpha() = default;

void CPriceAlpha::black_scholes_option_price_greeks(CIndex &index, CYield &yield, CProduct &product) {

    std::string option_type = product.m_option_type;
    double S = index.m_spot;
    double X = product.m_strike;
    double r = yield.m_riskfree;
    double q = index.m_dividend;
    double sigma = index.m_vol;
    double T = product.m_maturity;

    int iop;

    iop = -1;
    if (option_type == "Call" || option_type == "call") { iop = 1; }

    m_price = black_scholes_option_price(iop, S, X, r, q, sigma, T);
    black_scholes_option_greeks(iop, S, X, r, q, sigma, T);
}

void CPriceAlpha::simulation_european_option_price_greeks(CIndex &index, CYield &yield, CProduct &product, int n_sim) {
    std::string option_type = product.m_option_type;
    double S = index.m_spot;
    double X = product.m_strike;
    double r = yield.m_riskfree;
    double q = index.m_dividend;
    double sigma = index.m_vol;
    double T = product.m_maturity;

    int iop;

    iop = -1;
    if (option_type == "Call" || option_type == "call") { iop = 1; }

    int i;
    double St, price, mudt, ssqrtdt;
    double deltasum, gammasum, vegasum, rhosum;
    double delta, gamma, vega, rho, theta;
    double z, delta_w, gamma_w, vega_w, rho_w, lprice, sigmat;

    price = 0.0;
    mudt = (r - q - 0.5 * sigma * sigma) * T;
    ssqrtdt = sigma * std::sqrt(T);
    sigmat = sigma * T;

    deltasum = gammasum = vegasum = rhosum = 0.0;

    for (i = 0; i < n_sim; ++i) {
        z = normdistrand_BoxMuller();
        St = S * std::exp(mudt + ssqrtdt * z);

        delta_w = z / (S * ssqrtdt);
        gamma_w = (z * z - z * ssqrtdt - 1.0) / (S * ssqrtdt * S * ssqrtdt);
        vega_w = (z * z - 1.0) / sigma - z * std::sqrt(T);
        rho_w = (z * std::sqrt(T) / sigma - T);

        lprice = MAX(iop * (St - X), 0.0);
        price += lprice;
        deltasum += lprice * delta_w;
        gammasum += lprice * gamma_w;
        vegasum += lprice * vega_w;
        rhosum += lprice * rho_w;
    }

    m_price = std::exp(-r * T) * price / n_sim;
    delta = std::exp(-r * T) * deltasum / n_sim;
    gamma = std::exp(-r * T) * gammasum / n_sim;
    vega = std::exp(-r * T) * vegasum / n_sim;
    rho = std::exp(-r * T) * rhosum / n_sim;

    theta = r * m_price - (r - q) * S * delta - 0.5 * sigma * sigma * S * S * gamma;

    vega /= 100.0; // sigma 1% 이동시 가격 영향
    rho /= 100.0; // r 1% 이동시 가격 영향
    theta /= 365.0; // t 1일 이동시 가격 영향

    std::cout << "Delta : " << delta << std::endl;
    std::cout << "Gamma : " << gamma << std::endl;
    std::cout << "Vega : " << vega << std::endl;
    std::cout << "Rho : " << rho << std::endl;
    std::cout << "Theta : " << theta << std::endl;
}

void
CPriceAlpha::binomial_tree_european_option_price_greeks(CIndex &index, CYield &yield, CProduct &product, int n_step) {

    std::string option_type = product.m_option_type;
    double S = index.m_spot;
    double X = product.m_strike;
    double r = yield.m_riskfree;
    double q = index.m_dividend;
    double sigma = index.m_vol;
    double T = product.m_maturity;

    int iop;

    iop = -1;
    if (option_type == "Call" || option_type == "call") { iop = 1; }

    double alpha = m_alpha;

    int i, j, ccnode;
    double dx, nu, dfactor, Pu, Pd, dt;
    double delta, gamma, vega, rho, theta;
    double smax, smin, nodemax, nodemin;
    double v22, vsp, vsm, vrp, vrm;

    dt = T / n_step;
    nu = r - q - 0.5 * sigma * sigma;
    dx = std::sqrt(sigma * sigma * dt + nu * nu * dt * dt);

    Pu = 0.5 + 0.5 * nu * dt / dx;
    Pd = 1.0 - Pu;

    smax = S * std::exp(nu * T + alpha * sigma * std::sqrt(T));
    smin = S * std::exp(nu * T - alpha * sigma * std::sqrt(T));
    ccnode = 0;

    nodemax = nodemin = S;
    do {
        ccnode++;
        nodemax = S * std::exp(ccnode * dx);
        nodemin = S * std::exp(-ccnode * dx);
    } while (nodemax < smax || nodemin > smin);

    double *Sp, *OCV, *NCV;

    Sp = new double[ccnode + 1];
    OCV = new double[ccnode + 1];
    NCV = new double[ccnode + 1];

    n_step += 2; // greeks 계산을 위해 2dt 추가

    if ((n_step - ccnode) % 2 == 0) {
        for (i = 0; i <= ccnode; ++i) {
            Sp[i] = S * std::exp((ccnode - 2 * i) * dx); // Tree 주가 생성
            OCV[i] = MAX(iop * (Sp[i] - X), 0.0); // 만기시점 옵션 Payoff
        }
    } else {
        for (i = 1; i <= ccnode; ++i) {
            Sp[i] = S * std::exp((ccnode - 2 * i + 1) * dx); // Tree 주가 생성
            OCV[i] = MAX(iop * (Sp[i] - X), 0.0); // 만기시점 옵션 Payoff
        }
    }

    dfactor = std::exp(-r * dt);

    for (i = n_step - 1; i >= 2; --i) { // t=2dt까지만 계산하면 됨
        if (i >= ccnode) {
            if ((i - ccnode) % 2 == 0) {
                NCV[0] = OCV[1] * dfactor;
                for (j = 1; j < ccnode; ++j) {
                    NCV[j] = (Pu * OCV[j] + Pd * OCV[j + 1]) * dfactor;
                }
                NCV[ccnode] = OCV[ccnode - 1] * dfactor;
            } else {
                for (j = 1; j <= ccnode; ++j) {
                    NCV[j] = (Pu * OCV[j - 1] + Pd * OCV[j]) * dfactor;
                }
            }
            for (j = 0; j <= ccnode; ++j) {
                OCV[j] = NCV[j];
            }
        } else {
            for (j = 0; j <= i; ++j) {
                OCV[j] = (Pu * OCV[j] + Pd * OCV[j + 1]) * dfactor;
            }
        }
        if (i == 4) { v22 = OCV[2]; }
    }

    double ds, dr;
    ds = dr = 0.0001;

    m_price = OCV[1]; // 옵션의 현재가치
    delta = (OCV[0] - OCV[2]) / (S * std::exp(2 * dt) - S * std::exp(-2 * dx));
    gamma = (((OCV[0] - OCV[1]) / (S * std::exp(2 * dx) - S)) - ((OCV[1] - OCV[2]) / (S - S * std::exp(-2 * dx)))) /
            ((S * std::exp(2 * dx) - S * std::exp(-2 * dx)) / 2.0);
    theta = (v22 - OCV[1]) / (2.0 * dt);
    n_step -= 2; // greeks 계산을 위해 추가한 2dt를 다시 빼줌

    vsp = ci_binomial_tree_european_option_space(iop, S, X, r, q, sigma + ds, T, n_step, alpha);
    vsm = ci_binomial_tree_european_option_space(iop, S, X, r, q, sigma - ds, T, n_step, alpha);
    vrp = ci_binomial_tree_european_option_space(iop, S, X, r + dr, q, sigma, T, n_step, alpha);
    vrm = ci_binomial_tree_european_option_space(iop, S, X, r - dr, q, sigma, T, n_step, alpha);

    vega = (vsp - vsm) / (2.0 * ds);
    rho = (vrp - vrm) / (2.0 * dr);

    delete[]Sp;
    delete[]OCV;
    delete[]NCV;

    vega /= 100.0; // sigma 1% 이동시 가격 영향
    rho /= 100.0; // r 1% 이동시 가격 영향
    theta /= 365.0; // t 1일 이동시 가격 영향

    std::cout << "Delta : " << delta << std::endl;
    std::cout << "Gamma : " << gamma << std::endl;
    std::cout << "Vega : " << vega << std::endl;
    std::cout << "Rho : " << rho << std::endl;
    std::cout << "Theta : " << theta << std::endl;
}

void
CPriceAlpha::trinomial_tree_european_option_price_greeks(CIndex &index, CYield &yield, CProduct &product, int n_step) {

    std::string option_type = product.m_option_type;
    double S = index.m_spot;
    double X = product.m_strike;
    double r = yield.m_riskfree;
    double q = index.m_dividend;
    double sigma = index.m_vol;
    double T = product.m_maturity;

    int iop;

    iop = -1;
    if (option_type == "Call" || option_type == "call") { iop = 1; }

    double alpha = m_alpha;

    int i, j, ccnode, nstep;

    double dx, nu, dfactor, Pu, Pm, Pd, dt;
    double smax, smin, nodemax, nodemin;
    double delta, gamma, vega, rho, theta;
    double v12, vsp, vsm, vrp, vrm;

    dt = T / n_step;
    dx = std::sqrt(3.0 * dt) * sigma;
    nu = r - q - 0.5 * sigma * sigma;

    Pu = 0.5 * ((sigma * sigma * dt + nu * nu * dt * dt) / (dx * dx) + nu * dt / dx);
    Pd = 0.5 * ((sigma * sigma * dt + nu * nu * dt * dt) / (dx * dx) - nu * dt / dx);
    Pm = 1 - Pu - Pd;

    smax = S * std::exp(nu * T + alpha * sigma * std::sqrt(T));
    smin = S * std::exp(nu * T - alpha * sigma * std::sqrt(T));
    ccnode = 0;

    nodemin = nodemax = S;

    do {
        ccnode++;
        nodemax = S * std::exp(ccnode * dx);
        nodemin = S * std::exp(-ccnode * dx);
    } while (nodemax < smax || nodemin > smin);

    double *Sp, *OCV, *NCV;

    nstep = 2 * ccnode + 1;
    Sp = new double[nstep];
    OCV = new double[nstep];
    NCV = new double[nstep];

    n_step += 1; // greeks 계산을 위해 2dt 추가


    for (i = 0; i < nstep; ++i) {
        Sp[i] = S * std::exp((ccnode - i) * dx); // Tree 주가 생성
        OCV[i] = MAX(iop * (Sp[i] - X), 0.0); // 만기시점 옵션 Payoff
    }

    dfactor = std::exp(-r * dt);

    for (i = n_step - 1; i >= 1; --i) { // t=1dt까지만 계산하면 됨

        if (i >= ccnode) {
            NCV[0] = (Pm * OCV[0] + Pd * OCV[1]) / (Pm + Pd) * dfactor;
            for (j = 1; j < nstep - 1; ++j) {
                NCV[j] = (Pu * OCV[j - 1] + Pm * OCV[j] + Pd * OCV[j + 1]) * dfactor;
            }
            NCV[nstep - 1] = (Pu * OCV[nstep - 2] + Pm * OCV[nstep - 1]) / (Pu + Pm) * dfactor;
            for (j = 0; j < nstep; ++j) { OCV[j] = NCV[j]; }
        } else {
            for (j = 0; j < 2 * i + 1; ++j) {
                OCV[j] = (Pu * OCV[j] + Pm * OCV[j + 1] + Pd * OCV[j + 2]) * dfactor;
            }
        }
        if (i == 2) { v12 = OCV[2]; }
    }

    double ds, dr;
    ds = dr = 0.0001;

    m_price = OCV[1]; // 옵션의 현재가치

    delta = (OCV[0] - OCV[2]) / (S * std::exp(dt) - S * std::exp(-dx));
    gamma = (((OCV[0] - OCV[1]) / (S * std::exp(dx) - S))
             - ((OCV[1] - OCV[2]) / (S - S * std::exp(-dx)))) /
            ((S * std::exp(dx) - S * std::exp(-dx)) / 2.0);

    theta = (v12 - OCV[1]) / dt;

    n_step -= 1; // greeks 계산을 위해 추가한 1dt를 다시 빼줌

    vsp = ci_trinomial_tree_european_option_space(iop, S, X, r, q, sigma + ds, T, n_step, alpha);
    vsm = ci_trinomial_tree_european_option_space(iop, S, X, r, q, sigma - ds, T, n_step, alpha);
    vrp = ci_trinomial_tree_european_option_space(iop, S, X, r + dr, q, sigma, T, n_step, alpha);
    vrm = ci_trinomial_tree_european_option_space(iop, S, X, r - dr, q, sigma, T, n_step, alpha);

    vega = (vsp - vsm) / (2.0 * ds);
    rho = (vrp - vrm) / (2.0 * dr);

    delete[]Sp;
    delete[]OCV;
    delete[]NCV;

    vega /= 100.0; // sigma 1% 이동시 가격 영향
    rho /= 100.0; // r 1% 이동시 가격 영향
    theta /= 365.0; // t 1일 이동시 가격 영향

    std::cout << "Delta : " << delta << std::endl;
    std::cout << "Gamma : " << gamma << std::endl;
    std::cout << "Vega : " << vega << std::endl;
    std::cout << "Rho : " << rho << std::endl;
    std::cout << "Theta : " << theta << std::endl;
}

void CPriceAlpha::implicit_fdm_european_option_price_greeks(CIndex &index, CYield &yield, CProduct &product,
                                                            int n_step, int n_spot) {
    std::string option_type = product.m_option_type;
    double S = index.m_spot;
    double X = product.m_strike;
    double r = yield.m_riskfree;
    double q = index.m_dividend;
    double sigma = index.m_vol;
    double T = product.m_maturity;

    int iop = -1;
    if (option_type == "Call" || option_type == "call") { iop = 1; }

    double alpha = m_alpha;

    int i, j, ccnode, nspot, t;
    double dt, dx, edx;
    double nu, Pu, Pd, Pm;
    double smax, smin, nodemax, nodemin;
    double delta, gamma, vega, rho, theta;
    double oldv, vsp, vsm, vrp, vrm;

    dt = T / n_step;
    dx = sigma * std::sqrt(3.0 * dt);
    edx = std::exp(dx);


    nu = r - q - 0.5 * sigma * sigma;
    Pu = -0.5 * dt * (sigma * sigma / (dx * dx) + nu / dx);
    Pm = 1.0 + dt * (sigma * sigma / (dx * dx) + r);
    Pd = -0.5 * dt * (sigma * sigma / (dx * dx) - nu / dx);

    smax = S * std::exp(nu * T + alpha * sigma * std::sqrt(T));
    smin = S * std::exp(nu * T - alpha * sigma * std::sqrt(T));
    ccnode = 0;

    nodemin = nodemax = S;

    do {
        ccnode++;
        nodemax = S * std::exp(ccnode * dx);
        nodemin = S * std::exp(-ccnode * dx);
    } while (nodemax < smax || nodemin > smin);

    nspot = 2 * ccnode + 1;
    if ((nspot % 2) == 0) { nspot = nspot + 1; }
    auto *Sp = new double[nspot];
    auto *CV = new double[nspot];


    for (i = 0; i < nspot; ++i) { Sp[i] = S * std::exp((ccnode - i) * dx); } // Tree 주가 생성

    double **smatrix, **tmp;
    double *known_value, *unknown_value;

    smatrix = new double *[nspot - 2]; // row
    tmp = new double *[nspot - 2];

    known_value = new double[nspot - 2];
    unknown_value = new double[nspot - 2];

    for (i = 0; i < nspot - 2; ++i) {
        smatrix[i] = new double[nspot - 2]; // column
        tmp[i] = new double[nspot - 2];
    }

    for (i = 0; i < nspot - 2; ++i) { for (j = 0; j < nspot - 2; ++j) { smatrix[i][j] = 0.0; }}

    double bu, bd;

    if (iop == 1) { // 콜옵션
        bu = Sp[0] - Sp[1]; // S가 충분히 크면 Delta = 1 이므로
        bd = 0.0; // / S가 충분히 작으면  Delta = 0 이므로
    } else { // 풋옵션
        bu = 0.0; // S가 충분히 크면 Delta = 0 이므로
        bd = -(Sp[nspot - 2] - Sp[nspot - 1]); // S가 충분히 크면 Delta = -1 이므로
    }

    // 0~n_spot - 3은 n_spot - 2 개
    smatrix[0][0] = Pm;
    smatrix[0][1] = Pd;

    for (i = 1; i < nspot - 3; ++i) {
        smatrix[i][i - 1] = Pu;
        smatrix[i][i] = Pm;
        smatrix[i][i + 1] = Pd;
    }

    smatrix[nspot - 3][nspot - 4] = Pu;
    smatrix[nspot - 3][nspot - 3] = Pm;

    // 만기시점에서의 경계조건 부여
    for (i = 0; i < nspot; ++i) { CV[i] = MAX(iop * (Sp[i] - X), 0); }

    // 만기에서부터 Backward Induction 실행하여 현재가치 계산
    for (t = n_step - 1; t >= 0; --t) {
        for (i = 0; i < nspot - 2; ++i) {
            for (j = 0; j < nspot - 2; ++j) { tmp[i][j] = smatrix[i][j]; }
            known_value[i] = CV[i + 1];
        }

        CV[0] = bu; // 경계조건 적용
        CV[n_spot - 1] = bd; // 경계조건 적용
        known_value[0] -= Pu * CV[0];
        known_value[nspot - 3] -= Pd * CV[n_spot - 1];

        // t-1 시점 미지값 계산
        tridiagonal_elimination(tmp, known_value, unknown_value, nspot - 2);
        for (i = 0; i < nspot - 2; ++i) { CV[i + 1] = unknown_value[i]; }

        if (t == 1) { oldv = CV[(nspot - 1) / 2]; }
    }

    double ds, dr;
    ds = dr = 0.0001;
    m_price = CV[(nspot - 1) / 2];


    delta = (CV[(nspot - 1) / 2 - 1] - CV[(nspot - 1) / 2 + 1]) / (S * std::exp(dt) - S * std::exp(-dx));

    gamma = (((CV[(nspot - 1) / 2 - 1] - CV[(nspot - 1) / 2]) / (S * std::exp(dx) - S))
             - ((CV[(nspot - 1) / 2] - CV[(nspot - 1) / 2 + 1])
                / (S - S * std::exp(-dx)))) / ((S * std::exp(dx) - S * std::exp(-dx)) / 2.0);

    theta = (oldv - CV[(nspot - 1) / 2]) / dt;


    vsp = ci_implicit_fdm_european_option_space(iop, S, X, r, q, sigma + ds, T, n_step, n_spot, alpha);
    vsm = ci_implicit_fdm_european_option_space(iop, S, X, r, q, sigma - ds, T, n_step, n_spot, alpha);
    vrp = ci_implicit_fdm_european_option_space(iop, S, X, r + dr, q, sigma, T, n_step, n_spot, alpha);
    vrm = ci_implicit_fdm_european_option_space(iop, S, X, r - dr, q, sigma, T, n_step, n_spot, alpha);

    vega = (vsp - vsm) / (2.0 * ds);
    rho = (vrp - vrm) / (2.0 * dr);

    delete[] Sp;
    delete[] CV;

    for (i = 0; i < nspot - 2; ++i) {
        delete[] smatrix[i];
        delete[] tmp[i];
    }

    delete[] smatrix;
    delete[] tmp;
    delete[] known_value;
    delete[] unknown_value;

    vega /= 100.0; // sigma 1% 이동시 가격 영향
    rho /= 100.0; // r 1% 이동시 가격 영향
    theta /= 365.0; // t 1일 이동시 가격 영향

    std::cout << "Delta : " << delta << std::endl;
    std::cout << "Gamma : " << gamma << std::endl;
    std::cout << "Vega : " << vega << std::endl;
    std::cout << "Rho : " << rho << std::endl;
    std::cout << "Theta : " << theta << std::endl;
}

double
CPriceAlpha::black_scholes_option_price(int iop, double S, double X, double r, double q, double sigma, double T) {
    double d1, d2, price;

    d1 = (std::log(S / X) + (r - q + (sigma * sigma) / 2) * T) / (sigma * std::sqrt(T));
    d2 = d1 - sigma * std::sqrt(T);

    price = iop * (S * std::exp(-q * T) * N(iop * d1) - X * std::exp(-r * T) * N(iop * d2)); // Black-scholes

    return price;
}

void CPriceAlpha::black_scholes_option_greeks(int iop, double S, double X, double r, double q, double sigma, double T) {
    double d1, d2;
    double delta, gamma, vega, rho, theta;

    d1 = (std::log(S / X) + (r - q + (sigma * sigma) / 2.0) * T) / (sigma * std::sqrt(T));
    d2 = d1 - sigma * std::sqrt(T);

    delta = iop * std::exp(-q * T) * N(iop * d1);
    gamma = n(d1) * std::exp(-q * T) / (S * sigma * std::sqrt(T));
    vega = S * n(d1) * std::exp(-q * T) * std::sqrt(T);
    rho = iop * X * T * std::exp(-r * T) * N(iop * d2);
    theta = -S * n(d1) * std::exp(-q * T) * sigma / (2 * std::sqrt(T))
            + iop * (q * S * std::exp(-q * T)
                     * N(iop * d1) - r * X *
                                     std::exp(-r * T) * N(iop * d2));

    vega /= 100.0; // sigma 1% 이동시 가격 영향
    rho /= 100.0; // r 1% 이동시 가격 영향
    theta /= 365.0; // t 1일 이동시 가격 영향

    std::cout << "Delta : " << delta << std::endl;
    std::cout << "Gamma : " << gamma << std::endl;
    std::cout << "Vega : " << vega << std::endl;
    std::cout << "Rho : " << rho << std::endl;
    std::cout << "Theta : " << theta << std::endl;
}

double
CPriceAlpha::ci_binomial_tree_european_option_space(int iop, double S, double X, double r, double q, double sigma,
                                                    double T, int n_step, double alpha) {
    int i, j, ccnode;
    double dx, nu, dfactor, Pu, Pd, dt;
    double smax, smin, nodemax, nodemin;
    double price;


    dt = T / n_step;
    nu = r - q - 0.5 * sigma * sigma;
    dx = std::sqrt(sigma * sigma * dt + nu * nu * dt * dt);

    Pu = 0.5 + 0.5 * nu * dt / dx;
    Pd = 1.0 - Pu;

    smax = S * std::exp(nu * T + alpha * sigma * std::sqrt(T));
    smin = S * std::exp(nu * T - alpha * sigma * std::sqrt(T));
    ccnode = 0;


    nodemax = nodemin = S;

    do {
        ccnode++;
        nodemax = S * std::exp(ccnode * dx);
        nodemin = S * std::exp(-ccnode * dx);

    } while (nodemax < smax || nodemin > smin);

    double *Sp;
    double *OCV, *NCV;

    Sp = new double[ccnode + 1];
    OCV = new double[ccnode + 1];
    NCV = new double[ccnode + 1];


    if ((n_step - ccnode) % 2 == 0) {
        for (i = 0; i <= ccnode; ++i) {
            Sp[i] = S * std::exp((ccnode - 2 * i) * dx); // Tree 주가 생성
            OCV[i] = MAX(iop * (Sp[i] - X), 0.0); // 만기 시점 옵션 payoff
        }
    } else {
        for (i = 1; i <= ccnode; ++i) {
            Sp[i] = S * std::exp((ccnode - 2 * i + 1) * dx); // Tree 주가 생성
            OCV[i] = MAX(iop * (Sp[i] - X), 0.0); // 만기 시점 옵션 payoff
        }
    }

    dfactor = std::exp(-r * dt);
    // Tree 의 만기에서 부터 Backward Induction 실행하여 현재가치 계산
    for (i = n_step - 1; i >= 0; --i) {
        if (i >= ccnode) {
            if ((i - ccnode) % 2 == 0) {
                NCV[0] = OCV[1] * dfactor;

                for (j = 1; j < ccnode; ++j) {
                    NCV[j] = (Pu * OCV[j] + Pd * OCV[j + 1]) * dfactor;
                }

                NCV[ccnode] = OCV[ccnode - 1] * dfactor;
            } else {
                for (j = 1; j <= ccnode; ++j) {
                    NCV[j] = (Pu * OCV[j - 1] + Pd * OCV[j]) * dfactor;
                }
            }
            for (j = 0; j <= ccnode; ++j) { OCV[j] = NCV[j]; }
        } else {
            for (j = 0; j <= i; ++j) {
                OCV[j] = (Pu * OCV[j] + Pd * OCV[j + 1]) * dfactor;
            }
        }
    }

    price = OCV[0]; // 옵션의 현재가치

    delete[]Sp;
    delete[]OCV;
    delete[]NCV;

    return price;
}

double
CPriceAlpha::ci_trinomial_tree_european_option_space(int iop, double S, double X, double r, double q, double sigma,
                                                     double T, int n_step, double alpha) {
    int i, j, ccnode, nstep;
    double dx, nu, dfactor, Pu, Pd, dt, Pm;
    double smax, smin, nodemax, nodemin;
    double price;

    dt = T / n_step;
    dx = sigma * std::sqrt(3.0 * dt);
    nu = r - q - 0.5 * sigma * sigma;

    Pu = 0.5 * ((sigma * sigma * dt + nu * nu * dt * dt) / (dx * dx) + nu * dt / dx);
    Pd = 0.5 * ((sigma * sigma * dt + nu * nu * dt * dt) / (dx * dx) - nu * dt / dx);
    Pm = 1 - Pu - Pd;

    smax = S * std::exp(nu * T + alpha * sigma * std::sqrt(T));
    smin = S * std::exp(nu * T - alpha * sigma * std::sqrt(T));
    ccnode = 0;

    nodemin = nodemax = S;

    do {
        ccnode++;
        nodemax = S * std::exp(ccnode * dx);
        nodemin = S * std::exp(-ccnode * dx);
    } while (nodemax < smax || nodemin > smin);

    double *Sp;
    double *OCV, *NCV;

    nstep = 2 * ccnode + 1;

    Sp = new double[nstep];
    OCV = new double[nstep];
    NCV = new double[nstep];


    for (i = 0; i < nstep; ++i) {
        Sp[i] = S * std::exp((ccnode - i) * dx); // Tree 주가 생성
        OCV[i] = MAX(iop * (Sp[i] - X), 0.0); // 만기 페이오프
    }

    dfactor = std::exp(-r * dt);

    // Tree 의 만기에서 부터 Backward Induction 실행하여 현재가치 계산
    for (i = n_step - 1; i >= 0; --i) {
        if (i >= ccnode) {
            NCV[0] = (Pm * OCV[0] + Pd * OCV[1]) / (Pm + Pd) * dfactor;
            for (j = 1; j < nstep - 1; ++j) {
                NCV[j] = (Pu * OCV[j - 1] + Pm * OCV[j] + Pd * OCV[j + 1]) * dfactor;
            }

            NCV[nstep - 1] = (Pu * OCV[nstep - 2] + Pm * OCV[nstep - 1]) / (Pu + Pm) * dfactor;

            for (j = 0; j < nstep; ++j) { OCV[j] = NCV[j]; }
        } else {
            for (j = 0; j < 2 * i + 1; ++j) {
                OCV[j] = (Pu * OCV[j] + Pm * OCV[j + 1] + Pd * OCV[j + 2]) * dfactor;
            }
        }
    }

    price = OCV[0]; // 옵션의 현재가치

    delete[]Sp;
    delete[]OCV;
    delete[]NCV;

    return price;
}

double CPriceAlpha::ci_implicit_fdm_european_option_space(int iop, double S, double X, double r, double q, double sigma,
                                                          double T, int n_step, int n_spot, double alpha) {
    int i, j, ccnode, nspot, t;
    double dx, nu, Pu, Pd, dt, Pm, edx;
    double smax, smin, nodemax, nodemin;

    dt = T / n_step;
    dx = sigma * std::sqrt(3.0 * dt);
    edx = std::exp(dx);

    nu = r - q - 0.5 * sigma * sigma;
    Pu = -0.5 * dt * (sigma * sigma / (dx * dx) + nu / dx);
    Pm = 1.0 + dt * (sigma * sigma / (dx * dx) + r);
    Pd = -0.5 * dt * (sigma * sigma / (dx * dx) - nu / dx);

    smax = S * std::exp(nu * T + alpha * sigma * std::sqrt(T));
    smin = S * std::exp(nu * T - alpha * sigma * std::sqrt(T));
    ccnode = 0;

    nodemin = nodemax = S;

    do {
        ccnode++;
        nodemax = S * std::exp(ccnode * dx);
        nodemin = S * std::exp(-ccnode * dx);
    } while (nodemax < smax || nodemin > smin);

    nspot = 2 * ccnode + 1;

    if ((nspot % 2) == 0) { nspot = nspot + 1; }

    auto *Sp = new double[nspot];
    auto *CV = new double[nspot];

    for (i = 0; i < nspot; ++i) { Sp[i] = S * std::exp((ccnode - i) * dx); } // Tree 주가 생성


    double **smatrix, **tmp;
    double *known_value, *unknown_value;

    smatrix = new double *[nspot - 2]; // row
    tmp = new double *[nspot - 2];

    known_value = new double[nspot - 2];
    unknown_value = new double[nspot - 2];

    for (i = 0; i < nspot - 2; ++i) {
        smatrix[i] = new double[nspot - 2]; // column
        tmp[i] = new double[nspot - 2];
    }

    for (i = 0; i < nspot - 2; ++i) { for (j = 0; j < nspot - 2; ++j) { smatrix[i][j] = 0.0; }}

    double bu, bd;

    if (iop == 1) { // 콜옵션
        bu = Sp[0] - Sp[1]; // S가 충분히 크면 Delta = 1 이므로
        bd = 0.0; // / S가 충분히 작으면  Delta = 0 이므로
    } else { // 풋옵션
        bu = 0.0; // S가 충분히 크면 Delta = 0 이므로
        bd = -(Sp[nspot - 2] - Sp[nspot - 1]); // S가 충분히 크면 Delta = -1 이므로
    }

    // 0~n_spot - 3은 n_spot - 2 개
    smatrix[0][0] = Pm;
    smatrix[0][1] = Pd;

    for (i = 1; i < nspot - 3; ++i) {
        smatrix[i][i - 1] = Pu;
        smatrix[i][i] = Pm;
        smatrix[i][i + 1] = Pd;
    }

    smatrix[nspot - 3][nspot - 4] = Pu;
    smatrix[nspot - 3][nspot - 3] = Pm;

    // 만기시점에서의 경계조건 부여
    for (i = 0; i < nspot; ++i) { CV[i] = MAX(iop * (Sp[i] - X), 0); }

    // 만기에서부터 Backward Induction 실행하여 현재가치 계산
    for (t = n_step - 1; t >= 0; --t) {
        for (i = 0; i < nspot - 2; ++i) {
            for (j = 0; j < nspot - 2; ++j) { tmp[i][j] = smatrix[i][j]; }
            known_value[i] = CV[i + 1];
        }

        CV[0] = bu; // 경계조건 적용
        CV[n_spot - 1] = bd; // 경계조건 적용
        known_value[0] -= Pu * CV[0];
        known_value[nspot - 3] -= Pd * CV[n_spot - 1];

        // t-1 시점 미지값 계산
        tridiagonal_elimination(tmp, known_value, unknown_value, nspot - 2);
        for (i = 0; i < nspot - 2; ++i) { CV[i + 1] = unknown_value[i]; }
    }

    double price;

    price = CV[(nspot - 1) / 2];
    delete[]Sp;
    delete[]CV;

    for (i = 0; i < nspot - 2; ++i) {
        delete[] smatrix[i];
        delete[] tmp[i];
    }

    delete[] smatrix;
    delete[] tmp;
    delete[] known_value;
    delete[] unknown_value;

    return price;
}




