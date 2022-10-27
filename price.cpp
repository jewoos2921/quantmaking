//
// Created by jewoo on 2022-10-23.
//
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include "index.h"
#include "yield.h"
#include "product.h"
#include "price.h"
#include "mathlib.h"
#include "sobol.h"

CPrice::CPrice() {
    m_price = 0;
}

CPrice::~CPrice() {
}

void CPrice::black_scholes_option_price(CIndex &index, CYield &yield, CProduct &product) {
    std::string option_type = product.m_option_type;
    double S = index.m_spot;
    double X = product.m_strike;
    double r = yield.m_riskfree;
    double q = index.m_dividend;
    double sigma = index.m_vol;
    double T = product.m_maturity;

    int iop;
    double d1, d2;

    iop = -1;
    if (option_type == "Call" || option_type == "call") {
        iop = 1;
    }
    d1 = (std::log(S / X) + (r - q + (sigma * sigma) / 2) * T) / (sigma * std::sqrt(T));
    d2 = d1 - sigma * std::sqrt(T);

    // 블랙 숄츠 모델
    m_price = iop * (S * std::exp(-q * T) * N(iop * d1) - X * std::exp(-r * T) * N(iop * d2));
}

void CPrice::simulation_european_option_price(CIndex &index, CYield &yield, CProduct &product, int n_sim) {
    std::string option_type = product.m_option_type;
    double S = index.m_spot;
    double X = product.m_strike;
    double r = yield.m_riskfree;
    double q = index.m_dividend;
    double sigma = index.m_vol;
    double T = product.m_maturity;

    int i, iop;
    double St, price, mudt, ssqrtdt;
    long idum = clock();

    price = 0.0;
    mudt = (r - q * 0.5 * sigma * sigma) * T;
    ssqrtdt = sigma * std::sqrt(T);

    iop = -1;
    if (option_type == "Call" || option_type == "call") {
        iop = 1;
    }

    std::srand(static_cast<unsigned >(std::time(nullptr))); // Random seed 초기화

    for (i = 0; i < n_sim; ++i) {
        St = S * std::exp(mudt + ssqrtdt * normdistrand_BoxMuller()); // 주가 생성
        price += MAX(iop * (St - X), 0.0); // 만기 payoff
    }

    m_price = (price / n_sim) * std::exp(-r * T); // 만기 payoff의 평균을 현가

}

void CPrice::binomial_tree_european_option_price(CIndex &index, CYield &yield, CProduct &product, int n_step) {
    std::string option_type = product.m_option_type;
    double S = index.m_spot;
    double X = product.m_strike;
    double r = yield.m_riskfree;
    double q = index.m_dividend;
    double sigma = index.m_vol;
    double T = product.m_maturity;

    int i, j, iop;
    double up, down, dfactor, Pu, Pd, dt;

    dt = T / n_step;
    up = std::exp(sigma * std::sqrt(dt));
    down = 1.0 / up;

    Pu = (std::exp((r - q) * dt) - down) / (up - down); // up 확률
    Pd = 1.0 - Pu; // down 확률

    double *Sp, *Cv;

    Sp = new double[n_step + 1];
    Cv = new double[n_step + 1];

    // 트리 주가 생성
    Sp[0] = S * std::pow(up, n_step);
    for (i = 1; i <= n_step; ++i) {
        Sp[i] = Sp[i - 1] * down * down;
    }

    iop = -1;
    if (option_type == "Call" || option_type == "call") {
        iop = 1;
    }

    // 만기시점 옵션 Payoff
    for (i = 0; i <= n_step; ++i) {
        Cv[i] = MAX(iop * (Sp[i] - X), 0.0);
    }

    dfactor = std::exp(-r * dt);
    // Tree의 만기에서부터 Backward Induction 실행하여 현재가치 계싼
    for (i = n_step - 1; i >= 0; --i) {
        for (j = 0; j <= i; ++j) {
            Cv[j] = (Pu * Cv[j] + Pd * Cv[j + 1]) * dfactor;
        }
    }
    m_price = Cv[0]; // 옵션의 현재 가치
    delete[]Sp;
    delete[]Cv;

}

void
CPrice::lognormal_binomial_tree_european_option_price(CIndex &index, CYield &yield, CProduct &product, int n_step) {
    std::string option_type = product.m_option_type;
    double S = index.m_spot;
    double X = product.m_strike;
    double r = yield.m_riskfree;
    double q = index.m_dividend;
    double sigma = index.m_vol;
    double T = product.m_maturity;

    int i, j, iop;
    double dx, nu, dfactor, Pu, Pd, dt;


    dt = T / n_step;
    nu = r - q - 0.5 * sigma * sigma;
    dx = std::sqrt(sigma * sigma * dt + nu * nu * dt * dt);

    Pu = 0.5 + 0.5 * nu * dt / dx; // up 확률
    Pd = 1.0 - Pu; // down 확률

    double *Sp, *Cv;

    Sp = new double[n_step + 1];
    Cv = new double[n_step + 1];

    for (i = 0; i <= n_step; ++i) {
        Sp[i] = S * std::exp((n_step - 2 * i) * dx); // Tree 주가 생성
    }

    iop = -1;
    if (option_type == "Call" || option_type == "call") {
        iop = 1;
    }

    for (i = 0; i < n_step; ++i) {
        Cv[i] = MAX(iop * (Sp[i] - X), 0.0); // 만기 페이오프
    }
    dfactor = std::exp(-r * dt);

    // Tree의 만기에서부터 Backward Induction 실행하여 현재가치 계산
    for (i = n_step - 1; i >= 0; --i) {
        for (j = 0; j <= i; ++j) {
            Cv[j] = (Pu * Cv[j] + Pd * Cv[j + 1]) * dfactor;
        }
    }

    m_price = Cv[0]; // 옵션의 현재 가치

    delete[]Sp;
    delete[]Cv;
}

void CPrice::trinomial_tree_european_option_price(CIndex &index, CYield &yield, CProduct &product, int n_step) {
    std::string option_type = product.m_option_type;
    double S = index.m_spot;
    double X = product.m_strike;
    double r = yield.m_riskfree;
    double q = index.m_dividend;
    double sigma = index.m_vol;
    double T = product.m_maturity;

    int i, j, iop, nstep;
    double dx, nu, dfactor, Pu, Pm, Pd, dt;
    dt = T / n_step;
    dx = sigma * std::sqrt(3.0 * dt);
    nu = r - q - 0.5 * sigma * sigma;

    Pu = 0.5 * ((sigma * sigma * dt + nu * nu * dt * dt) / (dx * dx) + nu * dt / dx); // up 확률
    Pd = 0.5 * ((sigma * sigma * dt + nu * nu * dt * dt) / (dx * dx) - nu * dt / dx); // down 확률
    Pm = 1 - Pu - Pd; // middle 확률

    double *Sp, *Cv;

    nstep = 2 * n_step + 1;

    Sp = new double[nstep];
    Cv = new double[nstep];


    for (i = 0; i < nstep; ++i) {
        Sp[i] = S * std::exp((n_step - i) * dx); // Tree 주가 생성
    }

    iop = -1;
    if (option_type == "Call" || option_type == "call") {
        iop = 1;
    }

    for (i = 0; i < nstep; ++i) {
        Cv[i] = MAX(iop * (Sp[i] - X), 0.0); // 만기 페이오프
    }

    dfactor = std::exp(-r * dt);

    // Tree의 만기에서부터 Backward Induction 실행하여 현재가치 계산
    for (i = n_step - 1; i >= 0; --i) {
        for (j = 0; j <= 2 * i + 1; ++j) {
            Cv[j] = (Pu * Cv[j] + Pm * Cv[j + 1] + Pd * Cv[j + 2]) * dfactor;
        }
    }

    m_price = Cv[0]; // 옵션의 현재 가치

    delete[]Sp;
    delete[]Cv;
}

void
CPrice::implicit_fdm_european_option_price(CIndex &index, CYield &yield, CProduct &product, int n_step, int n_spot) {
    std::string option_type = product.m_option_type;
    double S = index.m_spot;
    double X = product.m_strike;
    double r = yield.m_riskfree;
    double q = index.m_dividend;
    double sigma = index.m_vol;
    double T = product.m_maturity;
    int i, j, t, iop;

    if ((n_spot % 2) == 0) {
        n_spot = n_spot + 1;
    }

    double *Sp, *Cv;

    Sp = new double[n_spot];
    Cv = new double[n_spot];

    double dt, dx, edx;
    dt = T / n_step;
    dx = sigma * std::sqrt(3.0 * dt);
    edx = std::exp(dx);

    Sp[0] = S * std::exp((n_spot - 1) / 2 * dx);
    for (i = 1; i < n_spot; ++i) {
        Sp[i] = Sp[i - 1] / edx;
    }
    Sp[(n_spot - 1) / 2] = S;

    iop = -1;
    if (option_type == "Call" || option_type == "call") {
        iop = 1;
    }

    double Pu, Pm, Pd, nu;

    nu = r - q - 0.5 * sigma * sigma;
    Pu = -0.5 * dt * (sigma * sigma / (dx * dx) + nu / dx);
    Pm = 1.0 + dt * (sigma * sigma / (dx * dx) + r);
    Pd = -0.5 * dt * (sigma * sigma / (dx * dx) - nu / dx);

    double **smatrix, **tmp, *known_value, *unknown_value;
    smatrix = new double *[n_spot - 2];  // row create
    tmp = new double *[n_spot - 2]; // row create
    known_value = new double[n_spot - 2];
    unknown_value = new double[n_spot - 2];
    for (i = 0; i < n_spot - 2; ++i) {
        smatrix[i] = new double[n_spot - 2]; // column create
        tmp[i] = new double[n_spot - 2]; // column create
    }

    for (i = 0; i < n_spot - 2; ++i) {
        for (j = 0; j < n_spot - 2; ++j) {
            smatrix[i][j] = 0.0;
        }
    }

    double bd, bu;
    if (iop == 1) { // 콜옵션
        bu = Sp[0] - Sp[1]; // S가 충분히 크면 Delta = 1 이므로
        bd = 0.0; // / S가 충분히 작으면  Delta = 0 이므로
    } else { // 풋옵션
        bu = 0.0; // S가 충분히 크면 Delta = 0 이므로
        bd = -(Sp[n_spot - 2] - Sp[n_spot - 1]); // S가 충분히 크면 Delta = -1 이므로
    }

    // 0~n_spot - 3은 n_spot - 2 개
    smatrix[0][0] = Pm;
    smatrix[0][1] = Pd;

    for (i = 1; i < n_spot - 3; ++i) {
        smatrix[i][i - 1] = Pu;
        smatrix[i][i] = Pm;
        smatrix[i][i + 1] = Pd;
    }

    smatrix[n_spot - 3][n_spot - 4] = Pu;
    smatrix[n_spot - 3][n_spot - 3] = Pm;

    // 만기시점에서의 경계조건 부여
    for (i = 0; i < n_spot; ++i) {
        Cv[i] = MAX(iop * (Sp[i] - X), 0);
    }

    // 만기에서부터 Backward Induction 실행하여 현재가치 계산
    for (t = n_step - 1; t > 0; --t) {
        for (i = 0; i < n_spot - 2; ++i) {
            for (j = 0; j < n_spot - 2; ++j) {
                tmp[i][j] = smatrix[i][j];
            }

            known_value[i] = Cv[i + 1];
        }
        Cv[0] = bu; // 경계조건 적용
        Cv[n_spot - 1] = bd; // 경계조건 적용
        known_value[0] -= Pu * Cv[0];
        known_value[n_spot - 3] -= Pd * Cv[n_spot - 1];

        // t-1 시점 미지값 계산
        tridiagonal_elimination(tmp, known_value, unknown_value, n_spot - 2);
        for (i = 0; i < n_spot - 2; ++i) {
            Cv[i + 1] = unknown_value[i];
        }
    }

    m_price = Cv[(n_spot - 1) / 2];
    delete[]Sp;
    delete[]Cv;

    for (i = 0; i < n_spot - 2; ++i) {
        delete[] smatrix[i];
        delete[] tmp[i];
    }

    delete[] smatrix;
    delete[] tmp;
    delete[] known_value;
    delete[] unknown_value;
}

/*
 * 스프레드 옵션
*/
void CPrice::spread_option_price(CIndex &index1, CIndex &index2, CYield &yield, CProduct &product) {
    std::string option_type = product.m_option_type;
    double rho = product.m_correlation;
    double X = product.m_strike;
    double T = product.m_maturity;

    double S1 = index1.m_spot;
    double q1 = index1.m_dividend;
    double sigma1 = index1.m_vol;

    double S2 = index2.m_spot;
    double q2 = index2.m_dividend;
    double sigma2 = index2.m_vol;
    double r = yield.m_riskfree;

    int iop;
    double F, S, sigma, d1, d2;

    F = S2 * std::exp(-q2 * T) / (S2 * std::exp(-q2 * T)
                                  + X * std::exp(-r * T));

    sigma = std::sqrt(sigma1 * sigma1 + std::pow((sigma2 * F), 2)
                      - 2 * rho * sigma1 * sigma2 * F);

    S = S1 * std::exp(-q1 * T) / (S2 * std::exp(-q2 * T)
                                  + X * std::exp(-r * T));

    d1 = (std::log(S) + 0.5 * sigma * sigma * T) / (sigma * std::sqrt(T));
    d2 = d1 - sigma * std::sqrt(T);

    iop = -1;
    if (option_type == "Call" || option_type == "call") {
        iop = 1;
    }

    m_price = (S2 * std::exp(-q2 * T)
               + X * std::exp(-r * T)
                 * iop *
                 (S * N(iop * d1) - N(iop * d2)));
}

void
CPrice::simulation_spread_option_price(CIndex &index1, CIndex &index2, CYield &yield, CProduct &product, int n_sim) {

    std::string option_type = product.m_option_type;
    double rho = product.m_correlation;
    double X = product.m_strike;
    double T = product.m_maturity;
    double S1 = index1.m_spot;
    double q1 = index1.m_dividend;
    double sigma1 = index1.m_vol;

    double S2 = index2.m_spot;
    double q2 = index2.m_dividend;
    double sigma2 = index2.m_vol;
    double r = yield.m_riskfree;


    int i, iop;
    double rv1, rv2, epsilon1, epsilon2;
    double St1, St2, price, mu1t, mu2t;
    double s1sqrtt, s2sqrtt, sqrtrho;


    std::srand(static_cast<unsigned >(std::time(nullptr))); // Random seed 초기화

    price = 0.0;

    mu1t = (r - q1 - 0.5 * sigma1 * sigma1) * T;
    mu2t = (r - q2 - 0.5 * sigma2 * sigma2) * T;

    s1sqrtt = sigma1 * std::sqrt(T);
    s2sqrtt = sigma2 * std::sqrt(T);
    sqrtrho = std::sqrt(1.0 - rho * rho);

    iop = -1;
    if (option_type == "Call" || option_type == "call") {
        iop = 1;
    }

    for (i = 0; i < n_sim; ++i) {
        rv1 = normdistrand_BoxMuller();
        rv2 = normdistrand_BoxMuller();

        epsilon1 = rv1;
        epsilon2 = rho * rv1 + sqrtrho * rv2;

        St1 = S1 * std::exp(mu1t + s1sqrtt * epsilon1); // 주가 생성
        St2 = S2 * std::exp(mu2t + s2sqrtt * epsilon2); // 주가 생성

        price += MAX(iop * ((St1 - St2) - X), 0.0); // 만기 payoff
    }

    m_price = (price / n_sim) * std::exp(-r * T); // 만기 payoff의 평균을 현가
}

void CPrice::binomial_tree_spread_option_price(CIndex &index1, CIndex &index2,
                                               CYield &yield, CProduct &product,
                                               int n_step) {
    std::string option_type = product.m_option_type;
    double rho = product.m_correlation;
    double X = product.m_strike;
    double T = product.m_maturity;
    double S1 = index1.m_spot;
    double q1 = index1.m_dividend;
    double sigma1 = index1.m_vol;

    double S2 = index2.m_spot;
    double q2 = index2.m_dividend;
    double sigma2 = index2.m_vol;
    double r = yield.m_riskfree;

    int i, k, j, iop;

    double dx1, dx2, nu1, nu2;
    double p1, p2, p3, p4;
    double dfactor, dt;

    dt = T / n_step;
    dx1 = sigma1 * std::sqrt(dt);
    dx2 = sigma2 * std::sqrt(dt);

    nu1 = r - q1 - 0.5 * sigma1 * sigma1;
    nu2 = r - q2 - 0.5 * sigma2 * sigma2;

    double *Sp1, *Sp2, **CV;

    Sp1 = new double[n_step + 1];
    Sp2 = new double[n_step + 1];
    CV = new double *[n_step + 1];

    for (i = 0; i < n_step + 1; ++i) {
        CV[i] = new double[n_step + 1];
    }

    // Tree 주가 생성
    Sp1[0] = S1 * std::exp(n_step * dx1);
    Sp2[0] = S2 * std::exp(n_step * dx2);

    for (i = 1; i <= n_step; ++i) {
        Sp1[i] = Sp1[i - 1] * std::exp(-2 * dx1);
        Sp2[i] = Sp2[i - 1] * std::exp(-2 * dx2);
    }

    p1 = (dx1 * dx2 +
          (dx2 * nu1 + dx1 * nu2 + rho * sigma1 * sigma2) * dt) / (4 * dx1 * dx2); // s1 up, d2 up
    p2 = (dx1 * dx2 +
          (dx2 * nu1 - dx1 * nu2 - rho * sigma1 * sigma2) * dt) / (4 * dx1 * dx2); // s1 up, d2 down
    p3 = (dx1 * dx2 -
          (dx2 * nu1 + dx1 * nu2 - rho * sigma1 * sigma2) * dt) / (4 * dx1 * dx2); // s1 down, d2 down
    p4 = (dx1 * dx2 -
          (dx2 * nu1 - dx1 * nu2 + rho * sigma1 * sigma2) * dt) / (4 * dx1 * dx2); // s1 down, d2 up

    iop = -1;
    if (option_type == "Call" || option_type == "call") {
        iop = 1;
    }

    // 만기시점 옵션 payoff
    for (i = 0; i <= n_step; ++i) {
        for (j = 0; j <= n_step; ++j) {
            CV[i][j] = MAX(iop * ((Sp1[i] - Sp2[j]) - X), 0.0);
        }
    }

    dfactor = std::exp(-r * dt);

    // Tree의 만기에서부터 Backward Induction 실행하여 현재가치 계산
    for (i = n_step - 1; i >= 0; --i) {
        for (j = 0; j <= i; ++j) {
            for (k = 0; k <= i; ++k) {
                CV[j][k] =
                        (p1 * CV[j][k] + p2 * CV[j][k + 1] +
                         p3 * CV[j + 1][k + 1] + p4 * CV[j + 1][k]) * dfactor;
            }
        }
    }

    m_price = CV[0][0]; // 옵션의 현재가치

    for (i = 0; i < n_step + 1; ++i) {
        delete[] CV[i];
    }

    delete[] CV;
    delete[] Sp1;
    delete[] Sp2;
}

void CPrice::trinomial_tree_spread_option_price(CIndex &index1, CIndex &index2, CYield &yield, CProduct &product,
                                                int n_step) {
    std::string option_type = product.m_option_type;
    double rho = product.m_correlation;
    double X = product.m_strike;
    double T = product.m_maturity;
    double S1 = index1.m_spot;
    double q1 = index1.m_dividend;
    double sigma1 = index1.m_vol;

    double S2 = index2.m_spot;
    double q2 = index2.m_dividend;
    double sigma2 = index2.m_vol;
    double r = yield.m_riskfree;

    int i, k, j, iop, nstep;

    double dx1, dx2, nu1, nu2;
    double p1, p2, p3, p4, p5;
    double dfactor, dt;

    dt = T / n_step;
    dx1 = sigma1 * std::sqrt(3.0 * dt);
    dx2 = sigma2 * std::sqrt(3.0 * dt);

    nu1 = r - q1 - 0.5 * sigma1 * sigma1;
    nu2 = r - q2 - 0.5 * sigma2 * sigma2;

    nstep = 2 * n_step + 1;

    double *Sp1, *Sp2, **CV;
    Sp1 = new double[nstep];
    Sp2 = new double[nstep];
    CV = new double *[nstep];

    for (i = 0; i < nstep; ++i) {
        CV[i] = new double[nstep];
    }

    // Tree 주가 생성
    Sp1[0] = S1 * std::exp(n_step * dx1);
    Sp2[0] = S2 * std::exp(n_step * dx2);

    for (i = 1; i <= n_step; ++i) {
        Sp1[i] = Sp1[i - 1] * std::exp(-dx1);
        Sp2[i] = Sp2[i - 1] * std::exp(-dx2);
    }

    double f1, f2, f3, f4, f5, f6;
    f1 = nu1 * dt / dx1;
    f2 = nu2 * dt / dx2;
    f3 = sigma1 * sigma1 * dt / (dx1 * dx1);
    f4 = nu1 * nu1 * dt * dt / (dx1 * dx1);
    f5 = rho * sigma1 * sigma2 * dt / (dx1 * dx2);
    f6 = nu1 * nu2 * dt * dt / (dx1 * dx2);


    p1 = 0.25 * (f1 + f2 + f3 + f4 + f5 + f5); // s1 up, s2 up
    p2 = 0.25 * (f1 - f2 + f3 + f4 - f5 - f5); // s1 up, s2 down
    p3 = 0.25 * (-f1 - f2 + f3 + f4 + f5 + f5); // s1 down, s2 down
    p4 = 0.25 * (-f1 + f2 + f3 + f4 - f5 - f5); // s1 down, s2 up
    p5 = 1.0 - f3 - f4;                         // s1 middle, s2 middle

    iop = -1;
    if (option_type == "Call" || option_type == "call") {
        iop = 1;
    }

    // 만기시점 옵션 payoff
    for (i = 0; i < nstep; ++i) {
        for (j = 0; j < nstep; ++j) {
            CV[i][j] = MAX(iop * ((Sp1[i] - Sp2[j]) - X), 0.0);
        }
    }

    dfactor = std::exp(-r * dt);

    // Tree의 만기에서부터 Backward Induction 실행하여 현재가치 계산
    for (i = n_step - 1; i >= 0; --i) {
        for (j = 0; j < 2 * i + 1; ++j) {
            for (k = 0; k < 2 * i + 1; ++k) {
                CV[j][k] = (p1 * CV[j][k] + p2 * CV[j][k + 2] +
                            p3 * CV[j + 2][k + 2] + p4 * CV[j + 2][k]
                            + p5 * CV[j + 1][k + 1]) * dfactor;
            }
        }
    }

    m_price = CV[0][0];

    for (i = 0; i < nstep; ++i) {
        delete[] CV[i];
    }

    delete[] CV;
    delete[] Sp1;
    delete[] Sp2;

}

void
CPrice::implicit_fdm_spread_option_price(CIndex &index1, CIndex &index2, CYield &yield, CProduct &product,
                                         int n_step, int n_spot) {
    std::string option_type = product.m_option_type;
    double rho = product.m_correlation;
    double X = product.m_strike;
    double T = product.m_maturity;
    double S1 = index1.m_spot;
    double q1 = index1.m_dividend;
    double sigma1 = index1.m_vol;

    double S2 = index2.m_spot;
    double q2 = index2.m_dividend;
    double sigma2 = index2.m_vol;
    double r = yield.m_riskfree;

    int i, k, j, iop, t;
    double *Sp1, *Sp2, **oldCV, **newCV;

    if ((n_spot % 2) == 0) {
        n_spot = n_spot + 1;
    }

    Sp1 = new double[n_spot];
    Sp2 = new double[n_spot];
    oldCV = new double *[n_spot];
    newCV = new double *[n_spot];

    for (i = 0; i < n_spot; ++i) {
        oldCV[i] = new double[n_spot];
        newCV[i] = new double[n_spot];
    }


    double dt, cl, mu1, mu2, dx1, dx2, edx1, edx2;

    cl = 4.0;
    mu1 = r - q1 - 0.5 * sigma1 * sigma1;
    mu2 = r - q2 - 0.5 * sigma2 * sigma2;
    dt = T / n_step;

    dx1 = sigma1 * std::sqrt(3.0 * dt);
    dx2 = sigma2 * std::sqrt(3.0 * dt);

    edx1 = std::exp(dx1);
    edx2 = std::exp(dx2);

    Sp1[0] = S1 * std::exp((n_spot - 1) / 2 * dx1); // 기초 자산의 노드별 가격 설정
    for (i = 1; i < n_spot; ++i) {
        Sp1[i] = Sp1[i - 1] / edx1;
    }
    Sp1[(n_spot - 1) / 2] = S1;


    Sp2[0] = S2 * std::exp((n_spot - 1) / 2 * dx2); // 기초 자산의 노드별 가격 설정
    for (i = 1; i < n_spot; ++i) {
        Sp2[i] = Sp2[i - 1] / edx2;
    }
    Sp2[(n_spot - 1) / 2] = S2;

    iop = -1;
    if (option_type == "Call" || option_type == "call") {
        iop = 1;
    }

    double pu1, pm1, pd1, pmm1, pu2, pm2, pd2, pmm2;

    pu1 = -0.5 * (mu1 / dx1 + sigma1 * sigma1 / (dx1 * dx1)) * dt;
    pm1 = 1.0 + sigma1 * sigma1 * dt / (dx1 * dx1) + 0.5 * r * dt;
    pd1 = 0.5 * (mu1 / dx1 - sigma1 * sigma1 / (dx1 * dx1)) * dt;
    pmm1 = 0.125 * rho * sigma1 * sigma2 * dt / (dx1 * dx2);


    pu2 = -0.5 * (mu2 / dx2 + sigma2 * sigma2 / (dx2 * dx2)) * dt;
    pm2 = 1.0 + sigma2 * sigma2 * dt / (dx2 * dx2) + 0.5 * r * dt;
    pd2 = 0.5 * (mu2 / dx2 - sigma2 * sigma2 / (dx2 * dx2)) * dt;
    pmm2 = pmm1;

    // n_spot 개 미지수중에 2개의 경계조선을 적용하므로 실제 미지수는 n_spot-2개 가 된다.
    double **smatrix1, **smatrix2, **tmp, *known_value, *unknown_value;

    smatrix1 = new double *[n_spot - 2];  // row create
    smatrix2 = new double *[n_spot - 2];  // row create
    tmp = new double *[n_spot - 2]; // row create
    known_value = new double[n_spot - 2];
    unknown_value = new double[n_spot - 2];

    for (i = 0; i < n_spot - 2; ++i) {
        smatrix1[i] = new double[n_spot - 2]; // column create
        smatrix2[i] = new double[n_spot - 2]; // column create
        tmp[i] = new double[n_spot - 2]; // column create
    }


    for (i = 0; i < n_spot - 2; ++i) {
        for (j = 0; j < n_spot - 2; ++j) {
            smatrix1[i][j] = 0.0;
            smatrix2[i][j] = 0.0;
        }
    }


    // 0~nnodel-3는 nnodel-2 개
    // 1 자산이 미지수 이므로 1자산의 격자수
    smatrix1[0][0] = pm1 + 2.0 * pu1; // Gamma = 0 경계 조건
    smatrix1[0][1] = pd1 - pu1; // Gamma = 0 경계 조건

    for (i = 1; i < n_spot - 3; ++i) {
        smatrix1[i][i - 1] = pu1;
        smatrix1[i][i] = pm1;
        smatrix1[i][i + 1] = pd1;
    }
    smatrix1[n_spot - 3][n_spot - 4] = pu1 - pd1;// Gamma = 0 경계 조건
    smatrix1[n_spot - 3][n_spot - 3] = pm1 + 2.0 * pd1;// Gamma = 0 경계 조건

    // 2 자산이 미지수 이므로 2자산의 격자수
    smatrix2[0][0] = pm2 + 2.0 * pu2; // Gamma = 0 경계 조건
    smatrix2[0][1] = pd2 - pu2; // Gamma = 0 경계 조건

    for (i = 1; i < n_spot - 3; ++i) {
        smatrix2[i][i - 1] = pu2;
        smatrix2[i][i] = pm2;
        smatrix2[i][i + 1] = pd2;
    }

    smatrix2[n_spot - 3][n_spot - 4] = pu2 - pd2;// Gamma = 0 경계 조건
    smatrix2[n_spot - 3][n_spot - 3] = pm2 + 2.0 * pd2;// Gamma = 0 경계 조건


    // 만기시점에서의 경계조건 부여
    // 만기시점 옵션
    for (i = 0; i < n_spot; ++i) {
        for (j = 0; j < n_spot; ++j) {
            oldCV[i][j] = MAX(iop * ((Sp1[i] - Sp2[j]) - X), 0.0);
        }
    }

    // 만기에서의 Backward Induction 실행하여 현재가치 계산
    for (t = n_step - 1; t >= 0; --t) {
        // 1st phase
        for (i = 1; i < n_spot - 1; ++i) {
            for (j = 0; j < n_spot - 2; ++j) {
                for (k = 0; k < n_spot - 2; ++k) {
                    tmp[j][k] = smatrix1[j][k];
                }

                known_value[j] =
                        oldCV[j + 1][i] + pmm1 *
                                          (oldCV[j][i - 1] + oldCV[j + 2][i + 1]
                                           - oldCV[j][i + 1] - oldCV[j + 2][i - 1]);
            }

            tridiagonal_elimination(tmp, known_value, unknown_value, n_spot - 2);
            for (j = 0; j < n_spot - 2; ++j) {
                newCV[j + 1][i] = unknown_value[i];
            }
        }

        // 경계 조선을 적용한 값 대입
        for (i = 0; i < n_spot; ++i) {
            newCV[i][0] = (2.0 * newCV[i][1] - newCV[i][2]);
            newCV[i][n_spot - 1] = (2.0 * newCV[i][n_spot - 2] - newCV[i][n_spot - 3]);
        }

        // 2nd phase
        for (i = 1; i < n_spot - 1; ++i) {
            for (j = 0; j < n_spot - 2; ++j) {
                for (k = 0; k < n_spot - 2; ++k) {
                    tmp[j][k] = smatrix2[j][k];
                }
                known_value[j] =
                        newCV[i][j + 1] + pmm2 *
                                          (newCV[i - 1][j] + newCV[i + 1][j + 2]
                                           - newCV[i + 1][j] - newCV[i - 1][j + 2]);
            }
            tridiagonal_elimination(tmp, known_value, unknown_value, n_spot - 2);
            for (j = 0; j < n_spot - 2; ++j) {
                oldCV[i][j + 1] = unknown_value[j];
            }
        }

        for (i = 0; i < n_spot; ++i) {
            oldCV[0][i] = (2.0 * oldCV[1][i] - oldCV[2][i]);
            oldCV[n_spot - 1][i] = (2.0 * oldCV[n_spot - 2][i] - oldCV[n_spot - 3][i]);
        }
    }

    m_price = oldCV[(n_spot - 1) / 2][(n_spot - 1) / 2];
    for (i = 0; i < n_spot; ++i) {
        delete[] oldCV[i];
        delete[] newCV[i];
    }

    delete[]newCV;
    delete[]oldCV;
    delete[]Sp2;
    delete[]Sp1;

    for (i = 0; i < n_spot - 2; ++i) {
        delete[] smatrix1[i];
        delete[] smatrix2[i];
        delete[] tmp[i];
    }

    delete[] smatrix1;
    delete[] smatrix2;
    delete[] tmp;
    delete[] known_value;
    delete[] unknown_value;
}

void CPrice::simulation_european_option_price(CIndex &index, CYield &yield, CProduct &product, unsigned long n_sim,
                                              double *rn) {
    std::string option_type = product.m_option_type;
    double S = index.m_spot;
    double X = product.m_strike;
    double r = yield.m_riskfree;
    double q = index.m_dividend;
    double sigma = index.m_vol;
    double T = product.m_maturity;

    unsigned long i;
    int iop;
    double St, price, mudt, ssqrtdt;
    double *pv;
    pv = new double[n_sim];

    price = 0.0;
    mudt = (r - q - 0.5 * sigma * sigma) * T;
    ssqrtdt = sigma * std::sqrt(T);
    iop = -1;

    if (option_type == "Call" || option_type == "call") {
        iop = 1;
    }

    for (i = 0; i < n_sim; ++i) {
        rn[i] = normdistrand_BoxMuller();
        St = S * std::exp(mudt + ssqrtdt * rn[i]); // 주가생성
        pv[i] = MAX(iop * (St - X), 0.0); // 만기 payoff
        price += pv[i];
    }

    m_price = (price / n_sim) * std::exp(-r * T); // 만기 payoff의 평균을 현가
    mean_stddev_error(n_sim, pv);

    delete[] pv;
}

void CPrice::antithetic_variation_simulation_european_option_space(CIndex &index, CYield &yield, CProduct &product,
                                                                   unsigned long n_sim, double *rn) {
    std::string option_type = product.m_option_type;
    double S = index.m_spot;
    double X = product.m_strike;
    double r = yield.m_riskfree;
    double q = index.m_dividend;
    double sigma = index.m_vol;
    double T = product.m_maturity;

    unsigned long i;
    int iop;
    double St, price, mudt, ssqrtdt, prn, nrn;
    double *pv;

    pv = new double[n_sim];
    price = 0.0;
    mudt = (r - q - 0.5 * sigma * sigma) * T;
    ssqrtdt = sigma * std::sqrt(T);
    iop = -1;

    if (option_type == "Call" || option_type == "call") {
        iop = 1;
    }

    for (i = 0; i < n_sim; ++i) {
        St = S * std::exp(mudt + ssqrtdt * rn[i]); // 주가생성
        prn = MAX(iop * (St - X), 0.0); // 만기 payoff

        St = S * std::exp(mudt - ssqrtdt * rn[i]); // 주가생성
        nrn = MAX(iop * (St - X), 0.0); // 만기 payoff

        pv[i] = (prn + nrn) / 2.0;
        price += pv[i];
    }

    m_price = (price / n_sim) * std::exp(-r * T); // 만기 payoff의 평균을 현가
    mean_stddev_error(n_sim, pv);

    delete[] pv;

}

void CPrice::control_variation_simulation_european_option_space(CIndex &index, CYield &yield, CProduct &product,
                                                                unsigned long n_sim, double *rn) {
    std::string option_type = product.m_option_type;
    double S = index.m_spot;
    double X = product.m_strike;
    double r = yield.m_riskfree;
    double q = index.m_dividend;
    double sigma = index.m_vol;
    double T = product.m_maturity;

    unsigned long i;
    int iop;
    double price, mudt, ssqrtdt, meanS, meanP;
    double *pv, *sv;
    pv = new double[n_sim];
    sv = new double[n_sim];

    price = 0.0;
    meanP = 0.0;
    meanS = 0.0;

    mudt = (r - q - 0.5 * sigma * sigma) * T;
    ssqrtdt = sigma * std::sqrt(T);
    iop = -1;

    if (option_type == "Call" || option_type == "call") {
        iop = 1;
    }

    for (i = 0; i < n_sim; ++i) {
        sv[i] = S * std::exp(mudt + ssqrtdt * rn[i]); // 주가생성
        pv[i] = MAX(iop * (sv[i] - X), 0.0); // 만기 payoff
        meanS += sv[i];
        meanP += pv[i];
    }

    meanS /= n_sim;
    meanP /= n_sim;

    double variance, covariance, beta, ES;
    variance = covariance = 0.0;
    for (i = 0; i < n_sim; ++i) {
        variance += (sv[i] - meanS) * (sv[i] - meanS);
        covariance += (sv[i] - meanS) * (pv[i] - meanP);
    }

    variance /= (n_sim - 1);
    covariance /= (n_sim - 1);

    beta = covariance / variance;
    ES = S * std::exp((r - q) * T);

    for (i = 0; i < n_sim; ++i) {
        pv[i] = pv[i] + beta * (ES - sv[i]);
        price += pv[i];
    }

    m_price = (price / n_sim) * std::exp(-r * T); // 만기 payoff의 평균을 현가
    mean_stddev_error(n_sim, pv);

    delete[] pv;
    delete[] sv;
}

void CPrice::importance_sampling_simulation_european_option_space(CIndex &index, CYield &yield, CProduct &product,
                                                                  unsigned long n_sim, double *rn) {
    std::string option_type = product.m_option_type;
    double S = index.m_spot;
    double X = product.m_strike;
    double r = yield.m_riskfree;
    double q = index.m_dividend;
    double sigma = index.m_vol;
    double T = product.m_maturity;

    unsigned long i;
    int iop;
    double St, price, mudt, ssqrtdt, c, likeR;
    double *pv;

    pv = new double[n_sim];
    price = 0.0;
    c = 0.4; // 0.1~0.7에서 선택
    mudt = (r + c - q - 0.5 * sigma * sigma) * T;
    ssqrtdt = sigma * std::sqrt(T);
    iop = -1;

    if (option_type == "Call" || option_type == "call") {
        iop = 1;
    }

    for (i = 0; i < n_sim; ++i) {
        St = S * std::exp(mudt + ssqrtdt * rn[i]); // 주가생성
        likeR = std::exp(-0.5 * c * c * T / (sigma * sigma) - c * rn[i] / sigma); // likelihood ratio

        pv[i] = MAX(iop * (St - X), 0.0) * likeR; // 만기 payoff
        price += pv[i];
    }

    m_price = (price / n_sim) * std::exp(-r * T); // 만기 payoff의 평균을 현가
    mean_stddev_error(n_sim, pv);
    delete[] pv;
}

void CPrice::halton_sequence_simulation_european_option_space(CIndex &index, CYield &yield, CProduct &product,
                                                              unsigned long n_sim) {
    std::string option_type = product.m_option_type;
    double S = index.m_spot;
    double X = product.m_strike;
    double r = yield.m_riskfree;
    double q = index.m_dividend;
    double sigma = index.m_vol;
    double T = product.m_maturity;

    unsigned long i;
    int iop;
    double St, price, mudt, ssqrtdt;
    double *pv, *rn;

    unsigned D = 1;
    unsigned N = n_sim;

    rn = new double[N];
    pv = new double[N];

    unsigned x;
    auto *halton = new double[N];

    for (x = 0; x < N; ++x) {
        halton[x] = halton_sequence_(x, 2);
        rn[x] = inverse_normal_cumulative_distribution_function(halton[x]);
    }

    price = 0.0;
    mudt = (r - q - 0.5 * sigma * sigma) * T;
    ssqrtdt = sigma * std::sqrt(T);

    iop = -1;
    if (option_type == "Call" || option_type == "call") {
        iop = 1;
    }

    for (i = 0; i < N; ++i) {
        St = S * std::exp(mudt + ssqrtdt * rn[i]); // 주가생성
        pv[i] = MAX(iop * (St - X), 0.0); // 만기 payoff
        price += pv[i];
    }

    m_price = (price / N) * std::exp(-r * T); // 만기 payoff의 평균을 현가

    mean_stddev_error(n_sim, pv);

    delete[]pv;
    delete[]rn;
    delete[]halton;
}

void CPrice::sobol_sequence_simulation_european_option_space(CIndex &index, CYield &yield, CProduct &product,
                                                             unsigned long n_sim) {
    std::string option_type = product.m_option_type;
    double S = index.m_spot;
    double X = product.m_strike;
    double r = yield.m_riskfree;
    double q = index.m_dividend;
    double sigma = index.m_vol;
    double T = product.m_maturity;

    unsigned long i;
    int iop;
    double St, price, mudt, ssqrtdt;
    double *pv, *rn;

    unsigned N;

    for (i = 0; i < n_sim; ++i) {
        N = static_cast<int>(std::pow(2.0, static_cast<double>(i))) - 1; // 2^n -1가 인상적 ..1023, 2047, 4095, ...
        if (N > n_sim) break;
    }

    unsigned D = 1; // 2^n -1가 인상적

    rn = new double[N];
    pv = new double[N];

    unsigned x;

    double **sobol = sobol_points(N, D, "new-joe-kuo-6.21201.txt");

    for (x = 0; x < N; ++x) {
        rn[x] = inverse_normal_cumulative_distribution_function(sobol[x][0]);
    }

    price = 0.0;
    mudt = (r - q - 0.5 * sigma * sigma) * T;
    ssqrtdt = sigma * std::sqrt(T);

    iop = -1;
    if (option_type == "Call" || option_type == "call") {
        iop = 1;
    }
    for (i = 0; i < N; ++i) {
        St = S * std::exp(mudt + ssqrtdt * rn[i]); // 주가생성
        pv[i] = MAX(iop * (St - X), 0.0); // 만기 payoff
        price += pv[i];
    }

    m_price = (price / N) * std::exp(-r * T); // 만기 payoff의 평균을 현가

    mean_stddev_error(n_sim, pv);
    for (x = 0; x < N; ++x) {
        delete[] sobol[x];
    }

    delete[]pv;
    delete[]rn;
    delete[]sobol;

}

void CPrice::ci_binomial_tree_european_option_space(CIndex &index, CYield &yield, CProduct &product, int n_step,
                                                    double alpha) {
    std::string option_type = product.m_option_type;
    double S = index.m_spot;
    double X = product.m_strike;
    double r = yield.m_riskfree;
    double q = index.m_dividend;
    double sigma = index.m_vol;
    double T = product.m_maturity;

    int i, j, iop, ccnode;
    double dx, nu, dfactor, Pu, Pd, dt;
    double smax, smin, nodemax, nodemin;

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

    iop = -1;
    if (option_type == "Call" || option_type == "call") {
        iop = 1;
    }

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

    m_price = OCV[0]; // 옵션의 현재가치

    delete[]Sp;
    delete[]OCV;
    delete[]NCV;

}

void CPrice::ci_trinomial_tree_european_option_space(CIndex &index, CYield &yield, CProduct &product, int n_step,
                                                     double alpha) {
    std::string option_type = product.m_option_type;
    double S = index.m_spot;
    double X = product.m_strike;
    double r = yield.m_riskfree;
    double q = index.m_dividend;
    double sigma = index.m_vol;
    double T = product.m_maturity;

    int i, j, iop, ccnode, nstep;
    double dx, nu, dfactor, Pu, Pd, dt, Pm;
    double smax, smin, nodemax, nodemin;

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

    for (i = 0; i < nstep; ++i) { Sp[i] = S * std::exp((ccnode - i) * dx); } // Tree 주가 생성

    iop = -1;
    if (option_type == "Call" || option_type == "call") { iop = 1; }

    for (i = 0; i < nstep; ++i) {
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

    m_price = OCV[0]; // 옵션의 현재가치

    delete[]Sp;
    delete[]OCV;
    delete[]NCV;
}

void
CPrice::ci_implicit_fdm_european_option_space(CIndex &index, CYield &yield, CProduct &product, int n_step, int n_spot,
                                              double alpha) {
    std::string option_type = product.m_option_type;
    double S = index.m_spot;
    double X = product.m_strike;
    double r = yield.m_riskfree;
    double q = index.m_dividend;
    double sigma = index.m_vol;
    double T = product.m_maturity;

    int i, j, iop, ccnode, nspot, t;
    double dx, nu, Pu, Pd, dt, Pm;
    double smax, smin, nodemax, nodemin;

    dt = T / n_step;
    dx = sigma * std::sqrt(3.0 * dt);

    nu = r - q - 0.5 * sigma * sigma;
    Pu = -0.5 * dt *(sigma * sigma / (dx * dx) + nu / dx);
    Pm = 1.0 + dt * (sigma * sigma / (dx * dx) + r);
    Pd = -0.5 * dt *(sigma * sigma / (dx * dx) - nu / dx);

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

    iop = -1;
    if (option_type == "Call" || option_type == "call") { iop = 1; }

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
    for (t = n_step - 1; t > 0; --t) {
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

    m_price = CV[(nspot - 1) / 2];
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
}
