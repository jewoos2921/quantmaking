//
// Created by jewoo on 2022-10-05.
//
#include <iostream>
#include <cmath>
#include "option.h"
#include "yield.h"
#include "index.h"

#define PI 3.14159265358979323846264338327950288419716939937510582097

void COption::get_strike() {
    std::cout << "행사 가격: ";
    std::cin >> m_strike;
}

void COption::get_maturity() {
    std::cout << "잔존 만기(년): ";
    std::cin >> m_maturity;
}

double COption::N(double z) {
    double a1, a2, a3, a4, a5, r, c;

    a1 = 0.31938153;
    a2 = -0.356563782;
    a3 = 1.781477937;
    a4 = -1.821255978;
    a5 = 1.330274429;
    r = 0.2316419;
    c = 1 / std::sqrt(2 * PI);

    if (z > 100.0) { return 1.0; }
    else if (z < -100.0) { return 0.0; }


    double x = std::fabs(z);
    double k = 1.0 / (1.0 + r * x);
    double b = c * std::exp((-z * z) / 2.0);
    double nv = ((((a5 * k + a4) * k + a3) * k + a2) * k + a1) * k;
    nv = 1.0 - b * nv;
    if (z < 0.0) { nv = 1.0 - nv; }
    return nv;
}

void COption::european_calloption_price(CIndex &cIndex, CYield &cYield) {
    double d1, d2, spot, vol, dividend, riskfree;

    spot = cIndex.m_spot;
    riskfree = cYield.m_riskfree;
    dividend = cIndex.m_dividend;
    vol = cIndex.m_vol;

    d1 = (std::log(spot / m_strike) + (riskfree - dividend + (vol * vol) / 2.0) * m_maturity)
         / (vol * std::sqrt(m_maturity));
    d2 = d1 - vol * std::sqrt(m_maturity);
    m_optionprice = spot * std::exp(-dividend * m_maturity) * N(d1)
                    - m_strike * std::exp(-riskfree * m_maturity) * N(d2);

}
