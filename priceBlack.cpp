//
// Created by jewoo on 2022-11-11.
//

#include <iostream>
#include <cmath>
#include "cbdt.h"
#include "priceBlack.h"
#include "mathlib.h"

CPriceBlack::CPriceBlack() {
    m_coupon_frequency = 1;
    m_ncoupon = 1;
    m_notioanl = 100.0;
    m_coupon_rate = 0.05;
    m_bond_maturity = 5.0;

    m_option_position = 1;
    m_strike = 98.0;
    m_option_maturity = 3.0;
}

CPriceBlack::~CPriceBlack() = default;

CPriceBlack::CPriceBlack(double notional, double coupon_rate, double bond_maturity, int coupon_frequency) {
    m_coupon_frequency = coupon_frequency;
    m_notioanl = notional;
    m_coupon_rate = coupon_rate;
    m_bond_maturity = bond_maturity;
    m_ncoupon = static_cast<int>(m_bond_maturity * m_coupon_frequency);
}

double CPriceBlack::coupon_bond_price(CCBDT &bdt) {
    int i, j, nnode, ithc;
    double price, coupon;
    double *bprice, *coupon_time;

    nnode = static_cast<int>(m_bond_maturity / bdt.m_dt);
    if (nnode > bdt.m_nnode) {
        std::cout << "BDT 트리를 보다 조밀하게 생성하세요" << std::endl;
        exit(0);
    }

    coupon = m_notioanl * m_coupon_rate / m_coupon_frequency;
    bprice = new double[nnode + 1];
    coupon_time = new double[m_ncoupon];

    for (i = 0; i < m_ncoupon; ++i) {
        coupon_time[i] = static_cast<double>(i + 1) / m_coupon_frequency;
    }

    for (j = 0; j < nnode + 1; ++j) {
        bprice[j] = m_notioanl + coupon; // 채권만기시 원금 지급
    }

    double pu, pd;
    pu = bdt.m_pu;
    pd = 1.0 - pu;
    ithc = m_ncoupon - 2;


    for (i = nnode - 1; i >= 0; --i) {
        for (j = 0; j <= i; ++j) {
            // 채권가격 할인
            bprice[j] = (pu * bprice[j] + pd * bprice[j + 1]) * bdt.m_ndf[i][j];
        }
        if (coupon_time[ithc] == i * bdt.m_dt) {
            for (j = 0; j <= i; ++j) {
                bprice[j] += coupon; // 이표지급
            }
            ithc--;
        }
    }

    price = bprice[0];

    delete[] bprice;
    delete[] coupon_time;

    std::cout << "coupon bond price: " << price << std::endl;
    return price;
}

double CPriceBlack::coupon_bond_option_price(CCBDT &bdt) {
    int i, j, nonode, ithc, nbnode;
    double price, coupon;
    double *bprice, *coupon_time, *oprice;

    nbnode = static_cast<int>(m_bond_maturity / bdt.m_dt);
    nonode = static_cast<int>(m_option_maturity / bdt.m_dt);

    if (m_option_maturity > m_bond_maturity) {
        std::cout << "옵션만기가 채권만기보다 깁니다." << std::endl;
        exit(0);
    }

    if (nbnode > bdt.m_nnode) {
        std::cout << "BDT 트리를 보다 조밀하게 생성하세요" << std::endl;
        exit(0);
    }

    coupon = m_notioanl * m_coupon_rate / m_coupon_frequency;
    bprice = new double[nbnode + 1];
    coupon_time = new double[m_ncoupon];


    for (i = 0; i < m_ncoupon; ++i) {
        coupon_time[i] = static_cast<double>(i + 1) / m_coupon_frequency;
    }

    for (j = 0; j < nbnode + 1; ++j) {
        bprice[j] = m_notioanl + coupon; // 채권만기시 원금 지급
    }

    double pu, pd;
    pu = bdt.m_pu;
    pd = 1.0 - pu;
    ithc = m_ncoupon - 2;


    for (i = nbnode - 1; i >= nonode; --i) {
        for (j = 0; j <= i; ++j) {
            // 채권가격 할인
            bprice[j] = (pu * bprice[j] + pd * bprice[j + 1]) * bdt.m_ndf[i][j];
        }

        if (coupon_time[ithc] == i * bdt.m_dt) {
            for (j = 0; j <= i; ++j) {
                bprice[j] += coupon; // 이표지급
            }
            ithc--;
        }
    }

    oprice = new double[nonode + 1]; // 현재시점을 고려하여 +1;
    for (j = 0; j < nonode + 1; ++j) {
        // 옵션만기시점 Payoff
        oprice[j] = MAX(m_option_position * (bprice[j] - m_strike), 0);
    }

    for (i = nbnode - 1; i >= 0; --i) {
        for (j = 0; j <= i; ++j) {
            // 옵션가격 할인
            oprice[j] = (pu * bprice[j] + pd * bprice[j + 1]) * bdt.m_ndf[i][j];
        }
    }

    price = oprice[0];

    delete[] bprice;
    delete[] oprice;
    delete[] coupon_time;

    std::cout << "coupon bond option price: " << price << std::endl;
    return price;
}

double CPriceBlack::callable_coupon_bond_price(CCBDT &bdt) {
    int i, j, nonode, ithc, nbnode;
    double price, coupon;
    double *bprice, *coupon_time;

    nbnode = static_cast<int>(m_bond_maturity / bdt.m_dt);
    nonode = static_cast<int>(m_option_maturity / bdt.m_dt);

    if (m_option_maturity > m_bond_maturity) {
        std::cout << "옵션만기가 채권만기보다 깁니다." << std::endl;
        exit(0);
    }

    if (nbnode > bdt.m_nnode) {
        std::cout << "BDT 트리를 보다 조밀하게 생성하세요" << std::endl;
        exit(0);
    }

    coupon = m_notioanl * m_coupon_rate / m_coupon_frequency;
    bprice = new double[nbnode + 1];
    coupon_time = new double[m_ncoupon];


    for (i = 0; i < m_ncoupon; ++i) {
        coupon_time[i] = static_cast<double>(i + 1) / m_coupon_frequency;
    }

    for (j = 0; j < nbnode + 1; ++j) {
        bprice[j] = m_notioanl + coupon; // 채권만기시 원금 지급
    }

    double pu, pd;
    pu = bdt.m_pu;
    pd = 1.0 - pu;
    ithc = m_ncoupon - 2;


    for (i = nbnode - 1; i >= 0; --i) {
        for (j = 0; j <= i; ++j) {
            // 채권가격 할인
            bprice[j] = (pu * bprice[j] + pd * bprice[j + 1]) * bdt.m_ndf[i][j];
        }

        if (coupon_time[ithc] == i * bdt.m_dt) {
            for (j = 0; j <= i; ++j) {
                bprice[j] += coupon; // 이표지급
            }
            ithc--;
        }

        if (i == nonode) {
            if (m_option_position == 1) { // callbale bond
                for (j = 0; j < nonode + 1; ++j) {
                    bprice[j] = MIN(bprice[j], m_strike); // 옵션 만기 시점 페이오프
                }
            } else {
                for (j = 0; j < nonode + 1; ++j) {
                    // puttable bond
                    bprice[j] = MAX(bprice[j], m_strike);
                }
            }
        }
    }


    price = bprice[0];

    delete[] bprice;

    delete[] coupon_time;

    std::cout << "callalbe coupon bond option price: " << price << std::endl;
    return price;
}
