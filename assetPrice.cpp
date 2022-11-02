//
// Created by jewoo on 2022-11-02.
//
#include <iostream>
#include <cmath>
#include "bdt.h"
#include "cbdt.h"
#include "assetPrice.h"
#include "mathlib.h"

AssetCPrice::AssetCPrice() {
    m_coupon_frequency = 1;
    m_n_coupon = 1;
    m_notional = 100;
    m_coupon_rate = 0.0;
    m_bond_maturity = 5.0;
}

AssetCPrice::~AssetCPrice() = default;

AssetCPrice::AssetCPrice(double notional, double coupon_rate, double bond_maturity, int coupon_frequency) {
    m_coupon_frequency = coupon_frequency;
    m_notional = notional;
    m_coupon_rate = coupon_rate;
    m_bond_maturity = bond_maturity;
    m_n_coupon = static_cast<int>(m_bond_maturity * m_coupon_frequency);
}

// 채권 가격 계산
double AssetCPrice::zero_coupon_bond_price(CBDT &bdt) {
    int i, j, nnode;
    double price;
    double *bprice;

    nnode = static_cast<int>(m_bond_maturity / bdt.m_dt);
    if (nnode > bdt.m_nnode) {
        std::cout << "BDT Tree를 보다 조밀하게 생성하세요" << std::endl;
        exit(0);
    }
    bprice = new double[nnode + 1];
    for (j = 0; j < nnode + 1; ++j) {
        bprice[j] = m_notional; // 채권 만기시 원금 지급
    }

    double pd, pu;
    pu = bdt.m_pu;
    pd = 1.0 - pu;

    for (i = nnode - 1; i >= 0; --i) {
        for (j = 0; j <= i; ++j) {
            bprice[j] = (pu * bprice[j] + pd * bprice[j + 1]) * bdt.m_ndf[i][j]; // 채권가격 할인
        }
    }
    price = bprice[0];

    delete[] bprice;
    std::cout << "Zero coupon bond price : " << price << std::endl;
    return price;
}

double AssetCPrice::zero_coupon_bond_price_with_vol(CCBDT &bdt) const {
    int i, j, nnode;
    double price;
    double *bprice;

    nnode = static_cast<int>(m_bond_maturity / bdt.m_dt);
    if (nnode > bdt.m_nnode) {
        std::cout << "BDT Tree를 보다 조밀하게 생성하세요" << std::endl;
        exit(0);
    }
    bprice = new double[nnode + 1];
    for (j = 0; j < nnode + 1; ++j) {
        bprice[j] = m_notional; // 채권 만기시 원금 지급
    }

    double pd, pu;
    pu = bdt.m_pu;
    pd = 1.0 - pu;

    for (i = nnode - 1; i >= 0; --i) {
        for (j = 0; j <= i; ++j) {
            bprice[j] = (pu * bprice[j] + pd * bprice[j + 1]) * bdt.m_ndf[i][j]; // 채권가격 할인
        }
    }
    price = bprice[0];

    delete[] bprice;
    std::cout << "Zero coupon bond price : " << price << std::endl;
    return price;
}

