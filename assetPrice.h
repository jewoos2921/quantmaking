//
// Created by jewoo on 2022-11-02.
//

#pragma once

class CBDT;

class CCBDT;

class AssetCPrice {
public:
    AssetCPrice();

    ~AssetCPrice();

    AssetCPrice(double notional, double coupon_rate, double bond_maturity,
                int coupon_frequency);

public:
    int m_coupon_frequency; // 쿠폰 지급 주기
    int m_n_coupon; // 쿠폰 수량
    double m_notional; // 원금
    double m_coupon_rate; // 쿠폰이자율
    double m_bond_maturity; // 채권 만기

public:
    double zero_coupon_bond_price(CBDT &bdt); // 채권 가격 계산

    double zero_coupon_bond_price_with_vol(CCBDT &bdt) const;
};