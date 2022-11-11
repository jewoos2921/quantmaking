//
// Created by jewoo on 2022-11-09.
//

#pragma once

class CHWTree;

class CHWTPrice {
public:
    CHWTPrice();

    CHWTPrice(double notional, double coupon_rate, double bond_maturity,
              int coupon_frequency);

    ~CHWTPrice();

public:
    int m_coupon_frequency; // 쿠폰지급주기
    int m_ncoupon; // 쿠폰 수량

    double m_notional; // 원금
    double m_coupon_rate; // 쿠폰 이자율
    double m_bond_maturity; // 채권만기

public:
    double zero_coupon_bond_price(CHWTree &hw); // 채권가격 계산

    double zero_coupon_bond_price_black_karasinski_short_rate_tree(CHWTree &hw);
};