//
// Created by jewoo on 2022-11-11.
//

#pragma once

class CCBDT;

class CPriceBlack {
public:
    CPriceBlack();

    ~CPriceBlack();

    CPriceBlack(double notional, double coupon_rate, double bond_maturity,
                int coupon_frequency);

public:
    int m_coupon_frequency; // 쿠폰 지급 주기
    int m_ncoupon;  // 쿠폰 수량

    double m_notioanl;  // 원금
    double m_coupon_rate;   // 쿠폰이자율
    double m_bond_maturity; // 채권만기

    int m_option_position;  // 1:콜, -1:풋
    double m_strike;    // 옵션 행사 가격
    double m_option_maturity;   // 옵션 만기

public:
    double coupon_bond_price(CCBDT &bdt);   // 채권가격 계산

    double coupon_bond_option_price(CCBDT &bdt);    // 옵션가격

    double callable_coupon_bond_price(CCBDT &bdt);  // 옵션부채권가격

};