//
// Created by jewoo on 2022-10-21.
//
#include <iostream>
#include <utility>
#include "product.h"

CProduct::CProduct() {
    m_strike = 105.0;
    m_maturity = 1.0;
    m_option_type = "call";
    m_correlation = 0.5;
}

CProduct::CProduct(std::string option_type, double strike, double maturity) {
    m_strike = strike;
    m_maturity = maturity;
    m_option_type = std::move(option_type);
}

CProduct::~CProduct() {

}

void CProduct::get_strike() {
    std::cout << "행가 가격 : ";
    std::cin >> m_strike;
}

void CProduct::get_maturity() {
    std::cout << "잔존 만기(년) : ";
    std::cin >> m_maturity;
}

CProduct::CProduct(std::string option_type, double strike, double maturity, double correlation) {
    m_strike = strike;
    m_maturity = maturity;
    m_option_type = std::move(option_type);
    m_correlation = correlation;
}
