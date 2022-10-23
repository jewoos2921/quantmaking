//
// Created by jewoo on 2022-10-05.
//
#include <iostream>
#include "yield.h"

CYield::CYield() {
    m_riskfree = 0.05;
}

CYield::~CYield() {

}

CYield::CYield(double riskfree) {
    m_riskfree = riskfree;
}

void CYield::get_riskfree() {
    std::cout << "무위험 이자율(%) : ";
    std::cin >> m_riskfree;
    m_riskfree /= 100.0;
}
