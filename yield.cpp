//
// Created by jewoo on 2022-10-05.
//
#include <iostream>
#include "yield.h"

void CYield::get_riskfree() {
    std::cout << "무위험 이자율(%): ";
    std::cin >> m_riskfree;
    m_riskfree /= 100.0;
}
