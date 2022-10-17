//
// Created by jewoo on 2022-10-05.
//
#include <iostream>
#include "yield.h"

CYield::CYield() {
    std::cout << "CYield 디폴트 생성자" << std::endl;
    m_riskfree = 0.05;
}

CYield::~CYield() {
    std::cout << "CYield 소멸자" << std::endl;
}

CYield::CYield(double riskfree) {
    std::cout << "CYield 생성자" << std::endl;
    m_riskfree = riskfree;
}
