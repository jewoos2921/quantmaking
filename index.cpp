//
// Created by jewoo on 2022-10-05.
//

#include <iostream>
#include "index.h"


CIndex::CIndex() {
    std::cout << "CIndex 디폴트 생성자" << std::endl;
    m_spot = 100.0;
    m_vol = 0.2;
    m_dividend = 0.03;
}

CIndex::~CIndex() {
    std::cout << "CIndex 소멸자" << std::endl;
}

CIndex::CIndex(double spot, double vol, double dividend) {
    std::cout << "CIndex 생성자" << std::endl;
    m_spot = spot;
    m_vol =vol;
    m_dividend = dividend;
}
