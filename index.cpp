//
// Created by jewoo on 2022-10-05.
//

#include <iostream>
#include "index.h"


CIndex::CIndex() {
//    std::cout << "CIndex 디폴트 생성자" << std::endl;
    m_spot = 100.0;
    m_vol = 0.2;
    m_dividend = 0.03;
}

CIndex::~CIndex() {
//    std::cout << "CIndex 소멸자" << std::endl;
}

CIndex::CIndex(double spot, double vol, double dividend) {
//    std::cout << "CIndex 생성자" << std::endl;
    m_spot = spot;
    m_vol = vol;
    m_dividend = dividend;
}

void CIndex::get_spot() {
    std::cout << "기초 자산 가격: ";
    std::cin >> m_spot;
}

void CIndex::get_Vol() {
    std::cout << "변동성(%): ";
    std::cin >> m_vol;
    m_vol /= 100.0;
}

void CIndex::get_dividend() {
    std::cout << "배당률(%): ";
    std::cin >> m_dividend;
    m_dividend /= 100.0;
}

