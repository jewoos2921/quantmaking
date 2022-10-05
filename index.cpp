//
// Created by jewoo on 2022-10-05.
//

#include <iostream>
#include "index.h"

void CIndex::get_spot() {
    std::cout << "기초 자산 가격: ";
    std::cin >> m_spot;
}

void CIndex::get_dividend() {
    std::cout << "배당률(%): ";
    std::cin >> m_dividend;
    m_dividend /= 100.0;
}

void CIndex::get_vol() {
    std::cout << "변동성(%): ";
    std::cin >> m_vol;
    m_vol /= 100.0;
}
