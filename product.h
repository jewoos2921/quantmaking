//
// Created by jewoo on 2022-10-21.
//

#pragma once

#include <string>

class CProduct {
public:
    double m_strike;
    double m_maturity;
    std::string m_option_type;

public:
    CProduct();

    CProduct(std::string option_type,
             double strike, double maturity);

    ~CProduct();

    void get_strike();

    void get_maturity();

};
