//
// Created by jewoo on 2022-10-05.
//

#pragma once

class CYield;

class CIndex;

class COption {
public:
    CIndex *index;
    CYield *yield;
    double m_strike;
    double m_maturity;
    double m_optionprice;

public:

    COption();

    COption(double strike, double maturity);

    ~COption();

    double N(double z);

    void european_calloption_price(CIndex &cIndex, CYield &cYield);

};