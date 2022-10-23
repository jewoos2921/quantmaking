//
// Created by jewoo on 2022-10-05.
//

#pragma once

class CYield {
public:
    double m_riskfree;
public:
    CYield();

    CYield(double riskfree);

    ~CYield();

    void get_riskfree();
};