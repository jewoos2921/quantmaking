//
// Created by jewoo on 2022-10-05.
//

#pragma once

class CIndex {
public:
    CIndex();

    CIndex(double spot, double vol, double dividend);

    ~CIndex();

public:
    double m_spot;
    double m_vol;
    double m_dividend;
};