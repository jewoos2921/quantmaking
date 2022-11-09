//
// Created by jewoo on 2022-11-09.
//

#pragma once

class CHWTree {
public:
    CHWTree();

    ~CHWTree();

public:
    int m_nnode; // tree 단계수
    int m_jmax; // Jmax 값
    double m_alpha; // mean reversion speed
    double m_sigma; // short rate의 변동성
    double m_dt; // tree 1 단계의 시간 간격

    double *m_pu{}; // probability of up
    double *m_pm{}; // probability of middle
    double *m_pd{}; // probability of down

    double **m_ndf{}; // node discount factor
    double **m_nrate{}; // node short rate

public:
    // Normal Distribution Trinimial Tree 생성
    void build_hull_white_short_rate_tree(int nrtime, double *time, double *rate);
};