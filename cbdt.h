//
// Created by jewoo on 2022-11-02.
//

#pragma once

class CCBDT {
public:
    CCBDT();

    CCBDT(int n);

    CCBDT(double bond_maturity, double dt);

    ~CCBDT();

public:

    int m_nnode{};    // bdt tree의 node 수
    double m_dt{};    // tree의 간격
    double m_pu{};    // up 확률
    double m_pd{};    // down 확률

    double *m_time{}; // node time
    double *m_rate{}; // zero rate
    double *m_vol{};  // volatility

    double **m_ndf{}; // node discount factor
    double **m_nQ{};  // 상태가격 (할인된 확률)
    double **m_nsrate{};  // short rate tree
    double *m_medianU{};  // short rate의 추정된 중위수
    double **m_nQd{}; // 하락상태의 상태가격
    double *m_nvol{}; // node volatility

    bool m_bvolfix{}; // 변동성 고정 calibration이면 true 아니면 false

public:
    // 변동성 고려하여 newton raphson 법으로 BDT tree 구하기
    void build_bdt_tree_vol_curve(int nrate, double *rtime,
                                  double *rate, int nvol, double *vtime, double *vol);

};
