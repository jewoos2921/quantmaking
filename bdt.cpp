//
// Created by jewoo on 2022-11-02.
//
#include <iostream>
#include <cmath>

#include "bdt.h"
#include "mathlib.h"

CBDT::CBDT() {
    m_nnode = 100;
    m_dt = 1.0;
    m_pu = .5;
    m_pd = .5;
    m_bvolfix = true;

}

CBDT::CBDT(int n) {
}

CBDT::CBDT(double bond_maturity, double dt) {
    m_dt = dt;
    m_nnode = static_cast<int>(bond_maturity / dt) + 1; // 채권 만기까지 Tree 생성

    m_pu = 0.5;
    m_pd = 0.5;

    m_ndf = new double *[m_nnode];
    m_nQ = new double *[m_nnode];
    m_nsrate = new double *[m_nnode];
    m_medianU = new double[m_nnode];

    for (int i = 0; i < m_nnode; ++i) {
        m_ndf[i] = new double[m_nnode];
        m_nQ[i] = new double[m_nnode];
        m_nsrate[i] = new double[m_nnode];
    }

    for (int i = 0; i < m_nnode; ++i) {
        for (int j = 0; j < m_nnode; ++j) {
            m_ndf[i][j] = 0.0;
            m_nQ[i][j] = 0.0;
            m_nsrate[i][j] = 0.0;
        }
        m_medianU[i] = 0.0;
    }

    m_bvolfix = true;
}

CBDT::~CBDT() {
    for (int i = 0; i < m_nnode; ++i) {
        delete[] m_ndf[i];
        delete[] m_nQ[i];
        delete[] m_nsrate[i];
    }

    delete[] m_ndf;
    delete[] m_nQ;
    delete[] m_nsrate;
    delete[] m_medianU;
    delete[] m_time;
    delete[] m_rate;
}

void CBDT::build_bdt_tree_fixed_vol(int nrate, double *rtime, double *rate, double vol) {
    std::cout << "변동성을 고정시킨 BDT 모형" << std::endl;
    int i, j, k, iter;
    double x1, x2, f, fx, sqdt, Error, Epsilon, PDB;
    Epsilon = 1.0e-15;

    m_bvolfix = true;

    m_time = new double[m_nnode];
    m_rate = new double[m_nnode];

    for (i = 0; i < m_nnode - 1; ++i) {
        m_time[i] = (i + 1) * m_dt;
        m_rate[i] = linear_interpolation(nrate, rtime, rate, m_time[i]);
    }

    m_medianU[0] = m_rate[0];
    m_nsrate[0][0] = m_medianU[0];
    m_ndf[0][0] = 1.0 / (1 + m_nsrate[0][0] * m_dt);
    m_nQ[0][0] = 1.0;
    sqdt = std::sqrt(m_dt);

    for (i = 1; i < m_nnode - 1; ++i) {
        m_nQ[i][0] = m_pu * m_nQ[i - 1][0] * m_ndf[i - 1][0];
        for (j = 1; j < i; ++j) {
            m_nQ[i][j] = m_pu * m_nQ[i - 1][j - 1] * m_ndf[i - 1][j - 1]
                         + m_pd * m_nQ[i - 1][j] * m_ndf[i - 1][j];
        }

        m_nQ[i][i] = m_pd * m_nQ[i - 1][i - 1] * m_ndf[i - 1][i - 1];

        PDB = 1.0 / std::pow((1.0 + m_rate[i]), m_time[i]);
        x1 = m_medianU[i - 1];
        iter = 0;
        do {
            f = 0.0;
            fx = 0.0;
            k = 0;
            for (j = 0; j <= i; ++j) {
                k = i - 2 * j;
                f += m_nQ[i][j] / (1.0 + x1 * std::exp(vol * k * sqdt) * m_dt);
                fx -= m_nQ[i][j] * std::pow((1.0 + x1 * std::exp(vol * k * sqdt) * m_dt), -2)
                      * std::exp(vol * k * sqdt) * m_dt;
            }
            x2 = x1 - (f - PDB) / (fx);
            Error = std::fabs(x2 - x1);
            x1 = x2;
            iter++;
        } while (Error > Epsilon && iter < 1000);

        m_medianU[i] = x1;
        k = 0;
        if (i < m_nnode) {
            for (j = 0; j <= i; ++j) {
                k = i - 2 * j;
                m_nsrate[i][j] = m_medianU[i] * std::exp(vol * k * sqdt); // SPOT m_nsrate[][]을 적용한 1기간의 할인율 (단리 계산)
                m_ndf[i][j] = 1.0 / (1.0 + m_nsrate[i][j] * m_dt);
            }
        }
    }
}
