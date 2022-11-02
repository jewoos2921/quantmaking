//
// Created by jewoo on 2022-11-02.
//
#include <iostream>
#include <cmath>

#include "cbdt.h"
#include "mathlib.h"

CCBDT::CCBDT() {
    m_nnode = 100;
    m_dt = 1.0;
    m_pu = .5;
    m_pd = .5;
    m_bvolfix = true;

}

CCBDT::CCBDT(int n) {
}

CCBDT::CCBDT(double bond_maturity, double dt) {
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

CCBDT::~CCBDT() {
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

    if (m_bvolfix != true) {
        for (int i = 0; i < m_nnode; ++i) {
            delete[] m_nQd[i];
        }
        delete[] m_nQd;
        delete[] m_vol;
        delete[] m_nvol;
    }
}

void CCBDT::build_bdt_tree_vol_curve(int nrate, double *rtime, double *rate, int nvol, double *vtime, double *vol) {
    std::cout << "변동성 기간구조를 고려한 BDT 모형 " << std::endl;

    int i, j, k, iter;
    double x1, y1, f, fx, sqdt;
    double Error, Epsilon;
    double iv;
    double *Pu, *Pd, PDB;

    Epsilon = 1.0e-15;

    m_bvolfix = false;

    m_nQd = new double *[m_nnode];
    m_nvol = new double[m_nnode];

    for (i = 0; i < m_nnode; ++i) {
        m_nQd[i] = new double[m_nnode];
    }

    for (i = 0; i < m_nnode; ++i) {
        for (j = 0; j < m_nnode; ++j) { m_nQd[i][j] = 0.0; }
    }

    m_time = new double[m_nnode];
    m_rate = new double[m_nnode];
    m_vol = new double[m_nnode];

    for (i = 0; i < m_nnode - 1; ++i) {
        m_time[i] = (i + 1) * m_dt;
        m_rate[i] = linear_interpolation(nrate, rtime, rate, m_time[i]);
        m_vol[i] = linear_interpolation(nvol, vtime, vol, m_time[i]);
    }

    Pu = new double[m_nnode];
    Pd = new double[m_nnode];

    for (i = 0; i < m_nnode; ++i) {
        Pu[i] = 0.0;
        Pd[i] = 0.0;
    }

    m_medianU[0] = m_rate[0];
    m_nsrate[0][0] = m_medianU[0];
    m_ndf[0][0] = 1.0 / (1 + m_nsrate[0][0] * m_dt);
    m_nQ[1][0] = 1.0;
    m_nQd[1][0] = 1.0;
    m_nvol[0] = m_vol[1];
    sqdt = std::sqrt(m_dt);

    for (i = 2; i < m_nnode; ++i) {
        f = 0.0;
        fx = 0.0;
        k = 0;
        x1 = 1.0;
        iter = 0;
        iv = m_vol[i - 1];

        PDB = 1.0 / std::pow((1.0 + m_rate[i - 1]), m_time[i - 1]);

        do {
            f = (x1 + std::pow(x1, std::exp(-2 * iv * sqdt))) / (2.0 * (1.0 + m_nsrate[0][0] * m_dt));
            fx = (1.0 / (2.0 * (1.0 + m_nsrate[0][0] * m_dt)))
                 * (1 + std::exp(-2.0 * iv * sqdt) * (std::pow(x1, std::exp(-2 * iv * sqdt) - 1)));
            x1 = x1 - (f - PDB) / fx;
            Error = std::fabs(f - PDB);
            iter++;
        } while (Error > Epsilon && iter < 1000);
        Pu[i] = x1;
        Pd[i] = std::pow(x1, std::exp(-2.0 * iv * sqdt));
    }

    // with newton_raphson method
    int npram, m, n;
    double **smatrix, *unknwon, *known;

    npram = 2;
    smatrix = new double *[npram];
    for (i = 0; i < npram; ++i) { smatrix[i] = new double[npram]; }

    unknwon = new double[npram];
    known = new double[npram];

    double g, fy, gx, gy;
    int ju, jd;
    double minerror, optx1, opty1;

    for (i = 1; i < m_nnode - 1; ++i) {
        iter = 0;
        x1 = m_rate[i];
        y1 = m_vol[i];

        minerror = 1.0e15;
        optx1 = x1;
        opty1 = y1;

        do {
            for (m = 0; m < npram; ++m) {
                for (n = 0; n < npram; ++n) { smatrix[m][n] = 0; }
                unknwon[m] = 0.0;
                known[m] = 0.0;
            }

            f = g = fx = fy = gx = gy = 0.0;

            // Gauss - newton 방법
            for (j = 0; j < i; ++j) {
                ju = i - 2 * j;

                f += m_nQ[i][j] / (1.0 + x1 * std::exp(y1 * ju * sqdt) * m_dt);

                fx -= m_nQ[i][j] * std::pow((1.0 + x1 * std::exp(y1 * ju * sqdt) * m_dt), -2)

                      * std::exp(y1 * ju * sqdt) * m_dt;
                fy -= m_nQ[i][j] * std::pow((1.0 + x1 * std::exp(y1 * ju * sqdt) * m_dt), -2)
                      * std::exp(y1 * ju * sqdt * m_dt) * std::exp(y1 * ju * sqdt);

                jd = i - 2 * (j + 1);

                g += m_nQ[i][j] / (1.0 + x1 * std::exp(y1 * jd * sqdt) * m_dt);

                gx -= m_nQ[i][j] * std::pow((1.0 + x1 * std::exp(y1 * jd * sqdt) * m_dt), -2)
                      * std::exp(y1 * jd * sqdt) * m_dt;

                gy -= m_nQ[i][j] * std::pow((1.0 + x1 * std::exp(y1 * jd * sqdt) * m_dt), -2)
                      * std::exp(y1 * jd * sqdt * m_dt) * std::exp(y1 * jd * sqdt);
            }

            smatrix[0][0] += (fx * fx + gx * gx);
            smatrix[0][1] += (fx * fy + gx * gy);
            smatrix[1][0] += (fy * fx + gy * gx);
            smatrix[1][1] += (fy * fy + gy * gy);
            known[0] -= (f - Pu[i + 1]) * fx + (g - Pd[i + 1]) * gx;
            known[1] -= (f - Pu[i + 1]) * fy + (g - Pd[i + 1]) * gy;

            gaussian_elimination(smatrix, known, unknwon, npram);

            x1 += unknwon[0];
            y1 += unknwon[1];

            f = g = 0.0;
            for (j = 0; j < i; ++j) {
                ju = i - 2 * j;
                f += m_nQ[i][j] / (1.0 + x1 * std::exp(y1 * ju * sqdt) * m_dt);
                jd = i - 2 * (j + 1);
                g += m_nQ[i][j] / (1.0 + x1 * std::exp(y1 * jd * sqdt) * m_dt);
            }

            Error = std::fabs(f - Pu[i + 1]) + std::fabs(g - Pd[i + 1]);

            if (Error < minerror) {
                minerror = Error;
                optx1 = x1;
                opty1 = y1;
            }
            iter++;
        } while (Error > Epsilon && iter < 1000);

        m_medianU[i] = optx1;
        m_nvol[i] = opty1;

        k = 0;
        if (i < m_nnode) {
            for (j = 0; j <= i; ++j) {
                k = i - 2 * j;
                m_nsrate[i][j] = m_medianU[i] * std::exp(m_vol[i] * k * sqdt);
                m_ndf[i][j] = 1.0 / (1.0 + m_nsrate[i][j] * m_dt);
            }
        }

        if (i < m_nnode - 1) {
            m_nQ[i + 1][0] = m_pu * m_nQ[i][0] * m_ndf[i][0];
            m_nQd[i + 1][0] = m_pu * m_nQd[i][0] * m_ndf[i][1];
            for (j = 1; j < i; ++j) {
                m_nQ[i + 1][j] = m_pu * m_nQ[i][j - 1] * m_ndf[i][j - 1]
                                 + m_pd * m_nQ[i][j] * m_ndf[i][j];

                m_nQd[i + 1][j] = m_pu * m_nQd[i][j - 1] * m_ndf[i][j]
                                  + m_pd * m_nQd[i][j] * m_ndf[i][j + 1];
            }
            m_nQ[i + 1][i] = m_pd * m_nQ[i][i - 1] * m_ndf[i][i - 1];
            m_nQd[i + 1][i] = m_pd * m_nQd[i][i - 1] * m_ndf[i][i];
        }
    }

    delete[] Pu;
    delete[] Pd;
}
