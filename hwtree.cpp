//
// Created by jewoo on 2022-11-09.
//
#include <iostream>
#include <cmath>
#include "hwtree.h"
#include "mathlib.h"

CHWTree::CHWTree() {
    m_nnode = 10;
    m_jmax = 5;
    m_alpha = 0.1;
    m_sigma = 0.01;
    m_dt = 1.0;
}

CHWTree::~CHWTree() {
    int i;

    // probability
    delete[] m_pu;
    delete[] m_pm;
    delete[] m_pd;

    for (i = 0; i < m_nnode; ++i) {
        delete[] m_ndf[i];
        delete[] m_nrate[i];
    }
}


void CHWTree::build_hull_white_short_rate_tree(int nrtime, double *time, double *rate) {
    std::cout << "Hull-White Model" << std::endl;

    /*
     원래는
     2        ******************************************
     1       *******************************************
     0      ********************************************
    -1       *******************************************
    -2        ******************************************
     이런 모습이여야 하나, 배열이 0,0에서 시작됨에 따라
     0      ********************************************
     1       *******************************************
     2       *******************************************
     3        ******************************************
     4        ******************************************
    으로 변화시켜 트리 생성했으며 j의 값은 원점을 기준으로 shift하였음
     * */

    /* 각노드에서 Tree 전개 모양
      1번 확장 모양: 일반 경우
         *
     *   *
         *
      2번 확장 모양: Jmax에 도달시
     *   *
         *
         *
      3번 확장 모양: Jmin에 도달시
         *
         *
     *   *
     와 같이 확장
     */

    int i, j, k, m, imax;
    double dR, sum;
    double *P, *itime, *irate;
    double **nQ;
    double alpha_rate;

    m_jmax = static_cast<int>((ceil(0.1835 / (m_alpha * m_dt)))); // Jmax 계산

    if (m_nnode - 1 <= m_jmax) {
        imax = 2 * (m_nnode - 1) + 1;
    } else {
        imax = 2 * m_jmax + 1;
    }

    P = new double[m_nnode]; // 할인채
    itime = new double[m_nnode];
    irate = new double[m_nnode];
    nQ = new double *[m_nnode];

    m_ndf = new double *[m_nnode];
    m_nrate = new double *[m_nnode];

    for (i = 0; i < m_nnode; ++i) {
        nQ[i] = new double[imax];
        m_ndf[i] = new double[imax];
        m_nrate[i] = new double[imax];
    }


    m_pu = new double[imax];
    m_pm = new double[imax];
    m_pd = new double[imax];

    for (i = 0; i < m_nnode; ++i) {
        for (j = 0; j < imax; ++j) { nQ[i][j] = 0.0; }
        P[i] = 0.0;
        itime[i] = (i + 1) * m_dt;
        irate[i] = linear_interpolation(nrtime, time, rate, itime[i]);
    }

    dR = m_sigma * sqrt(3.0 * m_dt);
    m_nrate[0][0] = 0.0;

    for (i = 1; i < m_nnode; ++i) {
        if (i <= m_jmax) {
            imax = 2 * i;
        } else {
            imax = 2 * m_jmax;
        }
        for (j = 0; j <= imax; ++j) {
            if (i <= m_jmax) {
                k = i - j;
            } else {
                k = m_jmax - j;
            }
            m_nrate[i][j] = k * dR;
        }
    }

    // 노드에서의 확률 계산
    if (m_nnode - 1 <= m_jmax) {
        // node 수가 jmax보다 작으므로 3방향으로 퍼져나감, 1번 확장 모양
        for (j = 0; j <= 2 * (m_nnode - 1); ++j) {
            k = m_nnode - 1 - j;
            m_pu[j] = 1.0 / 6.0 + (pow((m_alpha * k * m_dt), 2.0) - m_alpha * k * m_dt) / 2.0;
            m_pm[j] = 2.0 / 3.0 - pow((m_alpha * k * m_dt), 2.0);
            m_pd[j] = 1.0 / 6.0 + (pow((m_alpha * k * m_dt), 2.0) + m_alpha * k * m_dt) / 2.0;
        }
    } else {
        // node 수가 jmax보다 클때이므로 Jmax 또는 Jmin 시 수렴함
        for (j = 0; j <= 2 * m_jmax; ++j) {
            k = m_jmax - j;

            if (j == 0) { // Jmax 시점, 2번 확장 모형
                m_pu[j] = 7.0 / 6.0 + (pow((m_alpha * k * m_dt), 2.0)
                                       - 3.0 * m_alpha * k * m_dt) / 2.0;

                m_pm[j] = -1.0 / 3.0 - pow((m_alpha * k * m_dt), 2.0)
                          + 2.0 * m_alpha * k * m_dt;

                m_pd[j] = 1.0 / 6.0 + (pow((m_alpha * k * m_dt), 2.0)
                                       - m_alpha * k * m_dt) / 2.0;

            } else if (j == 2 * m_jmax) { // Jmin 시점, 3번 확장 모형
                m_pu[j] = 1.0 / 6.0 + (pow((m_alpha * k * m_dt), 2.0)
                                       + m_alpha * k * m_dt) / 2.0;

                m_pm[j] = -1.0 / 3.0 - pow((m_alpha * k * m_dt), 2.0)
                          - 2.0 * m_alpha * k * m_dt;

                m_pd[j] = 7.0 / 6.0 + (pow((m_alpha * k * m_dt), 2.0)
                                       + 3.0 * m_alpha * k * m_dt) / 2.0;

            } else { // 기타일반 시점, 1번 확장 모형
                m_pu[j] = 1.0 / 6.0 + (pow((m_alpha * k * m_dt), 2.0)
                                       - m_alpha * k * m_dt) / 2.0;

                m_pm[j] = 2.0 / 3.0 - pow((m_alpha * k * m_dt), 2.0);

                m_pd[j] = 1.0 / 6.0 + (pow((m_alpha * k * m_dt), 2.0)
                                       + m_alpha * k * m_dt) / 2.0;

            }
        }
    }

    for (i = 0; i < m_nnode; ++i) {
        P[i] = exp(-irate[i] * itime[i]);
    }

    nQ[0][0] = 1.0;
    alpha_rate = irate[0];
    m_nrate[0][0] += alpha_rate;
    m_ndf[0][0] = exp(-m_nrate[0][0] * m_dt);

    // 노드에서의 상태가격 계산
    for (i = 1; i < m_nnode; ++i) {
        // node 수가 jamx보다 작은 경우, 1번 확장 모양
        if (i <= m_jmax) {
            imax = 2 * i - 2; // j, j-1, j+2까지 고려해주므로 -2로 적용
            for (j = 0; j <= imax; ++j) {
                if (m_nnode - 1 <= m_jmax) {
                    m = m_nnode - i + j;
                } else {
                    m = m_jmax - i + 1 + j;
                }

                // (i, j)노드, (i, j+1)노드 , (i, j+2)노드에서 상태 가격
                nQ[i][j] += nQ[i - 1][j] * m_pu[m] * exp(-m_nrate[i - 1][j] * m_dt);
                nQ[i][j + 1] += nQ[i - 1][j] * m_pm[m] * exp(-m_nrate[i - 1][j] * m_dt);
                nQ[i][j + 2] += nQ[i - 1][j] * m_pd[m] * exp(-m_nrate[i - 1][j] * m_dt);
            }
        } else {
            imax = 2 * m_jmax;

            for (j = 0; j <= imax; ++j) {
                if (j == 0) {
                    k = 0;
                } else if (j == imax) {
                    k = -2;
                } else {
                    k = -1;
                }
                m = j;
                // (i, j)노드, (i, j+1)노드 , (i, j+2)노드에서 상태 가격
                nQ[i][j + k] += nQ[i - 1][j] * m_pu[m] * exp(-m_nrate[i - 1][j] * m_dt);
                nQ[i][j + k + 1] += nQ[i - 1][j] * m_pm[m] * exp(-m_nrate[i - 1][j] * m_dt);
                nQ[i][j + k + 2] += nQ[i - 1][j] * m_pd[m] * exp(-m_nrate[i - 1][j] * m_dt);
            }
        }

        if (i <= m_jmax) {
            imax = 2 * i;
        } else {
            imax = 2 * m_jmax;
        }
        sum = 0.0;

        for (j = 0; j <= imax; ++j) { sum += nQ[i][j] * exp(-m_nrate[i][j] * m_dt); }

        alpha_rate = log(sum / P[i]) / m_dt; // 보정값 계산

        for (j = 0; j <= imax; ++j) {
            m_nrate[i][j] += alpha_rate; // 각노드에 보정값 적용
            if (m_nrate[i][j] < 0) {
                m_nrate[i][j] = 0.0;
            }
            m_ndf[i][j] = exp(-m_nrate[i][j] * m_dt);
        }
    }

    for (i = 0; i < m_nnode; ++i) {
        delete[] nQ[i];
    }

    delete[] nQ;
    delete[] P;
    delete[] itime;
    delete[] irate;
}
