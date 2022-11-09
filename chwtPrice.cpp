//
// Created by jewoo on 2022-11-09.
//

#include <iostream>
#include <cmath>
#include "hwtree.h"
#include "chwtPrice.h"
#include "mathlib.h"

CHWTPrice::CHWTPrice() {
    m_coupon_frequency = 1;
    m_ncoupon = 1;
    m_notional = 100.0;
    m_coupon_rate = 0.05;
    m_bond_maturity = 5.0;
}

CHWTPrice::CHWTPrice(double notional, double coupon_rate, double bond_maturity, int coupon_frequency) {
    m_coupon_frequency = coupon_frequency;
    m_notional = notional;
    m_coupon_rate = coupon_rate;
    m_bond_maturity = bond_maturity;
    m_ncoupon = static_cast<int>(m_bond_maturity * m_coupon_frequency);
}

CHWTPrice::~CHWTPrice() = default;

double CHWTPrice::zero_coupon_bond_price(CHWTree &hw) {
    int i, j, k, m, nnode, imax, jmax;
    double price;
    double *oldprice, *newbprice;

    nnode = static_cast<int>(m_bond_maturity / hw.m_dt);

    jmax = hw.m_jmax;

    if (nnode <= jmax) {
        imax = 2 * (nnode - 1) + 1;
    } else {
        imax = 2 * jmax + 1;
    }

    if (nnode > hw.m_nnode) {
        std::cout << "HULL and WHITE 트리를 보다 조밀하게 생성하세요" << std::endl;
        exit(0);
    }

    oldprice = new double[imax];
    newbprice = new double[imax];

    for (j = 0; j < imax; ++j) {
        oldprice[j] = m_notional; // 채권만기시 원금지급
    }

    for (i = nnode - 1; i >= 0; i--) {
        if (j <= jmax) { // 현재 node 가 Jmax보다 작을 경우, 즉 수렴하지 않은 경우
            imax = 2 * i;
        } else {
            // 현재 node 가 Jmax보다 클 경우, 즉 수렴하지 할 경우
            imax = 2 * jmax;
        }

        for (j = 0; j <= imax; ++j) {
            if (i <= jmax) {
                k = i - j; // 현재 node 가 Jmax보다 작을 경우
            } else {
                k = jmax - j; // 현재 node 가 Jmax보다 클 경우
            }

            // 노드의 확률 적용시 node의 수가 Jmax보다 작을 경우
            // 위에서 구한 k와 관계 없음
            if (nnode - 1 <= jmax) {
                m = nnode - 1 - i + j;
            } else {
                if (i <= jmax) {
                    // 노드의 확률 적용시 현재 node가 jMax보다 작은 경우
                    // 위에서 구한 k적용
                    m = jmax - k;
                } else {
                    // 노드의 확률 적용시 현재 node가 jMax보다 큰 경우, m=j적용
                    m = j;
                }
            }

            // 이전 노드 (i)의 3상태(j,j+1,j+2)의 채권가격을 노드의 확률을 고려하여 할인
            if (k == jmax) {
                // 2번 확장 모양 경우
                newbprice[j] =
                        (hw.m_pu[m] * oldprice[j] +
                         hw.m_pm[m] * oldprice[j + 1] +
                         hw.m_pd[m] * oldprice[j + 2]) * hw.m_ndf[i][j];
            } else if (k == -jmax) {
                // 3번 확장 모양 경우
                newbprice[j] = (hw.m_pu[m] * oldprice[j - 2]
                                + hw.m_pm[m] * oldprice[j - 1]
                                + hw.m_pd[m] * oldprice[j]) * hw.m_ndf[i][j];
            } else {
                if (i >= jmax) {
                    // 현재 node가 jMax보다 큰 경우 1번 확장인 경우
                    newbprice[j] = (hw.m_pu[m] * oldprice[j - 1]
                                    + hw.m_pm[m] * oldprice[j]
                                    + hw.m_pd[m] * oldprice[j + 1]) * hw.m_ndf[i][j];
                } else {
                    // 현재 node가 jMax보다 작은 경우 1번 확장인 경우
                    newbprice[j] = (hw.m_pu[m] * oldprice[j]
                                    + hw.m_pm[m] * oldprice[j + 1]
                                    + hw.m_pd[m] * oldprice[j + 2]) * hw.m_ndf[i][j];
                }
            }
        }

        for (j = 0; j <= imax; ++j) { oldprice[j] = newbprice[j]; }
    }

    price = newbprice[0];

    delete[] oldprice;
    delete[] newbprice;

    std::cout << "zero coupon bond price : " << price << std::endl;
    return price;

}
