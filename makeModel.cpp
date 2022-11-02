//
// Created by jewoo on 2022-11-02.
//

#include <iostream>

#include <iomanip>
#include "vasicek.h"
#include "cir.h"
#include "makeModel.h"
#include "assetPrice.h"
#include "bdt.h"
#include "cbdt.h"
#include "mathlib.h"

void make_CVasicek() {
    std::cout << std::setprecision(20);

    double alpha, theta, sigma, ini_shor_r;
    alpha = 0.1;
    theta = 0.1;
    sigma = 0.02;
    ini_shor_r = 0.1;

    CVasicek vasicek(alpha, theta, sigma);

    double today, bondmaturity;
    double notional, bondvalue, spotrate;

    today = 0.0;
    bondmaturity = 5.0;
    notional = 100.;

    bondvalue = vasicek.zero_bond_value(ini_shor_r, today, bondmaturity) * notional;
    std::cout << "Zero coupon bond value : " << bondvalue << std::endl;

    spotrate = vasicek.spot_rate(ini_shor_r, today, bondmaturity);
    std::cout << "spot rate : " << spotrate << std::endl;

    double coupon_rate;
    int coupon_frequency;
    double optionmatuity, strike;
    double zerooption, couponbondoption;

    strike = 98.0;
    coupon_rate = 0.1;
    coupon_frequency = 2;
    optionmatuity = 3.0;

    zerooption = vasicek.zero_bond_option_value(-1, notional, ini_shor_r, today,
                                                optionmatuity, bondmaturity, strike);
    std::cout << "Zero coupon bond option value : " << zerooption << std::endl;


    couponbondoption = vasicek.fixed_coupon_bond_option_value(-1, notional,
                                                              coupon_rate, coupon_frequency, ini_shor_r,
                                                              today, optionmatuity,
                                                              bondmaturity, strike);
    std::cout << "fixed coupon bond option value : " << couponbondoption << std::endl;


}

void make_CCirk() {
    std::cout << std::setprecision(20);

    double alpha, theta, sigma, ini_shor_r;
    alpha = 0.1;
    theta = 0.1;
    sigma = 0.02;
    ini_shor_r = 0.1;

    CCIR ccir(alpha, theta, sigma);

    double today, bondmaturity;
    double notional, bondvalue, spotrate;

    today = 0.0;
    bondmaturity = 5.0;
    notional = 100.;

    bondvalue = ccir.zero_bond_value(ini_shor_r, today, bondmaturity) * notional;
    std::cout << "Zero coupon bond value : " << bondvalue << std::endl;

    spotrate = ccir.spot_rate(ini_shor_r, today, bondmaturity);
    std::cout << "spot rate : " << spotrate << std::endl;

    double coupon_rate;
    int coupon_frequency;
    double optionmatuity, strike;
    double zerooption, couponbondoption;

    strike = 98.0;
    coupon_rate = 0.1;
    coupon_frequency = 2;
    optionmatuity = 3.0;

    zerooption = ccir.zero_bond_option_value(-1, notional, ini_shor_r, today,
                                             optionmatuity, bondmaturity, strike);
    std::cout << "Zero coupon bond option value : " << zerooption << std::endl;


    couponbondoption = ccir.fixed_coupon_bond_option_value(-1, notional,
                                                           coupon_rate, coupon_frequency, ini_shor_r,
                                                           today, optionmatuity,
                                                           bondmaturity, strike);
    std::cout << "fixed coupon bond option value : " << couponbondoption << std::endl;
}

void make_assets() {
    std::cout << std::setprecision(20);


    int nrtime;
    double fixed_vol;
    double *rtime, *rate;
    nrtime = 10;
    fixed_vol = 0.15;

    rtime = new double[nrtime];
    rate = new double[nrtime];

    rtime[0] = 1.0;
    rate[0] = 0.0627;
    rtime[1] = 2.0;
    rate[1] = 0.0630;
    rtime[2] = 3.0;
    rate[2] = 0.0641;
    rtime[3] = 4.0;
    rate[3] = 0.0651;
    rtime[4] = 5.0;
    rate[4] = 0.0661;
    rtime[5] = 6.0;
    rate[5] = 0.0672;
    rtime[6] = 7.0;
    rate[6] = 0.0682;
    rtime[7] = 8.0;
    rate[7] = 0.0693;
    rtime[8] = 9.0;
    rate[8] = 0.0704;
    rtime[9] = 10.0;
    rate[9] = 0.0713;

    int coupon_frequency = 1;
    double bond_maturity = 10.0;
    double node_dt = 1.0;
    double notional = 100.0;
    double coupon_rate = 0.0;

    CBDT bdt_vf(bond_maturity, node_dt);
    bdt_vf.build_bdt_tree_fixed_vol(nrtime, rtime, rate, fixed_vol);

    double bond_price;
    AssetCPrice vf_zc_price(notional, coupon_rate, bond_maturity, coupon_frequency);
    bond_price = vf_zc_price.zero_coupon_bond_price(bdt_vf); // 채권 가격 계산

    delete[] rate;
    delete[] rtime;
}

void make_assets2() {
    std::cout << std::setprecision(20);


    int nrtime;
    double fixed_vol;
    double *rtime, *rate;
    nrtime = 10;


    rtime = new double[nrtime];
    rate = new double[nrtime];

    rtime[0] = 1.0;
    rate[0] = 0.0627;
    rtime[1] = 2.0;
    rate[1] = 0.0630;
    rtime[2] = 3.0;
    rate[2] = 0.0641;
    rtime[3] = 4.0;
    rate[3] = 0.0651;
    rtime[4] = 5.0;
    rate[4] = 0.0661;
    rtime[5] = 6.0;
    rate[5] = 0.0672;
    rtime[6] = 7.0;
    rate[6] = 0.0682;
    rtime[7] = 8.0;
    rate[7] = 0.0693;
    rtime[8] = 9.0;
    rate[8] = 0.0704;
    rtime[9] = 10.0;
    rate[9] = 0.0713;

    double *vtime, *vol;
    vtime = new double[nrtime];
    vol = new double[nrtime];

    vtime[0] = 1.0;
    vol[0] = 0.17;
    vtime[1] = 2.0;
    vol[1] = 0.16;
    vtime[2] = 3.0;
    vol[2] = 0.15;
    vtime[3] = 4.0;
    vol[3] = 0.14;
    vtime[4] = 5.0;
    vol[4] = 0.13;
    vtime[5] = 6.0;
    vol[5] = 0.12;
    vtime[6] = 7.0;
    vol[6] = 0.11;
    vtime[7] = 8.0;
    vol[7] = 0.10;
    vtime[8] = 9.0;
    vol[8] = 0.09;
    vtime[9] = 10.0;
    vol[9] = 0.08;

    int coupon_frequency = 1;
    double bond_maturity = 10.0;
    double node_dt = 1.0;
    double notional = 100.0;
    double coupon_rate = 0.0;

    CCBDT bdt_vf(bond_maturity, node_dt);
    bdt_vf.build_bdt_tree_vol_curve(nrtime, rtime, rate, nrtime, vtime, vol);

    double bond_price;
    AssetCPrice vf_zc_price(notional, coupon_rate, bond_maturity, coupon_frequency);
    bond_price = vf_zc_price.zero_coupon_bond_price_with_vol(bdt_vf); // 채권 가격 계산

    delete[] rate;
    delete[] rtime;
}