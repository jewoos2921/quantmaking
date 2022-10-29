//
// Created by jewoo on 2022-10-29.
//

#pragma once

void nelson_siegel_zerorate(int cal_type, int nztime, double *ztime, double *zrate,
                            double *nsrate);

void get_ns_parameter_with_zerorate(int nztime, double *ztime, double *zrate,
                                    int nguess, double *guess);

void get_ns_parameter_with_bondprice(int nztime, double *ztime, double *zrate,
                                     int nguess, double *guess);

double price_ns_discountfactor(double *guess, double time);

double price_ns_discountfactor(double *guess, double time, int j, double h);

double ns_zerorate_function(double *parameter, double time);

double ns_zerorate_function(double *parameter, double time, int j, double h);
