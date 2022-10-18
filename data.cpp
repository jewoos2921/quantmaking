//
// Created by jewoo on 2022-10-18.
//
#include <iostream>
#include <fstream>
#include "data.h"

void
read_market_vol(char *filename, int n_strike, double *v_strike, int n_maturity, double *v_maturity, double **m_vol) {
    std::ifstream inf(filename, std::ios::in);
    if (!inf.is_open()) {
        std::cout << filename << "파일을 열 수 없 습니다." << std::endl;
        inf.clear();
        exit(0);
    }
    char text[20];
    inf >> text; // "만기:" 읽기

    for (int j = 0; j < n_maturity; ++j) {
        inf >> text;
        v_maturity[j] = std::atof(text);
    }

    inf >> text; // "행사가격:" 읽기
    for (int i = 0; i < n_strike; ++i) {
        inf >> text;
        v_strike[i] = std::atof(text);
    }

    inf >> text; // "변동성:" 읽기
    for (int i = 0; i < n_strike; ++i) {
        for (int j = 0; j < n_maturity; ++j) {
            inf >> text;
            m_vol[i][j] = std::atof(text) / 100.0;
        }
    }
    inf.close();
}

void save_vol(char *filename, int n_strike, double *v_strike, int n_maturity, double *v_maturity, double **vol) {
    std::ofstream outf(filename, std::ios::out | std::ios::trunc);

    // ios::trunc 사용시 기존 데이터 삭제후 추가  
    // ios::app 사용시 데이터 추가 가능
    if (!outf.is_open()) {
        std::cout << filename << "파일을 열 수 없 습니다." << std::endl;
        outf.clear();
        exit(0);
    }
    outf << "만기:" << " ";
    for (int j = 0; j < n_maturity; ++j) {
        outf << v_maturity[j] << " ";
    }
    outf << std::endl;

    outf << "행사가격" << " ";
    for (int i = 0; i < n_strike; ++i) {
        outf << v_strike[i] << " ";
    }
    outf << std::endl;
    outf << "변동성: " << " ";
    for (int i = 0; i < n_strike; ++i) {
        for (int j = 0; j < n_maturity; ++j) {
            outf << vol[i][j] * 100.0 << " ";
        }
        outf << std::endl;
        outf << " " << "";
    }
    outf.close();
}
