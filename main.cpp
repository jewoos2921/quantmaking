#include <iostream>
#include <cmath>
#include <iomanip>
#include "index.h"
#include "yield.h"
#include "option.h"
#include "mathlib.h"
#include "nonlinearsolver.h"

void make_gaussian_elimination();

void make_tridiagonal_elimination();

void choleski();

void make_bisection_method_and_newton();

void make_gauss_newton();

int main() {
//    make_gaussian_elimination();
//    printf("\n");
//    make_tridiagonal_elimination();
//    printf("\n");
//    choleski();
    make_bisection_method_and_newton();
    printf("\n");
    make_gauss_newton();
}


void make_european_calloption() {
    CIndex index(110., 0.25, 0.03); // 기초 자산, 변동성, 배당률
    CYield yield(0.055); // 무위험이자율
    COption option(105.0, 1.0); // 행사가격, 잔존 만기


    option.european_calloption_price(index, yield);

    std::cout << "European Call option 가격: " << option.m_optionprice << std::endl;
}

void make_gaussian_elimination() {
    // 4원 1차 연립 방정식 해 찾기
    int i, neq;
    double **smatrix;
    double *known_value;
    double *unknown_value;

    neq = 4;
    smatrix = new double *[neq];
    for (i = 0; i < neq; ++i) {
        smatrix[i] = new double[neq];
    }
    known_value = new double[neq];
    unknown_value = new double[neq];

    smatrix[0][0] = 1;
    smatrix[0][1] = 1;
    smatrix[0][2] = 0;
    smatrix[0][3] = 3;

    smatrix[1][0] = 2;
    smatrix[1][1] = 1;
    smatrix[1][2] = -1;
    smatrix[1][3] = 1;

    smatrix[2][0] = 3;
    smatrix[2][1] = -1;
    smatrix[2][2] = -1;
    smatrix[2][3] = 2;

    smatrix[3][0] = -1;
    smatrix[3][1] = 2;
    smatrix[3][2] = 3;
    smatrix[3][3] = 1;

    known_value[0] = 4;
    known_value[1] = 1;
    known_value[2] = -3;
    known_value[3] = 4;

    gaussian_elimination(smatrix, known_value, unknown_value, neq);

    for (i = 0; i < neq; ++i) {
        std::cout << "x_" << i + 1 << " : " << unknown_value[i] << std::endl;
    }

    for (i = 0; i < neq; ++i) {
        delete smatrix[i];
    }
    delete[] smatrix, known_value, unknown_value;
}

void make_tridiagonal_elimination() {
    // 4원 1차 연립 방정식 해 찾기
    int i, neq;
    double **smatrix;
    double *known_value;
    double *unknown_value;

    neq = 4;

    smatrix = new double *[neq];
    for (i = 0; i < neq; ++i) {
        smatrix[i] = new double[neq];
    }
    known_value = new double[neq];
    unknown_value = new double[neq];

    smatrix[0][0] = 2;
    smatrix[0][1] = -1;
    smatrix[0][2] = 0;
    smatrix[0][3] = 0;

    smatrix[1][0] = -1;
    smatrix[1][1] = 2;
    smatrix[1][2] = -1;
    smatrix[1][3] = 0;

    smatrix[2][0] = 0;
    smatrix[2][1] = -1;
    smatrix[2][2] = 2;
    smatrix[2][3] = -1;

    smatrix[3][0] = 0;
    smatrix[3][1] = 0;
    smatrix[3][2] = -1;
    smatrix[3][3] = 2;

    known_value[0] = 1;
    known_value[1] = 0;
    known_value[2] = 0;
    known_value[3] = 1;

    tridiagonal_elimination(smatrix, known_value, unknown_value, neq);

    for (i = 0; i < neq; ++i) {
        std::cout << "x_" << i + 1 << " : " << unknown_value[i] << std::endl;
    }

    for (i = 0; i < neq; ++i) {
        delete smatrix[i];
    }
    delete[] smatrix, known_value, unknown_value;
}

void choleski() {

    int ncol;
    double **A; // 행렬 선언
    double *UV; // 미지수, A x = b 에서 x
    double *KV; // 기지수, A x = b 에서 b

    ncol = 3;
    A = new double *[ncol];
    for (int i = 0; i < ncol; ++i) {
        new double[ncol];
    }
    UV = new double[ncol];
    KV = new double[ncol];

    A[0][0] = 4.0;
    A[0][1] = -4.0;
    A[0][2] = 8.0;

    A[1][0] = -4.0;
    A[1][1] = 5.0;
    A[1][2] = -11.0;

    A[2][0] = 1.0;
    A[2][1] = -11.0;
    A[2][2] = 32.0;

    KV[0] = 1.0;
    KV[1] = 1.5;
    KV[2] = 2;

    cholesky_solver(A, UV, KV, ncol);
    gaussian_elimination(A, KV, UV, ncol);
    std::cout << "가우스 소거법을 이용한 연립방정식의 해는 " << std::endl;
    for (int i = 0; i < ncol; ++i) {
        std::cout << "x_" << i + 1 << " : " << UV[i] << std::endl;
    }
    for (int i = 0; i < ncol; ++i) {
        delete[]A[i];
    }

    delete[] A;
    delete[] KV;
    delete[] UV;
}

void make_bisection_method_and_newton() {
    double solution;
    std::cout << std::setprecision(30); // 결과 값의 소수점 단위를 확대 표시
    solution = bisection_method();
    std::cout << "bisection method 계산된 해 : " << solution << std::endl;
    solution = newtonraphson_method();
    std::cout << "newtonraphson method 계산된 해 : " << solution << std::endl;
}

void make_gauss_newton() {
    int i, n_data;
    double *mx, *my;

    int n_parameter;
    double *parameter;

    std::cout << std::setprecision(5);
    n_data = 10;
    mx = new double[n_data];
    my = new double[n_data];

    for (i = 0; i < n_data; ++i) {
        mx[i] = 1.0 + 0.2 * i;
        my[i] = 1.561 * mx[i] - 0.725 * mx[i] * mx[i] - 0.532 * std::exp(-0.288 * mx[i]);
    }

    n_parameter = 4;
    parameter = new double[n_parameter];
    for (i = 0; i < n_parameter; ++i) {
        parameter[i] = 1.0; // 파라미터 초기화
    }

    gauss_newton_parameter_solver(n_parameter, mx, my, n_parameter, parameter, fx_function);
    // 비선형 파라미터 추정
    for (i = 0; i < n_parameter; ++i) {
        std::cout << "p[" << i << "] : " << parameter[i] << std::endl;
    }

    delete[] parameter;
    delete[] mx;
    delete[] my;
}

