cmake_minimum_required(VERSION 3.24)
project(quantmaking)

set(CMAKE_CXX_STANDARD 20)

add_executable(quantmaking
        main.cpp
        index.cpp
        yield.cpp
        option.cpp
        mathlib.cpp
        nonlinearsolver.cpp
        data.cpp
        mainNet.cpp
        product.cpp
        price.cpp
        makeOption.cpp
        sobol.cpp
        priceAlpha.cpp
        nelson_siegel.cpp
        vasicek.cpp
        makeModel.cpp
        cir.cpp
        bdt.cpp
        assetPrice.cpp
        cbdt.cpp
        hwtree.cpp
        chwtPrice.cpp
        priceBlack.cpp


        quantum/lib/xgates.cpp
        quantum/benchmarking/benchmark.cpp
        quantum/libq/apply.cpp
        quantum/libq/gates.cpp
        quantum/libq/qureg.cpp quantum/libq/gates_jit.cpp)
