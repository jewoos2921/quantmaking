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
        )
