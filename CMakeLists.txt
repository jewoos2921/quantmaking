cmake_minimum_required(VERSION 3.24)
project(quantmaking)

set(CMAKE_CXX_STANDARD 20)

add_executable(quantmaking main.cpp index.cpp index.h yield.cpp yield.h option.cpp option.h)
