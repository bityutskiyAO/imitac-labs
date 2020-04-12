cmake_minimum_required(VERSION 3.15)
project(laba_2)

set(CMAKE_CXX_STANDARD 14)

add_executable(laba_2 main.cpp lcg/LCG.cpp lcg/LCG.h distributions/Distribution.h distributions/Distribution.cpp distributions/DistributionsTest.cpp distributions/DistributionsTest.h)