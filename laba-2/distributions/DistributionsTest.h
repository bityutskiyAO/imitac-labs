//
// Created by 17491433 on 12.04.2020.
//

#ifndef LABA_2_DISTRIBUTIONSTEST_H
#define LABA_2_DISTRIBUTIONSTEST_H

#include "../lcg/LCG.h"
#include "Distribution.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>

class DistributionsTest {
private:
    Distribution* distribution;
public:
    DistributionsTest();
    double xiSquaredTest(const std::vector<int64_t>& v, const std::vector<double>& p, int64_t n);
    double freqBernuli(double p);
    double freqBenom(int64_t n, double p);
    double freqGeom(double p, int64_t k);
    double freqHGeom(int64_t n, int64_t N, int64_t K);
    double freqPuason(double rate, int64_t k);
    double freqLinear(double a, double b);
    double freqNormal(int64_t size, double M, double sigma, double space, double spaceNumber);
    double freqExp(int64_t size, double rate, double space, double spaceNumber);
};


#endif //LABA_2_DISTRIBUTIONSTEST_H
