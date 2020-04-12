//
// Created by 17491433 on 11.04.2020.
//

#ifndef LABA_2_DISTRIBUTION_H
#define LABA_2_DISTRIBUTION_H

#include "../lcg/LCG.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>

class Distribution {
private:
    LCG* lcg = new LCG();
public:
    Distribution();
    int64_t getRandomBernuliValue(double p);
    int64_t getRandomBenominalValue(int64_t n, double p);
    int64_t getRandomGeometricalValue(double p);
    int64_t getRandomHGeometricalValue(int64_t n, int64_t N, int64_t K);
    double getRandomLinearValue(double a, double b);
    double getRandomNormalValue(double M, double sigma);
    double getRandomPuasonValue(double rate);
    double getRandomExponentialValue(double rate);
    LCG* getLCG();
};


#endif //LABA_2_DISTRIBUTION_H
