//
// Created by 17491433 on 11.04.2020.
//

#include "Distribution.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>

LCG *Distribution::getLCG() {
    return lcg;
}

Distribution::Distribution() {

}

int64_t Distribution::getRandomBernuliValue(double p) {
    if(p >= lcg->getPsevdoRandom()) {
        return 1;
    } else {
        return 0;
    }
}

int64_t Distribution::getRandomBenominalValue(int64_t n, double p) {
    int64_t result(0);
    double value = lcg->getPsevdoRandom();
    while(value < p && value != 0 && result < n) {
        value = getLCG()->getPsevdoRandom();
        result++;
    }
    return result;
}

int64_t Distribution::getRandomGeometricalValue(double p) {
    return std::floor(getRandomExponentialValue(-std::log(1-p)));
}

int64_t Distribution::getRandomHGeometricalValue(int64_t n, int64_t N, int64_t K) {
    int64_t result = 0, p = static_cast<double>(K) / N;
    for (int i = 1; i <= n; ++i)
    {
        if (getRandomBernuliValue(p) && ++result == K)
            return result;
        getRandomBernuliValue(static_cast<double>(K - result) / static_cast<double>(N - i));
    }
    return result;
}

double Distribution::getRandomLinearValue(double a, double b) {
    return a + (b - a)*lcg->getPsevdoRandom();
}

double Distribution::getRandomNormalValue(double M, double sigma) {
    double sum(0);
    for (int i = 0; i < 12; ++i) {
        sum+=lcg->getPsevdoRandom();
    }
    return M + (sum - 6)* sigma;
}

double Distribution::getRandomPuasonValue(double rate) {
    double expRate = exp(-1*rate);
    double k = 0, prod = lcg->getPsevdoRandom();
    while (prod > expRate) {
        prod *= lcg->getPsevdoRandom();
        ++k;
    }
    return k;
}

double Distribution::getRandomExponentialValue(double rate) {
    return (-1.0/rate)*log(lcg -> getPsevdoRandom());
}
