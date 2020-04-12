//
// Created by 17491433 on 12.04.2020.
//

#include "DistributionsTest.h"
#include <cmath>

double DistributionsTest::xiSquaredTest(const std::vector<int64_t> &v, const std::vector<double> &p, int64_t n) {
    double result(0.0);
    for (int64_t i(0); i < v.size(); i += 1) {
        result += pow(v[i] - n * p[i], 2) / (n * p[i]);
    }
    return result;
}
double factorial(int64_t n) {
    if(n > 1) {
        return n*factorial(n-1);
    }
    return 1;
}

double sochet(int64_t k, int64_t n) {
    double result = factorial(n)/(factorial(k) * factorial(n-k));
    return result;
}

double DistributionsTest::freqBernuli(double p) {
    distribution = new Distribution();
    int64_t n = distribution->getLCG()->getSize();
    vector<int64_t> frequencyVector (2, 0);
    for (int64_t i = 0; i < n; ++i) {
        frequencyVector[distribution->getRandomBernuliValue(p)] ++;
    }
    std::vector<double> P(2, 0);
    P[0] = 1-p;
    P[1] = p;
    return xiSquaredTest(frequencyVector, P, n);
}

double DistributionsTest::freqBenom(int64_t n, double p) {
    distribution = new Distribution();
    int64_t size = distribution->getLCG()->getSize();
    vector<int64_t> frequencyVector(n + 1 ,0);
    for (int64_t i = 0; i < size; ++i) {
        frequencyVector[distribution->getRandomBenominalValue(n , p)] ++;
    }
    std::vector<double> P(n + 1, 0);
    for (int i = 0; i < P.size() ; ++i) {
        P[i] = sochet(i, n)*pow(p, i)*pow(1-p, n-i);
    }
    return xiSquaredTest(frequencyVector, P, n);
}

double DistributionsTest::freqGeom(double p, int64_t k) {
    distribution = new Distribution();
    int64_t size = distribution->getLCG()->getSize();
    vector<int64_t> frequencyVector(k+1, 0);
    for (int64_t i = 0; i < size; ++i) {
        int64_t rez = distribution->getRandomGeometricalValue(p);
        if(rez > k) {
            frequencyVector[k] = rez;
        } else {
            frequencyVector[rez] ++;
        }
    }
    std::vector<double> P(k + 1, 0);
    double last(0);
    for (int i = 0; i < P.size() - 1 ; ++i) {
        P[i] = p * pow(1-p, i);
        last+=P[i];
    }
    P[k] = 1 - last;
    return xiSquaredTest(frequencyVector, P, size);
}

double DistributionsTest::freqHGeom(int64_t n, int64_t N, int64_t K) {
    distribution = new Distribution();
    int64_t size = distribution->getLCG()->getSize();
    vector<int64_t> frequencyVector(N, 0);
    for (int64_t i = 0; i < size; ++i) {
        frequencyVector[distribution->getRandomHGeometricalValue(n, N, K)] ++;
    }
    std::vector<double> P(N, 0);
    for (int i = 0; i < P.size() ; ++i) {
        P[i] = sochet(i, K) * sochet(n - i, N - K) / sochet(n , N);
    }
    return xiSquaredTest(frequencyVector, P, size);
}

double DistributionsTest::freqPuason(double rate, int64_t k) {
    distribution = new Distribution();
    int64_t size = distribution->getLCG()->getSize();
    vector<int64_t> frequencyVector(k + 1, 0);
    for (int64_t i = 0; i < size; ++i) {
        int64_t rez = distribution->getRandomPuasonValue(rate);
        if(rez > k) {
            frequencyVector[k] = rez;
        } else {
            frequencyVector[rez] ++;
        }
    }
    std::vector<double> P(k + 1, 0);
    double last(0);
    for (int i = 0; i < P.size() - 1; ++i) {
        P[i] = pow(rate, i) * exp(-rate)/factorial(i);
        last+=P[i];
    }
    P[k] = 1 - last;
    return xiSquaredTest(frequencyVector, P, size);
}

double DistributionsTest::freqLinear(double a, double b) {
    return 0;
}

double DistributionsTest::freqNormal(int64_t size, double M, double sigma, double space, double spaceNumber) {
    std::vector<int64_t> frequencyVector(spaceNumber);
    double start = M - space * (spaceNumber / 2);
    for (size_t i = 0; i < size; i++) {
        double value = distribution->getRandomNormalValue(M, sigma);
        size_t ind = spaceNumber - 1;
        double cur_int = start;
        for (size_t i = 0; i < frequencyVector.size() - 1; i++)
        {
            if (value < space + cur_int && value > cur_int)
            {
                ind = i;
                break;
            }
            cur_int += space;
        }
        frequencyVector[ind]++;
    }

    double sample_mean = 0;
    double sample_std = 0;
    double cur_int = start + space / 2;
    for (size_t i = 0; i < frequencyVector.size(); i++) {
        sample_mean += frequencyVector[i]* cur_int;
        sample_std += frequencyVector[i] * cur_int * cur_int;
        cur_int += space;
    }
    sample_std = std::sqrt(sample_std/(size-1));
    sample_mean /= size;

    auto npj = std::vector<double>(spaceNumber);
    double m = space * size / sample_std;
    cur_int = start + space / 2;

    for (size_t i = 0; i < spaceNumber - 1; i++)
    {
        double z = (cur_int - sample_mean)/sample_std;
        npj[i] = m * (std::exp(-z*z/2) / std::sqrt(2*M_PI) );
        npj[spaceNumber - 1] += npj[i];
        cur_int += spaceNumber;
    }
    npj[spaceNumber - 1] = size - npj[spaceNumber - 1];


    return xiSquaredTest(frequencyVector, npj, size);
}

double DistributionsTest::freqExp(int64_t size, double rate, double space, double spaceNumber) {
    double value = distribution->getRandomExponentialValue(rate);
    std::vector<int64_t> frequencyVector(spaceNumber);
    size_t ind = spaceNumber - 1;
    for (size_t i = 0; i < size; i++) {
        for (size_t j = 0; j < frequencyVector.size() - 1; j++) {
            if (value < space * (j + 1) && value > space * j) {
                ind = j;
                break;
            }
        }
        frequencyVector[ind]++;
    }
    std::vector<double> P(spaceNumber);
    double last(0);
    for (size_t i = 0; i < spaceNumber - 1; i++) {
        P[i] = size * (1 - std::exp(-rate * space * (i + 1)));
        last+=P[i];
    }
    P[spaceNumber - 1] = 1 - last;
    return xiSquaredTest(frequencyVector, P, size);
}

DistributionsTest::DistributionsTest() {

}
