//
// Created by m1707760 on 17.03.2020.
//

#ifndef UNTITLED_LCG_H
#define UNTITLED_LCG_H

#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>

using namespace std;

class LCG {
private:
    int64_t m, a, b, currentXi, currentXni, capacity, Xo, n = pow(2, 20);
    int64_t d = 10;
public:
    LCG();
    double getPsevdoRandom();
    double getTransitionPsevdoRandom();
    int64_t modPow(int64_t a, int64_t n, int64_t m);
    void transitionProcedure();
    int64_t getConstant();
    int64_t getSize();
};


#endif //UNTITLED_LCG_H
