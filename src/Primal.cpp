//
// Created by Yuanyuan Qi on 8/1/24.
//

#include "Primal.h"
#include <algorithm>

void Primal::update(real begin, real end, int _k, func_base * func) {
    k = _k;
    for (int i = 0; i <= k; i++) {
        x[i] = begin + (end - begin) * i / k;
        y[i] = func->operator()(x[i]);
    }
    for (int i = 0; i < k; i++) {
        slope[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
    }
}

void Primal::update(const std::vector<real> &xx, const std::vector<real> &yy) {
    k = xx.size();
    for (int i = 0; i <= k; i++) {
        x[i] = xx[i];
        y[i] = yy[i];
    }
    for (int i = 0; i < k; i++) {
        slope[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
    }
}

real Primal::backtrace_delta(real d){
    int idx = binary_search_less(slope, k - 1, d);
    if (idx>0 && (d-slope[idx-1]) < Compare_eps) idx--;
    return x[idx];
}

#ifdef _DEBUG

bool Primal::self_check() const {
    for (int i = 0; i < k; i++) {
        if (x[i+1] < x[i]){
            return false;
        }
        if (abs((y[i+1]-y[i])/(x[i+1]-x[i]) -slope[i]) > Compare_eps){
            return false;
        }
    }
    return true;
}

#endif
