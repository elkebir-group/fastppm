//
// Created by Yuanyuan Qi on 8/1/24.
//

#include "Dual.h"
#include "Primal.h"

void Dual::dual_from_primal(const Primal &primal) {
    k = primal.k + 1;
    for (int i = 0; i < k - 1; i++) {
        x[i] = primal.slope.at(i);
        y[i] = primal.y.at(i) - primal.slope.at(i) * primal.x.at(i);
    }
    for (int i = 0; i < k; i++) {
        slope[i] = -primal.x.at(i);
    }
}

#ifdef _DEBUG

bool Dual::self_check() const {
    for (int i = 0; i < k - 2; i++) {
        if (x[i+1] <= x[i]){
            return false;
        }
        if (abs((y[i+1]-y[i])/(x[i+1]-x[i]) - slope[i+1]) > Compare_eps){
            return false;
        }
    }
    return true;
}

#endif