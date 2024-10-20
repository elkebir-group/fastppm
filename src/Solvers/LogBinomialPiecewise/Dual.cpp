//
// Created by Yuanyuan Qi on 8/1/24.
//

#include "Dual.h"
#include "Primal.h"

namespace LogBinomialPiecewiseLinearSolver {
void Dual::dual_from_primal(const Primal &primal) {
    k = primal.k + 1;
    if (k > slope.size()){
        slope.resize(k);
        x.resize(k-1);
        y.resize(k-1);
    }

    for (int i = 0; i < k - 1; i++) {
        x[i] = primal.slope[i];
        y[i] = primal.y[i] - primal.slope[i] * primal.x[i];
    }
    for (int i = 0; i < k; i++) {
        slope[i] = -primal.x[i];
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

};
