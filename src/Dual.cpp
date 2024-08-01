//
// Created by Yuanyuan Qi on 8/1/24.
//

#include "Dual.h"
#include "Primal.h"


PWL_open::PWL_open(int max_k) :
        k (0),
        x(max_k - 1),
        y(max_k - 1),
        slope(max_k) {
}

void PWL_open::dual_from_primal(const Primal &primal) {
    k = primal.k + 1;
    for (int i = 0; i < k - 1; i++) {
        x[i] = primal.slope[i];
        y[i] = primal.y[i] - primal.slope[i] * primal.x[i];
    }
    for (int i = 0; i < k; i++) {
        slope[i] = -primal.x[i];
    }
}

#ifdef _DEBUG

real PWL_open::operator()(real xx) const {
    int idx = std::lower_bound(x.begin(),x.begin()+k-1,xx)-x.begin();
    return y[idx]+(xx-x[idx])*slope[idx];
}


bool PWL_open::self_check() const {
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