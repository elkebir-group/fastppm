//
// Created by Yuanyuan Qi on 6/25/24.
//

#ifndef EFFICIENTLLHESTIMATOR_FUNC_H
#define EFFICIENTLLHESTIMATOR_FUNC_H
#include <cmath>
#include <vector>
#include <functional>

typedef double real;

#ifdef _DEBUG
#include <cassert>
const real assert_eps = 1.5e-6;
#endif

const real Compare_eps = 1e-7;
//const real
real log_eps(real x,real eps=1e-6,int s_n=3);

inline real func_llh(real para_var, real para_ref, real x) {
    return -para_var * log_eps(x) - para_ref * log_eps(1 - x);
}

#endif //EFFICIENTLLHESTIMATOR_FUNC_H
