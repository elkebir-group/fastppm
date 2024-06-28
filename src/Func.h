//
// Created by Yuanyuan Qi on 6/25/24.
//

#ifndef EFFICIENTLLHESTIMATOR_FUNC_H
#define EFFICIENTLLHESTIMATOR_FUNC_H
#include <cmath>
typedef double real;
real log_eps(real x,real eps=1e-6,int s_n=3);

struct Func {
    real var, ref;
    inline real operator()(real x) const {
        return var * log_eps(x) + ref * log_eps(1 - x);
    }
};

#endif //EFFICIENTLLHESTIMATOR_FUNC_H
