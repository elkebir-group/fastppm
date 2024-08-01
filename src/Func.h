//
// Created by Yuanyuan Qi on 6/25/24.
//

#ifndef EFFICIENTLLHESTIMATOR_FUNC_H
#define EFFICIENTLLHESTIMATOR_FUNC_H
#include <cmath>
#include <vector>
typedef double real;
const real r_zero=0;
const real r_one=1;
extern std::vector<std::pair<real,real> > helper; // n, k
extern std::vector<real> _d_y0;
const real Compare_eps = 1e-7;
#ifdef _DEBUG
const real assert_eps = 1.5e-6;
#endif
//const real
real log_eps(real x,real eps=1e-6,int s_n=3);

struct Func {
    real var, ref;
    inline real operator()(real x) const {
        return - var * log_eps(x) - ref * log_eps(1 - x);
    }
};

#endif //EFFICIENTLLHESTIMATOR_FUNC_H
