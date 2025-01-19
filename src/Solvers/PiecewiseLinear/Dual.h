//
// Created by Yuanyuan Qi on 8/1/24.
//

#ifndef EFFICIENTLLHESTIMATOR_DUAL_H
#define EFFICIENTLLHESTIMATOR_DUAL_H

#include "Func.h"
#include "Element.h"

namespace PiecewiseLinearSolver {
class Primal;
class Dual:public step_base{//x not bounded. // dual function
public:
    int k; // n of intervals;
    std::unordered_map<int, real> x;
    std::unordered_map<int, real> y;
    std::unordered_map<int, real> slope;

    void dual_from_primal(const Primal & primal);

    inline real backtrace_base (real x) const override{
        return 0;
    }

#ifdef _DEBUG
    bool self_check() const;
#endif
};
};

#endif //EFFICIENTLLHESTIMATOR_DUAL_H
