//
// Created by Yuanyuan Qi on 8/1/24.
//

#ifndef EFFICIENTLLHESTIMATOR_PRIMAL_H
#define EFFICIENTLLHESTIMATOR_PRIMAL_H

#include "Func.h"

#include <unordered_map>

namespace LogBinomialPiecewiseLinearSolver {
class Primal {// x Bounded on both end. //primal function
public:
    int k; // n of intervals;
    std::vector<real> x;
    std::vector<real> y;
    std::vector<real> slope;

    void update(real begin, real end, int _k, func_base * func);

    void update(const std::vector<real> & xx, const std::vector<real> & yy);

    real backtrace_delta(real d);

#ifdef _DEBUG
    bool self_check() const;
#endif
};
};
#endif //EFFICIENTLLHESTIMATOR_PRIMAL_H
