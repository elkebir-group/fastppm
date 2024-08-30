//
// Created by Yuanyuan Qi on 8/1/24.
//

#ifndef EFFICIENTLLHESTIMATOR_DUAL_H
#define EFFICIENTLLHESTIMATOR_DUAL_H

#include "Func.h"
#include <vector>

namespace LogBinomialPiecewiseLinearSolver {

class PWL_close;
typedef PWL_close Primal;

class PWL_open{//x not bounded. // dual function
public:
    PWL_open(int max_k);
    int k; // n of intervals;
    std::vector<real> x;
    std::vector<real> y;
    std::vector<real> slope;

    void dual_from_primal(const Primal & primal);

#ifdef _DEBUG
    real operator()(real x) const;
    bool self_check() const;
#endif
};

typedef PWL_open Dual;

};

#endif //EFFICIENTLLHESTIMATOR_DUAL_H
