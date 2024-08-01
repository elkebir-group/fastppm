//
// Created by Yuanyuan Qi on 8/1/24.
//

#ifndef EFFICIENTLLHESTIMATOR_PRIMAL_H
#define EFFICIENTLLHESTIMATOR_PRIMAL_H

#include "Func.h"

#include <vector>

class PWL_close {// x Bounded on both end. //primal function
public:
    PWL_close(int max_k);
    int k; // n of intervals;
    std::vector<real> x;
    std::vector<real> y;
    std::vector<real> slope;

    void update(real begin, real end, int _k, const Func & func);
    real backtrace(real d, real * debug = NULL);

#ifdef _DEBUG
    real operator()(real x) const;
    bool self_check() const;
#endif
};

typedef PWL_close Primal;

#endif //EFFICIENTLLHESTIMATOR_PRIMAL_H
