//
// Created by Yuanyuan Qi on 8/31/24.
//

#ifndef FASTPPM_ELEMENT_H
#define FASTPPM_ELEMENT_H

#include <vector>

namespace PiecewiseLinearSolver {
class step_base{// Dual func
public:
    virtual real backtrace_base(real x) const = 0;
};

class sum_base;
class state_base {//dp state function
public:
    virtual void base_case(const step_base *) = 0;
    virtual void optimize_with(const sum_base *, const step_base *) = 0;
};

class sum_base{// sum of dp state functions
public:
    virtual void sum(const std::vector<const state_base*> &) = 0;
    virtual real backtrace(const step_base *, real x) const = 0;
};
};

#endif //FASTPPM_ELEMENT_H
