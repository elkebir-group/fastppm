//
// Created by Yuanyuan Qi on 6/24/24.
//

#ifndef EFFICIENTLLHESTIMATOR_PWL_H
#define EFFICIENTLLHESTIMATOR_PWL_H

#include <unordered_map>

#include "Func.h"
#include "Element.h"

class Primal;

class Dual;

class State:public state_base,public sum_base{// x bounded on left side // left is always 0. i.e. x[0]=0 // state function
public:
    int k; // n of intervals;
    std::unordered_map<int, real> x;
    std::unordered_map<int, real> y;
    std::unordered_map<int, real> slope;
    static std::unordered_map<int, std::pair<real,real> > helper;
    static std::unordered_map<int, real> helper_x, helper_y;

    void base_case(const step_base * dual) override;
//    void base_case(const Dual & dual);
    void sum(const std::vector<const state_base *> & to_sum) override;
    void optimize_with(const sum_base * sum_func, const step_base * dual) override;
//    void optimize_with(const State & sum_func, const Dual & dual);
    real backtrace(const step_base * dual, real val) const override;

#ifdef _DEBUG
    real operator()(real x) const;
    bool self_check() const;
#endif
};

#endif //EFFICIENTLLHESTIMATOR_PWL_H
