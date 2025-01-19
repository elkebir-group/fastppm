//
// Created by Yuanyuan Qi on 6/12/24.
//

#ifndef EFFICIENTLLHESTIMATOR_SOLVER_H
#define EFFICIENTLLHESTIMATOR_SOLVER_H

#include <list>

#include "Func.h"
#include "State.h"
#include "Primal.h"
#include "Dual.h"
#include "Tree.h"

namespace PiecewiseLinearSolver {

class func_base;
class func_binom_llh;
class func_beta_binom_llh;

class Solver {
public:
    int n_intervals;
    const int &n;
    real dual_0;
//    std::unordered_map<int, func_binom_llh> funcs;
    std::unordered_map<int, Primal> primals;
    std::unordered_map<int, Dual> duals;
    std::unordered_map<int, State> states;
    std::unordered_map<int, State> sums;
    std::unordered_map<int, real> F;
    Tree T;

    Solver(int n_intervals);
    void init(const std::vector<int> &var, const std::vector<int> &ref,
                      const std::vector<std::list<int> > &llist, int root);
    virtual void init_range(std::unordered_map<int, real> & mid, real range, real eps=0.);
    real answer();
    void solve_F();
    real solve_iteratively(real frac=0.75, real obj=1e-6, real eps=0.);
    real solve(real eps=0.);
    real move(int from, int to);

    virtual func_base& func_i(int i) = 0;
};

class SolverBinomial : public Solver {
public:
    std::unordered_map<int, func_binom_llh> funcs;

    SolverBinomial(int n_intervals)
    : Solver(n_intervals)
    {
    }

    virtual func_base& func_i(int i) {
        return funcs[i];
    }
};

class SolverBetaBinomial : public Solver {
public:
    std::unordered_map<int, func_beta_binom_llh> funcs;
    real precision;

    SolverBetaBinomial(int n_intervals, double precision)
    : Solver(n_intervals)
    , precision(precision)
    {
        func_beta_binom_llh::precision = precision;
    }

    virtual void init_range(std::unordered_map<int, real> & mid, real range, real eps=0.);

    virtual func_base& func_i(int i) {
        return funcs[i];
    }
};

};
#endif //EFFICIENTLLHESTIMATOR_SOLVER_H
