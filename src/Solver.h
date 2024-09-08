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

class func_llh;

class Solver {
public:
    int n_intervals;
    const int &n;
    real dual_0;
    std::unordered_map<int, func_llh> funcs;
    std::unordered_map<int, Primal> primals;
    std::unordered_map<int, Dual> duals;
    std::unordered_map<int, State> states;
    std::unordered_map<int, State> sums;
    std::unordered_map<int, real> F;
    Tree T;

    Solver(int n_intervals);
    void init(const std::vector<int> &var, const std::vector<int> &ref,
              const std::vector<std::list<int> > &llist, int root);
    void init_range(std::unordered_map<int, real> & mid, real range);
    real answer();
    void solve_F();
    real main(real frac=0.75, real obj=1e-6);
    real move(int from, int to);
};


#endif //EFFICIENTLLHESTIMATOR_SOLVER_H
