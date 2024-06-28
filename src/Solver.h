//
// Created by Yuanyuan Qi on 6/12/24.
//

#ifndef EFFICIENTLLHESTIMATOR_SOLVER_H
#define EFFICIENTLLHESTIMATOR_SOLVER_H
#include <list>

#include "Func.h"
#include "PWL.h"

struct Node {
    Func func;
    PWL pwl;
    DPSTATE dpstate;
    std::vector<int> children;
    std::vector<DPSTATE*> c_state;
    Node(int n_intervals);
};

class Solver {
public:
    DPSTATE *r_dpstate;
    int root;
    std::vector<Node> nodes;
    Solver(const std::vector<int> & var, const std::vector<int> & ref, int n_intervals,
           const std::vector<std::list<int> > & llist, int root);
    void init_range(std::vector<real> & fl, std::vector<real> & fu);
    void init_range(real fl, real fu);
    void dfs(int node);
    real answer();
    void backtrace();
};


#endif //EFFICIENTLLHESTIMATOR_SOLVER_H
