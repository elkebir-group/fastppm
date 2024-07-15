//
// Created by Yuanyuan Qi on 6/12/24.
//

#ifndef EFFICIENTLLHESTIMATOR_SOLVER_H
#define EFFICIENTLLHESTIMATOR_SOLVER_H
#include <list>

#include "Func.h"
#include "DPSTATE.h"
#include "PWL.h"

struct Node {
    Func func;
    PWL pwl;
    int size;
    DPSTATE dpstate;
    std::vector<int> children;
    std::vector<DPSTATE*> c_state;
    Node(int n_intervals);
};

class Solver {
public:
    DPSTATE *r_dpstate;
    int root; real dual_0;
    int n_intervals;
    std::vector<Node> nodes;
    std::vector<real> dual_vars;
    std::vector<real> F;//
    Solver(const std::vector<int> & var, const std::vector<int> & ref, int n_intervals,
           const std::vector<std::list<int> > & llist, int root);
    void init_range(std::vector<real> & fl, std::vector<real> & fu);
    void init_range(real fl, real fu);
    void init_range(std::vector<real> & mid, real range);
    void dfs(int node);
    void dfs_BT(int node, real value);
    real answer();
    void backtrace();
    real main();
};


#endif //EFFICIENTLLHESTIMATOR_SOLVER_H
