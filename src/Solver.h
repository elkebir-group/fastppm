//
// Created by Yuanyuan Qi on 6/12/24.
//

#ifndef EFFICIENTLLHESTIMATOR_SOLVER_H
#define EFFICIENTLLHESTIMATOR_SOLVER_H
#include <list>

#include "Func.h"
#include "PWL.h"
#include "Primal.h"
#include "Dual.h"

class Solver {
public:
    int n_intervals, n, root;
    real dual_0;
    std::vector<std::pair<real,real> > helper;
    std::vector<Func> funcs;
    std::vector<Primal> primal;
    std::vector<Dual> dual;
    std::vector<PWL_Ropen> states;
    std::vector<PWL_Ropen> sum;

    std::vector<real> dual_vars;
    std::vector<real> F;//
    std::vector<std::vector<int> > children;
    std::vector<std::vector<PWL_Ropen*> > children_state;

    std::vector<real> BT_x, BT_y;

    Solver(int max_n, int k);
    void init(const std::vector<int> &var, const std::vector<int> &ref,
              const std::vector<std::list<int> > &llist, int root);
    void init_range(std::vector<real> & mid, real range);

    void dfs(int node);
    void dfs_BT(int node, real value);

    real answer();
    real main(real frac=0.75, real obj=1e-6);
};


#endif //EFFICIENTLLHESTIMATOR_SOLVER_H
