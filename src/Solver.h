//
// Created by Yuanyuan Qi on 6/12/24.
//

#ifndef EFFICIENTLLHESTIMATOR_SOLVER_H
#define EFFICIENTLLHESTIMATOR_SOLVER_H
#include <list>

#include "Func.h"
#include "PWL.h"

class Solver {
public:
    int n_intervals, n, root;
    real dual_0;
    std::vector<Func> funcs;
    std::vector<PWL_close> primal;
    std::vector<PWL_open> dual;
    std::vector<PWL_Ropen> states;
    std::vector<PWL_Ropen> sum;

    std::vector<real> dual_vars;
    std::vector<real> F;//
    std::vector<std::vector<int> > children;
    std::vector<std::vector<PWL_Ropen*> > children_state;
    Solver(int max_n, int k);
    void init(const std::vector<int> &var, const std::vector<int> &ref,
              const std::vector<std::list<int> > &llist, int root);
    void init_range(std::vector<real> & mid, real range);

    void dfs(int node);
    void dfs_BT(int node, real value);

    real answer();
    real main(real frac=0.75, real obj=1e-6);

    ////////debug purpose only
    std::vector<real> F_min,F_scm;
    std::vector<real> step,all_;
    void my_assertions(std::vector<real> &mid, real fu, int node);
};


#endif //EFFICIENTLLHESTIMATOR_SOLVER_H
