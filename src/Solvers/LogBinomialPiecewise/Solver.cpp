//
// Created by Yuanyuan Qi on 6/12/24.
//

#include <cstdio>
#include <functional>
#include <algorithm>

#include "Solver.h"
#include "Tree.h"
#include "Node.h"

namespace LogBinomialPiecewiseLinearSolver {

Solver::Solver(int n_intervals):n_intervals(n_intervals),n(T.size) {
}

void Solver::init(const std::vector<int> &var, const std::vector<int> &ref,
               const std::vector<std::list<int> > &llist, int root) {
    T.init(llist, root);
    for (int i = 0; i < n; i++) {
        T.nodes[i].step_ptr = &duals[i];
        T.nodes[i].state_ptr = &states[i];
        T.nodes[i].sum_of_children = &sums[i];
        funcs[i].var=(var[i]);
        funcs[i].ref=(ref[i]);
        F[i] = 0.5;
    }
}

void Solver::init_range(std::unordered_map<int, real> &mid, real fu, double eps){
    for (int i = 0; i < n; i++){
        if (funcs[i].var > 0) {
            primals[i].update(std::max(mid[i]-fu, eps),std::min(mid[i]+fu,1.),n_intervals,&funcs[i]);
        }
        else if (funcs[i].ref > 0) {
            primals[i].update(std::max(mid[i]-fu, 0.),std::min(mid[i]+fu,1. - eps),n_intervals,&funcs[i]);
        }
        else {
            primals[i].update(std::max(mid[i]-fu,0.),std::min(mid[i]+fu,1.),n_intervals,&funcs[i]);
        }
        duals[i].dual_from_primal(primals[i]);
    }
}

real Solver::answer() {
    int k = binary_search_greater(states[T.root].slope,states[T.root].k,1.);
    // this should just simply be zero all the time.
    real answer = states[T.root].y[k];
    dual_0 =  states[T.root].x[k];
    return answer;
}

void Solver::solve_F() {
    for (int i = 0; i < n; i++){
        if (i == T.root){
            F[i] = primals[i].backtrace_delta(T.BT_vals[i] - dual_0);
        } else {
            F[i] = primals[i].backtrace_delta(T.BT_vals[i] - T.BT_vals[T.nodes[i].parent]);
        }
    }
//    f primal[node].backtrace(dual_vars[node] - value );
}

real Solver::solve(real eps) {
     real  ans;
     init_range(F,0.51, eps);
     T.dfs(T.root);
     ans = answer();

     T.dfs_BT(T.root, dual_0);
     solve_F();
     return ans;
}

real Solver::solve_iteratively(real frac, real obj) {
    // real range = 1, ans;
    // init_range(F,0.51);
    // T.dfs(T.root);
    // ans = answer();
    // 
    // solve_F();
    // return ans;

    real range = 1, ans;
    init_range(F,0.51);

    do  {
        T.dfs(T.root);
        ans = answer();
        T.dfs_BT(T.root, dual_0);
        solve_F();
        range = range*frac;
        init_range(F,range/2);
    } while (range/n_intervals > obj);

    return ans;
}

//real Solver::move(int from, int to) {
//    T.move(from, to);
//    real ans = answer();
//    T.dfs_BT(T.root, dual_0);
//    solve_F();
//    return ans;
//}

};
