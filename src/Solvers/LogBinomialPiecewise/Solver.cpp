//
// Created by Yuanyuan Qi on 6/12/24.
//

#include "Solver.h"
#include <stack>
#include <cstdio>
#include <functional>
#include <algorithm>

namespace LogBinomialPiecewiseLinearSolver {

Solver::Solver(int max_n, int k):
F(max_n),
dual_vars(max_n),
n_intervals(k),
funcs(max_n),
primal(max_n,k),
dual(max_n,k+1),
sum(max_n,(max_n+1)*(k+1)),
states(max_n, (max_n+1)*(k+1)),
helper((max_n+1)*(k+1)),
children(max_n),
children_state(max_n),
BT_x((max_n+1)*(k+1)),
BT_y((max_n+1)*(k+1)),
loss_func(func_llh)
{
}

void Solver::init(const std::vector<int> &var, const std::vector<int> &ref,
               const std::vector<std::list<int> > &llist, int root) {
    n = var.size();
    this->root = root;
    for (int i = 0,idx; i < n; i++) {
        children[i] = std::vector<int>(llist[i].begin(),llist[i].end());
        children_state[i].resize(llist[i].size());
        idx=0;
        for (auto j:llist[i]) {
            children_state[i][idx++] = &states[j];
        }
        funcs[i] = {real(var[i]), real(ref[i])};
    }
    std::fill(F.begin(), F.end(), 0.5);
}

void Solver::init_range(std::vector<real> &mid, real fu){
    for (int i = 0; i < n; i++){
        primal[i].update(std::max(mid[i]-fu,0.),std::min(mid[i]+fu,1.),n_intervals,funcs[i], loss_func);
//        assert(primal[i].self_check());
        dual[i].dual_from_primal(primal[i]);
//        assert(dual[i].self_check());
    }
}

void Solver::dfs(int node) {
    if (children[node].empty()){
        states[node].base_case(dual[node]);
    }
    else {
        for (auto ch:children[node])
            dfs(ch);
        sum[node].sum(children_state[node], helper);
//        assert(sum[node].self_check());
        states[node].optimize_with(sum[node],dual[node]);
//        assert(states[node].self_check());
    }
}

void Solver::dfs_BT(int node, real value) {
    if (children[node].empty()) {
        dual_vars[node] = 0;
    }
    else {
        dual_vars[node] = sum[node].backtrace(dual[node],value,BT_x,BT_y);
    }
    real min_val = primal[node].backtrace(dual_vars[node] - value ), sum = 0;
    for (auto ch: children[node]) {
        dfs_BT(ch, dual_vars[node]);
    }
    F[node] = min_val;
}

real Solver::answer() {
    int k = std::lower_bound(states[root].slope.begin(), states[root].slope.begin()+states[root].k, 1,
                             std::greater<real>())
            - states[root].slope.begin(); // this should just simply be zero all the time.
    real answer = states[root].y[k];
    dual_0 =  states[root].x[k];
    return answer;
}


real Solver::main(real frac, real obj) {

    this->init_range(F,0.51);

    real range = 1, ans;
    do  {
            dfs(this->root);

            ans = answer();

            //printf("iteration: obj=%.12lf\n",ans); 

            dfs_BT(root, dual_0);

            range = range*frac;

            this->init_range(this->F,range/2);

        }while(range/n_intervals > obj);

    return ans;
}

Solver::~Solver() {
}

};
