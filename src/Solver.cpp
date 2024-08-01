//
// Created by Yuanyuan Qi on 6/12/24.
//

#include "Solver.h"
#include <stack>
#include <cassert>

Solver::Solver(int max_n, int k):
F(max_n),
dual_vars(max_n),
n_intervals(k),
funcs(max_n),
primal(max_n,k),
dual(max_n,k+1),
sum(max_n,(max_n+1)*(k+1)),
states(max_n, (max_n+1)*(k+1)),
children(max_n),
children_state(max_n),
F_min(max_n),
F_scm(max_n),
step(max_n),
all_(max_n)
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
        primal[i].update(std::max(mid[i]-fu,0.),std::min(mid[i]+fu,1.),n_intervals,funcs[i]);
        assert(primal[i].self_check());
        dual[i].dual_from_primal(primal[i]);
        assert(dual[i].self_check());
    }
}

void Solver::dfs(int node) {
//    if (node == 201){
//        system("pause");
//    }
    if (children[node].empty()){
        states[node].base_case(dual[node]);
    }
    else {
        for (auto ch:children[node])
            dfs(ch);
        sum[node].sum(children_state[node], helper);
        assert(sum[node].self_check());
        states[node].optimize_with(sum[node],dual[node]);
        assert(states[node].self_check());
    }
    if (abs(sum[node].slope[sum[node].k-1]-F_scm[node]) > Compare_eps){
        printf("%d : %lf %lf\n",node,sum[node].slope[sum[node].k-1], F_scm[node]);
        printf("ERROR\n");
    }
}

void Solver::dfs_BT(int node, real value) {
    real check1,check2;
    if (children[node].empty()) {
        dual_vars[node] = 0;
        check1 = dual[node](-value);
    }
    else {
        dual_vars[node] = sum[node].backtrace(dual[node],value, &check1);
    }
    check2 = states[node](value);
    real min_val = primal[node].backtrace(dual_vars[node] - value, &step[node] ), sum = 0;
    all_[node]=0;
    for (auto ch: children[node]) {
        dfs_BT(ch, dual_vars[node]);
        all_[node]+=all_[ch];
    }
    all_[node]+=step[node];
    F[node] = min_val;
    real a=std::min(std::min(check1,check2),all_[node]), b=std::max(std::max(check1,check2),all_[node]);
    if (b-a>1e-6){
        printf("%d: %lf %lf %lf\n",node, check1, check2, all_[node]);
        printf("ERROR\n");
    }
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

    my_assertions(this->F, 0.51, root);

    real range = 1, ans;
    do  {
            dfs(this->root);

            ans = answer();

            printf("iteration: obj=%.12lf\n",ans);

            dfs_BT(root, dual_0);

            range = range*frac;

            this->init_range(this->F,range/2);

            my_assertions(this->F, range/2, root);

        }while(range/n_intervals > obj);

    return ans;
}


void Solver::my_assertions(std::vector<real> &mid, real fu, int node){
    F_scm[node] = 0;
    for (auto ch: children[node]){
        my_assertions(mid,fu,ch);
        F_scm[node]+=std::max(F_min[ch],F_scm[ch]);
    }
    F_min[node] = std::max(0.,mid[node]-fu);
    if (mid[node]+fu < std::max(F_min[node],F_scm[node])){
        printf("Range is problematic!\n");
    }
}
