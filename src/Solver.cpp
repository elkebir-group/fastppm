//
// Created by Yuanyuan Qi on 6/12/24.
//

#include "Solver.h"
#include <iostream>

Node::Node(int n_intervals):
        pwl(n_intervals),
        dpstate()
{
}

Solver::Solver(const std::vector<int> &var, const std::vector<int> &ref, int n_intervals,
               const std::vector<std::list<int> > &llist, int root):
nodes(var.size(), n_intervals),
root(root)
{
    for (int i = 0,ch_idx; i < var.size(); i++){
        nodes[i].func.var = var[i];
        nodes[i].func.ref = ref[i];
        nodes[i].dpstate.init(llist[i].size(),n_intervals);
        nodes[i].c_state.resize(llist[i].size());
        nodes[i].children.resize(llist[i].size());
        ch_idx=0;
        for (auto ch:llist[i]){
            nodes[i].c_state[ch_idx] = & nodes[ch].dpstate;
            nodes[i].children[ch_idx] = ch;
            ch_idx++;
        }
    }
}

void Solver::init_range(std::vector<real> & fl, std::vector<real> & fu){
    for (int i = 0; i < nodes.size(); i++){
        nodes[i].pwl.update(fl[i],fu[i],nodes[i].func);
    }
}

void Solver::init_range(real fl, real fu){
    for (int i = 0; i < nodes.size(); i++){
        nodes[i].pwl.update(fl,fu,nodes[i].func);
    }
}

void Solver::dfs(int node) {
    for (auto ch:nodes[node].children)
        dfs(ch);
    nodes[node].dpstate.update(nodes[node].c_state,nodes[node].pwl);
}

real Solver::answer() {
    int k = std::lower_bound(nodes[root].dpstate.slope.begin(), nodes[root].dpstate.slope.end(), 1,
                             std::greater<real>())
            - nodes[root].dpstate.slope.begin();
    std::cout << "at " << k << std::endl;
    real answer = nodes[root].dpstate.y[0], cx = 0;
    for (int i = 1; i <= k; i++) {
        answer += nodes[root].dpstate.slope[i-1]* (nodes[root].dpstate.x[i]-nodes[root].dpstate.x[i-1]);
    }
    return answer;
}

void Solver::backtrace() {

}
