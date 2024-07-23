//
// Created by Yuanyuan Qi on 6/12/24.
//

#include "Solver.h"
#include <stack>
#include <iostream>

Node::Node(int n_intervals):
        pwl(n_intervals),
        dpstate(),
        size(0)
{
}

Solver::Solver(const std::vector<int> &var, const std::vector<int> &ref, int n_intervals,
               const std::vector<std::list<int> > &llist, int root):
nodes(var.size(), n_intervals),
dual_vars(var.size()),
F(var.size()),
n_intervals(n_intervals),
root(root)
{
    std::stack<int> Stack1, Stack2;
    Stack1.push(root); int top;
    while (!Stack1.empty()){
        top = Stack1.top();
        Stack1.pop(); Stack2.push(top);
        for (auto c:llist[top]){
            Stack1.push(c);
        }
    }
    while (!Stack2.empty()){
        top = Stack2.top();
        Stack2.pop();
        nodes[top].size = 1;
        for (auto c:llist[top]){
            nodes[top].size += nodes[c].size;
        }
    }

    for (int i = 0,ch_idx; i < var.size(); i++){
        nodes[i].func.var = var[i];
        nodes[i].func.ref = ref[i];
        nodes[i].dpstate.init(nodes[i].size,n_intervals);
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

void Solver::init_range(std::vector<real> &mid, real fu){
    for (int i = 0; i < nodes.size(); i++){
        nodes[i].pwl.update(mid[i]-fu,mid[i]+fu,nodes[i].func);
    }
}

void Solver::dfs(int node) {
    for (auto ch:nodes[node].children)
        dfs(ch);
    printf("-----node %d-----\n",node);
    nodes[node].dpstate.update(nodes[node].c_state,nodes[node].pwl, node ==5 || node==32);
}

void Solver::dfs_BT(int node, real value) {
//    if (nodes[node].children.empty()) dual_vars[node] = 0;
//    else
    dual_vars[node] = nodes[node].dpstate.backtrace(value);
    F[node] = 0;
    for (auto ch: nodes[node].children){
        dfs_BT(ch,dual_vars[node]);
        F[node] += F[ch];
    }
    F[node] = std::max(F[node],nodes[node].pwl.bt_f(dual_vars[node]-value));
}

real Solver::answer() {
    int k = std::lower_bound(nodes[root].dpstate.slope.begin(), nodes[root].dpstate.slope.end(), 1,
                             std::greater<real>())
            - nodes[root].dpstate.slope.begin();
    real answer = nodes[root].dpstate.y[0];
    for (int i = 1; i <= k; i++) {
        answer += nodes[root].dpstate.slope[i-1]* (nodes[root].dpstate.x[i]-nodes[root].dpstate.x[i-1]);
    }
    dual_0 = nodes[root].dpstate.x[k];
    return answer;
}

void Solver::backtrace() {
    dfs_BT(root, dual_0);
}

real Solver::main() {

    this->init_range(0,1);

    real ans;

    dfs(this->root);

    ans = answer();

//    backtrace();

    return ans;
}
