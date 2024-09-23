//
// Created by Yuanyuan Qi on 8/29/24.
//

#include "Tree.h"
#include "Node.h"
#include "Element.h"

namespace LogBinomialPiecewiseLinearSolver {

void Tree::init() {
    size = 0;
    root = -1;
}

void Tree::init(const std::vector<std::list<int> > & llist, int root) {
    size = llist.size();
    this->root = root;
    nodes[root].parent = -1;
    for (int i = 0; i < size; i++) {
        nodes[i].children = std::move(std::unordered_set<int>(llist[i].begin(), llist[i].end()));
        for (auto ch: nodes[i].children) {
            nodes[ch].parent = i;
        }
        nodes[i].set_of_all = &nodes;
    }
}

void Tree::dfs(int node) {
    if (nodes[node].children.empty()){
        nodes[node].base_case();
    } else {
        for (auto ch:nodes[node].children){
            dfs(ch);
        }
        nodes[node].sum();
        nodes[node].optimize();
    }
}

void Tree::dfs_BT(int node, real value) { // returns all dual value
    if (nodes[node].children.empty()) {
        BT_vals[node]=nodes[node].backtrace_base(value);
    } else {
        BT_vals[node]=nodes[node].backtrace(value);
    }
    for (auto ch:nodes[node].children){
        dfs_BT(ch, BT_vals[node]);
    }
}

void Tree::dfs_update(int node){
    if (nodes[node].children.empty()){
        nodes[node].base_case();
    } else {
        for (auto ch: nodes[node].children_to_update){
            dfs_update(ch);
        }
        nodes[node].children_to_update.clear();
        nodes[node].sum();
        nodes[node].optimize();
    }
}

void Tree::mark(int node) {
    for (int i = node; i != root; i = nodes[i].parent){
        nodes[nodes[i].parent].children_to_update.insert(i);
    }
}

//void Tree::move(int from, int to){
//    int parent_of_from = nodes[from].parent;
//    nodes[parent_of_from].children.erase(from);
//    nodes[to].children.insert(from);
//    nodes[from].parent = to;
//    mark(parent_of_from);
//    mark(to);
//    dfs_update(root);
//}
//
//void Tree::swap(int a, int d) {
//    if (a != root) nodes[nodes[a].parent].children.erase(a);
//    nodes[nodes[d].parent].children.erase(d);
//    std::swap(nodes[a].children, nodes[d].children);
//    if (a != root) nodes[nodes[a].parent].children.insert(d);
//    std::swap(nodes[a].parent, nodes[d].parent); // this has to happen before the next insert considering the case where a is parent of d
//    nodes[nodes[d].parent].children.insert(a);
//    if (a == root) root = d;
//    mark(a);
//    dfs_update(root);
//}
//
//void Tree::insert_as_leaf(int parent, step_base* step, sum_base* sum, state_base* state) {
//    nodes[size].parent = parent;
//    nodes[parent].children.insert(size);
//
//    nodes[size].set_of_all = &nodes;
//    nodes[size].step_ptr = step;
//    nodes[size].sum_of_children = sum;
//    nodes[size].state_ptr = state;
//    size++;
//    mark(size - 1);
//    dfs_update(root);
//}
//
//
//void Tree::insert_on_edge(int child, step_base *step, sum_base *sum, state_base *state) {
//    nodes[size].parent = nodes[child].parent;
//    if (child!=root) {
//        nodes[nodes[child].parent].children.erase(child);
//        nodes[nodes[child].parent].children.insert(size);
//    }
//    else{
//        root = size;
//    }
//    nodes[size].children.insert(child);
//    nodes[child].parent = size;
//
//    nodes[size].set_of_all = &nodes;
//    nodes[size].step_ptr = step;
//    nodes[size].sum_of_children = sum;
//    nodes[size].state_ptr = state;
//    size++;
//    mark(size - 1);
//    dfs_update(root);
//}

};
