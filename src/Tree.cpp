//
// Created by Yuanyuan Qi on 8/29/24.
//

#include "Tree.h"
#include "Node.h"
#include "Element.h"


void Tree::init() {
    size = 0;
    root = -1;
}

void Tree::init(const std::vector<std::list<int> > & llist, int root) {
    size = llist.size();
    this->root = root;
    nodes[root].parent=-1;
    for (int i = 0; i < size; i++) {
        nodes[i].children = std::move(std::unordered_set<int>(llist[i].begin(),llist[i].end()));
        for (auto ch:nodes[i].children){
            nodes[ch].parent = i;
        }
        nodes[i].set_of_all = & nodes;
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

void Tree::move(int from, int to){ // not defensive programming assuming from is never an ancestor of to
    int parent_of_from = nodes[from].parent;
    nodes[parent_of_from].children.erase(from);
    nodes[to].children.insert(from);
    nodes[from].parent = to;
    mark(parent_of_from);
    mark(to);
    dfs_update(root);
}
