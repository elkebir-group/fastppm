//
// Created by Yuanyuan Qi on 8/29/24.
//

#ifndef FASTPPM_TREE_H
#define FASTPPM_TREE_H

#include<vector>
#include<list>
#include<unordered_map>

#include "Func.h"
#include "Node.h"

class step_base;
class sum_base;
class state_base;

class Tree {
public:
    std::unordered_map<int, step_base*> steps; //dual
    std::unordered_map<int, sum_base*> sums;
    std::unordered_map<int, state_base*> states;
    std::unordered_map<int, Node> nodes;
    std::unordered_map<int, real> BT_vals; //multiple

    int size, root;

    void init();
    void init(const std::vector<std::list<int> > & llist, int root);

    void dfs(int node);

    void dfs_BT(int node, real value);

    void mark(int node);

    void dfs_update(int node);

    void move(int from, int to);
};


#endif //FASTPPM_TREE_H
