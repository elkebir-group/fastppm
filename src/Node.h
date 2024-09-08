//
// Created by Yuanyuan Qi on 8/29/24.
//

#ifndef FASTPPM_NODE_H
#define FASTPPM_NODE_H

#include<unordered_set>
#include<unordered_map>

#include "Func.h"

class step_base;
class state_base;
class sum_base;

class Node{
public:
    int parent;
    std::unordered_set<int> children;
    // please use already created space, do not delete
    step_base * step_ptr;
    state_base * state_ptr;
    sum_base * sum_of_children;
    std::unordered_map<int, Node> * set_of_all; //point this to the vector of nodes;
    void sum();
    void optimize() const;
    void base_case() const;
    real backtrace(real x) const;
    real backtrace_base (real x) const;
    std::unordered_set<int> children_to_update;
};

#endif //FASTPPM_NODE_H
