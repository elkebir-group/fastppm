//
// Created by Yuanyuan Qi on 8/29/24.
//

#include "Node.h"
#include "Element.h"

namespace PiecewiseLinearSolver {
void Node::sum() {
    std::vector<const state_base *> to_sum(children.size());
    int i=0;
    for (auto a:children){
        to_sum[i++] = (*set_of_all)[a].state_ptr;
    }
    sum_of_children -> sum(to_sum);
}

void Node::optimize() const {
    state_ptr -> optimize_with(sum_of_children, step_ptr);
}

void Node::base_case() const {
    state_ptr -> base_case(step_ptr);
}

real Node::backtrace(real x) const{
    return sum_of_children -> backtrace(step_ptr, x);
}

real Node::backtrace_base(real x) const{
    return step_ptr -> backtrace_base(x);
}
};
