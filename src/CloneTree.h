#ifndef CLONETREE_H
#define CLONETREE_H

#include "DiGraph.h"

#include <vector>
#include <unordered_map>
#include <stack>

inline std::vector<int> postorder_dfs(const digraph<int>& clone_tree, int root) {
    std::vector<int> postorder;
    std::set<int> visited;
    std::stack<int> stack;
    stack.push(root);
    while (!stack.empty()) {
        int vertex = stack.top();
        stack.pop();

        bool all_children_visited = true;
        for (auto child : clone_tree.successors(vertex)) {
            if (!visited.contains(child)) all_children_visited = false;
        }

        if (all_children_visited) {
            postorder.push_back(vertex);
            visited.insert(vertex);
        } else {
            stack.push(vertex);
            for (auto child : clone_tree.successors(vertex)) {
                stack.push(child);
            } 
        }
    }

    return postorder;
}

/* 
 * Computes the left inverse product F*B^{-1} where 
 * B is the n-clonal matrix B corresponding to T in 
 * O(mn) time where F is m x n.
 */
inline std::vector<std::vector<double>> left_inverse(
    const digraph<int>& clone_tree,
    const std::unordered_map<int, int>& vertex_map,
    const std::vector<std::vector<double>>& frequency_matrix
) {
    std::vector<std::vector<double>> usage_matrix = frequency_matrix;
    for (size_t j = 0; j < frequency_matrix.size(); j++) {
        for (size_t i = 0; i < frequency_matrix[j].size(); i++) {
            for (auto child : clone_tree.successors(vertex_map.at(i))) {
                size_t l = clone_tree[child].data;
                usage_matrix[j][i] -= frequency_matrix[j][l];
            }
        }
    }
    return usage_matrix;
}

/* 
 * Computes the left product F*B where B is the 
 * n-clonal matrix B corresponding to T in O(mn) 
 * time where F is m x n.
 */
inline std::vector<std::vector<double>> left_multiply(
    const digraph<int>& clone_tree,
    int root,
    const std::unordered_map<int, int>& vertex_map,
    const std::vector<std::vector<double>>& frequency_matrix
) {
    // root is in column coordinates
    auto postorder = postorder_dfs(clone_tree, vertex_map.at(root));

    std::vector<std::vector<double>> result = frequency_matrix;
    for (auto i : postorder) {
        int k = clone_tree[i].data;
        for (auto child : clone_tree.successors(i)) {
            size_t l = clone_tree[child].data;
            for (size_t j = 0; j < frequency_matrix.size(); j++) {
                result[j][k] += result[j][l];
            }
        }
    }

    return result;
}

#endif // CLONETREE_H
