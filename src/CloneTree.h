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
inline void left_inverse(
    const digraph<int>& clone_tree,
    const std::vector<int>& vertex_map,
    const std::vector<std::vector<float>>& frequency_matrix,
    std::vector<std::vector<float>>& result
) {
    for (size_t i = 0; i < result.size(); i++) {
        for (size_t j = 0; j < result[i].size(); j++) {
            result[i][j] = frequency_matrix[i][j];
        }
    }

    for (size_t j = 0; j < frequency_matrix.size(); j++) {
        for (size_t i = 0; i < frequency_matrix[j].size(); i++) {
            for (auto child : clone_tree.successors(vertex_map[i])) {
                size_t l = clone_tree[child].data;
                result[j][i] -= frequency_matrix[j][l];
            }
        }
    }
}

inline std::vector<std::vector<float>> left_inverse(
    const digraph<int>& clone_tree,
    const std::vector<int>& vertex_map,
    const std::vector<std::vector<float>>& frequency_matrix
) {
    std::vector<std::vector<float>> result(frequency_matrix.size(), std::vector<float>(frequency_matrix[0].size(), 0.0f));
    left_inverse(clone_tree, vertex_map, frequency_matrix, result);
    return result;
}

/* 
 * Computes the left product F*B where B is the 
 * n-clonal matrix B corresponding to T in O(mn) 
 * time where F is m x n.
 */
inline void left_multiply(
    const digraph<int>& clone_tree,
    int root,
    const std::vector<int>& vertex_map,
    const std::vector<std::vector<float>>& frequency_matrix,
    std::vector<std::vector<float>>& result
) {
    // root is in column coordinates
    auto postorder = postorder_dfs(clone_tree, vertex_map[root]);

    for (size_t i = 0; i < result.size(); i++) {
        for (size_t j = 0; j < result[i].size(); j++) {
            result[i][j] = frequency_matrix[i][j];
        }
    }

    for (auto i : postorder) {
        int k = clone_tree[i].data;
        for (auto child : clone_tree.successors(i)) {
            size_t l = clone_tree[child].data;
            for (size_t j = 0; j < frequency_matrix.size(); j++) {
                result[j][k] += result[j][l];
            }
        }
    }
}

inline std::vector<std::vector<float>> left_multiply(
    const digraph<int>& clone_tree,
    int root,
    const std::vector<int>& vertex_map,
    const std::vector<std::vector<float>>& frequency_matrix
) {
    std::vector<std::vector<float>> result(frequency_matrix.size(), std::vector<float>(frequency_matrix[0].size(), 0.0f));
    left_multiply(clone_tree, root, vertex_map, frequency_matrix, result);
    return result;
}

#endif // CLONETREE_H
