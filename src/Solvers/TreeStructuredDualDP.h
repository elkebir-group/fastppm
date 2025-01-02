#ifndef TREE_STRUCTURED_DUAL_DP_H
#define TREE_STRUCTURED_DUAL_DP_H

#include <vector>
#include <unordered_map>
#include <stack>

#include "DiGraph.h"

/* 
 * Solves the dual problem over a clone tree using the
 * abstract tree structured dual dynamic programming algorithm.
 * 
 * Requires the following methods to be defined for Representation:
 * - Representation::Representation(float frequency) 
 *      Computes J_i(.) for a leaf vertex.
 * - Representation Representation::operator+(const Representation& other)
 *      Adds J_i(.) and J_k(.) for two children vertices.
 * - Representation Representation::update_representation(float frequency)
 *      Computes J_i(.) <- min_{h_i(var, total) + J_i(.)} for 
 *      a vertex with a given frequency. 
 */
template <typename Representation>
void forward_solve(
    digraph<int>& clone_tree, 
    const std::vector<int>& vertex_map, 
    const std::vector<std::vector<float>>& frequency_matrix, 
    const std::vector<std::vector<float>>& weight_matrix, 
    int root,
    std::vector<Representation>& fs,
    std::vector<Representation>& gs,
    size_t j
) {
    size_t ncols = frequency_matrix[0].size();

    std::vector<bool> visited(ncols, false);

    std::stack<int> call_stack; // contains vertices to visit in column coordinates
    call_stack.push(root);
    while(!call_stack.empty()) {
        int i = call_stack.top();
        call_stack.pop();

        if (visited[vertex_map[i]]) continue;

        // If leaf, compute Representation and return.
        if (clone_tree.out_degree(vertex_map[i]) == 0) {
            fs[vertex_map[i]].reset_to_leaf(frequency_matrix[j][i], weight_matrix[j][i]);
            visited[vertex_map[i]] = true;
            continue;
        }
        
        // Recurse at children. 
        const auto& children = clone_tree.successors(vertex_map[i]);
        bool all_children_valid = true; 
        for (auto k : children) {
            if (!visited[k]) {
                if (all_children_valid) {
                    call_stack.push(i);
                }

                call_stack.push(clone_tree[k].data);
                all_children_valid = false;
            }
        }

        if (!all_children_valid) continue;

        auto &g = gs[vertex_map[i]];
        g.sum(children, fs);
        g.update_representation(frequency_matrix[j][i], weight_matrix[j][i], fs[vertex_map[i]]);
        visited[vertex_map[i]] = true;
    }
}

#endif
