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
    const std::unordered_map<int, int>& vertex_map, 
    const std::vector<std::vector<float>>& frequency_matrix, 
    const std::vector<std::vector<float>>& weight_matrix, 
    int root,
    std::unordered_map<int, Representation>& fs,
    std::unordered_map<int, Representation>& gs,
    size_t j
) {
    size_t ncols = frequency_matrix[0].size();

    std::vector<bool> visited(ncols, false);

    std::stack<int> call_stack; // contains vertices to visit in column coordinates
    call_stack.push(root);
    while(!call_stack.empty()) {
        int i = call_stack.top();
        call_stack.pop();

        if (visited[vertex_map.at(i)]) continue;

        // If leaf, compute Representation and return.
        if (clone_tree.out_degree(vertex_map.at(i)) == 0) {
            fs[vertex_map.at(i)] = std::move(Representation(frequency_matrix[j][i], weight_matrix[j][i]));
            visited[vertex_map.at(i)] = true;
            continue;
        }

        // Recurse at children. 
        bool all_children_valid = true; 
        for (auto k : clone_tree.successors(vertex_map.at(i))) {
            if (!visited[k]) {
                if (all_children_valid) {
                    call_stack.push(i);
                }

                call_stack.push(clone_tree[k].data);
                all_children_valid = false;
            }
        }

        if (!all_children_valid) continue;

        Representation g;
        for (auto k : clone_tree.successors(vertex_map.at(i))) {
            const Representation& f = fs[k];
            g = std::move(g + f); 
        }

        fs[vertex_map.at(i)] = std::move(g.update_representation(frequency_matrix[j][i], weight_matrix[j][i]));
        gs[vertex_map.at(i)] = std::move(g);

        visited[vertex_map.at(i)] = true;
    }
}

#endif
