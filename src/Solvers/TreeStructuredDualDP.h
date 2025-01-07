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
    const std::vector<int>& post_order, 
    const std::vector<int>& vertex_map, 
    const std::vector<std::vector<float>>& frequency_matrix, 
    const std::vector<std::vector<float>>& weight_matrix, 
    int root,
    std::vector<Representation>& fs,
    std::vector<Representation>& gs,
    std::vector<float>& cs_buffer,
    size_t j
) {
    size_t ncols = frequency_matrix[0].size();

    for (auto u : post_order) {
        int i = clone_tree[u].data;

        // If leaf, compute Representation and return.
        if (clone_tree.out_degree(u) == 0) {
            fs[vertex_map[i]].reset_to_leaf(frequency_matrix[j][i], weight_matrix[j][i]);
            continue;
        }
        
        // Recurse at children. 
        const auto& children = clone_tree.successors(u);
        auto &g = gs[u];
        g.sum(children, fs);
        g.update_representation(frequency_matrix[j][i], weight_matrix[j][i], fs[vertex_map[i]], cs_buffer);
    }
}

#endif
