#ifndef TREE_STRUCTURED_DUAL_DP_H
#define TREE_STRUCTURED_DUAL_DP_H

#include <vector>
#include <unordered_map>
#include <stack>

#include "../DiGraph.h"

/* 
 * Solves the dual problem over a clone tree using the
 * abstract tree structured dual dynamic programming algorithm.
 * 
 * Requires the following methods to be defined for Representation:
 * - Representation::Representation(double variant_reads, double total_reads) 
 *      Computes J_i(.) for a leaf vertex.
 * - Representation Representation::operator+(const Representation& other)
 *      Adds J_i(.) and J_k(.) for two children vertices.
 * - Representation Representation::update_representation(double variant_reads, double total_reads)
 *      Computes J_i(.) <- min_{h_i(var, total) + J_i(.)} for a vertex with variant_reads 
 *      and total_reads.
 */
template <typename Representation>
std::unordered_map<int, std::vector<Representation>> forward_solve(
    digraph<int>& clone_tree, 
    const std::unordered_map<int, int>& vertex_map, 
    const std::vector<std::vector<double>>& variant_reads, 
    const std::vector<std::vector<double>>& total_reads, 
    int root
) {
    size_t nrows = variant_reads.size();
    size_t ncols = variant_reads[0].size();

    // stores the representations for each vertex and sample
    std::unordered_map<int, std::vector<Representation>> fs;
    for (auto& p : vertex_map) {
        fs[p.first] = std::vector<Representation>(nrows);
    }

    std::vector<bool> visited(ncols, false);

    std::stack<int> call_stack; // contains vertices to visit in column coordinates
    call_stack.push(root);
    while(!call_stack.empty()) {
        int i = call_stack.top();
        call_stack.pop();

        if (visited[vertex_map.at(i)]) continue;

        // If leaf, compute Representation and return.
        if (clone_tree.out_degree(vertex_map.at(i)) == 0) {
            for (size_t j = 0; j < nrows; ++j) {
                fs[vertex_map.at(i)][j] = Representation(variant_reads[j][i], total_reads[j][i]);
            }

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

        for (size_t j = 0; j < nrows; ++j) {
            Representation g;
            for (auto k : clone_tree.successors(vertex_map.at(i))) {
                Representation f = fs[k][j];
                g = g + f; 
            }

            fs[vertex_map.at(i)][j] = g.update_representation(variant_reads[j][i], total_reads[j][i]);
        }

        visited[vertex_map.at(i)] = true;
    }

    return fs;
}

#endif
