#ifndef MAIN_H
#define MAIN_H

#include <vector>
#include <optional>
#include <unordered_map>

#include "DiGraph.h"

/* 
 * Represents the solver output.
 * 
 * Fields:
 *  - runtime: The time taken to solve the optimization problem in milliseconds.
 *  - objective: The value of the objective function after optimization.
 *  - usage_matrix: The usage matrix after optimization. This is an optional field.
 *  - frequency_matrix: The frequency matrix after optimization. This is an optional field.
 */
struct SolverResult {
    double runtime;
    double objective;
    std::optional<std::vector<std::vector<double>>> usage_matrix;
    std::optional<std::vector<std::vector<double>>> frequency_matrix;
};

SolverResult l2_solve(
    const std::unordered_map<int, int>& vertex_map,
    const std::vector<std::vector<int>>& variant_matrix,
    const std::vector<std::vector<int>>& total_matrix,
    const std::vector<std::vector<double>>& weight_matrix,
    const digraph<int>& clone_tree,
    size_t root
);

SolverResult log_binomial_solve(
    const std::unordered_map<int, int>& vertex_map,
    const std::vector<std::vector<int>>& variant_matrix,
    const std::vector<std::vector<int>>& total_matrix,
    const digraph<int>& clone_tree,
    size_t root,
    int K
);

SolverResult log_binomial_fixed_solve(
        const std::unordered_map<int, int>& vertex_map,
        const std::vector<std::vector<int>>& variant_matrix,
        const std::vector<std::vector<int>>& total_matrix,
        const digraph<int>& clone_tree,
        size_t root,
        int K
);

#endif
