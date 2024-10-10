#include <vector>
#include <fstream>
#include <iostream>
#include <chrono>

#include "Lib.h"
#include "CloneTree.h"
#include "DiGraph.h"
#include "Solvers/LogBinomialPiecewise/Solver.h"
#include "Solvers/L2/Solver.h"

/* 
 * Solves the optimization problem using the L2 loss function.
 */
SolverResult l2_solve(
    const std::unordered_map<int, int>& vertex_map,
    const std::vector<std::vector<int>>& variant_matrix,
    const std::vector<std::vector<int>>& total_matrix,
    const std::vector<std::vector<double>>& weight_matrix,
    const digraph<int>& clone_tree,
    size_t root
) {
    std::vector<std::vector<double>> frequency_matrix;
    for (size_t i = 0; i < variant_matrix.size(); i++) {
        std::vector<double> frequencies;
        for (size_t j = 0; j < variant_matrix[i].size(); j++) {
            double freq = total_matrix[i][j] == 0 ? 0 : static_cast<double>(variant_matrix[i][j]) / total_matrix[i][j];
            frequencies.push_back(freq);
        }
        frequency_matrix.push_back(frequencies);
    }

    L2Solver::Solver solver(clone_tree, vertex_map, frequency_matrix, weight_matrix, root);

    auto start = std::chrono::high_resolution_clock::now();
    solver.solve();
    auto end = std::chrono::high_resolution_clock::now();
    double runtime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    auto usage_matrix = left_inverse(clone_tree, vertex_map, solver.frequencies);
    return {runtime, solver.objective, usage_matrix, solver.frequencies};
}

SolverResult log_binomial_fixed_solve(
    const std::unordered_map<int, int>& vertex_map,
    const std::vector<std::vector<int>>& variant_matrix,
    const std::vector<std::vector<int>>& total_matrix,
    const digraph<int>& clone_tree,
    size_t root,
    int K
) {
    auto edges = clone_tree.edges();
    int n_clones = edges.size() + 1;
    std::vector<std::list<int>> link_list(n_clones);
    for (auto& [u, v] : edges) {
        link_list[clone_tree[u].data].push_back(clone_tree[v].data);
    }

    std::vector<std::vector<double>> frequency_matrix;

    auto start = std::chrono::high_resolution_clock::now();
    double objective = 0;
    for (size_t i = 0; i < variant_matrix.size(); i++) {
        std::vector<int> ref_vector(variant_matrix[i].size(), 0);
        for (size_t j = 0; j < variant_matrix[i].size(); j++) {
            ref_vector[j] = total_matrix[i][j] - variant_matrix[i][j];
        }

        LogBinomialPiecewiseLinearSolver::Solver solver(K);
        solver.init(variant_matrix[i], ref_vector, link_list, root);
        objective += solver.solve(1e-4);

        std::vector<double> frequencies(variant_matrix[i].size(), 0);
        for (size_t j = 0; j < variant_matrix[i].size(); j++) {
            frequencies[j] = solver.F[j];
        }

        frequency_matrix.push_back(frequencies);
    }
    auto end = std::chrono::high_resolution_clock::now();
    double runtime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    auto usage_matrix = left_inverse(clone_tree, vertex_map, frequency_matrix);
    return {runtime, objective, usage_matrix, frequency_matrix};
}
/* 
 * Solves the optimization problem using the binomial loss function
 * with Yuanyuan's progressive piecewise linear solver.
 */
SolverResult log_binomial_solve(
    const std::unordered_map<int, int>& vertex_map,
    const std::vector<std::vector<int>>& variant_matrix,
    const std::vector<std::vector<int>>& total_matrix,
    const digraph<int>& clone_tree,
    size_t root,
    int K
) {
    auto edges = clone_tree.edges();
    int n_clones = edges.size() + 1;
    std::vector<std::list<int>> link_list(n_clones);
    for (auto& [u, v] : edges) {
        link_list[clone_tree[u].data].push_back(clone_tree[v].data);
    }

    std::vector<std::vector<double>> frequency_matrix;

    auto start = std::chrono::high_resolution_clock::now();
    double objective = 0;
    for (size_t i = 0; i < variant_matrix.size(); i++) {
        std::vector<int> ref_vector(variant_matrix[i].size(), 0);
        for (size_t j = 0; j < variant_matrix[i].size(); j++) {
            ref_vector[j] = total_matrix[i][j] - variant_matrix[i][j];
        }

        LogBinomialPiecewiseLinearSolver::Solver solver(K);
        solver.init(variant_matrix[i], ref_vector, link_list, root);
        objective += solver.solve_iteratively(0.75, 1e-4); // TODO: make these parameters configurable
        
        std::vector<double> frequencies(variant_matrix[i].size(), 0);
        for (size_t j = 0; j < variant_matrix[i].size(); j++) {
            frequencies[j] = solver.F[j];
        }

        frequency_matrix.push_back(frequencies);
    }
    auto end = std::chrono::high_resolution_clock::now();
    double runtime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    
    auto usage_matrix = left_inverse(clone_tree, vertex_map, frequency_matrix);
    return {runtime, objective, usage_matrix, frequency_matrix};
}

