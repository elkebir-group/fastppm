#include <vector>
#include <fstream>
#include <iostream>
#include <chrono>

#include "Lib.h"
#include "Polynomial.h"
#include "CloneTree.h"
#include "DiGraph.h"
#include "Solvers/PiecewiseLinear/Solver.h"
#include "Solvers/L2/Solver.h"
#include "Solvers/SeparableADMM/Solver.h"

/* 
 * Solves the optimization problem using the L2 loss function.
 */
SolverResult l2_solve(
    const std::vector<int>& vertex_map,
    const std::vector<std::vector<int>>& variant_matrix,
    const std::vector<std::vector<int>>& total_matrix,
    const std::vector<std::vector<float>>& weight_matrix,
    const digraph<int>& clone_tree,
    size_t root
) {
    std::vector<std::vector<float>> frequency_matrix;
    for (size_t i = 0; i < variant_matrix.size(); i++) {
        std::vector<float> frequencies;
        for (size_t j = 0; j < variant_matrix[i].size(); j++) {
            float freq = total_matrix[i][j] == 0 ? 0 : static_cast<float>(variant_matrix[i][j]) / total_matrix[i][j];
            frequencies.push_back(freq);
        }
        frequency_matrix.push_back(frequencies);
    }

    L2Solver::Solver solver(clone_tree, vertex_map, frequency_matrix, weight_matrix, root);

    auto start = std::chrono::high_resolution_clock::now();
    solver.initialize();
    solver.solve();
    auto end = std::chrono::high_resolution_clock::now();
    double runtime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    auto usage_matrix = left_inverse(clone_tree, vertex_map, solver.frequencies);
    return {runtime, solver.objective, usage_matrix, solver.frequencies};
}

SolverResult log_binomial_admm_solve(
    const std::vector<int>& vertex_map,
    const std::vector<std::vector<int>>& variant_matrix,
    const std::vector<std::vector<int>>& total_matrix,
    const digraph<int>& clone_tree,
    size_t root
) {
    /* 
     * We pass in the log-binomial loss:
     *          L_i(f) = -(v_i log f + (d_i - v_i) log (1 - f)) 
     * to the SeparableADMM solver along with functions for 
     * computing the gradient and Hessian.
     */

    const std::function<double(double,int,int)> compute_obj = [](double freq, int var, int tot) {
        if (var == 0) {
            return -(tot - var) * log(1 - freq) ;
        } else if (var == tot) {
            return -var * log(freq);
        } else {
            return -var * log(freq) - (tot - var) * log(1 - freq);
        }
    };

    const std::function<double(double,double,int,int)> compute_minimizer = [](double rho, double w, int var, int tot) {
        if (var == 0) { 
            double a = -rho, b = rho + rho * w, c = tot - var - rho * w;
            double r1 = -1.0f, r2 = -1.0f;
            solve_quadratic(a, b, c, r1, r2);
            if (r1 >= 0 && r1 <= 1) {
                return r1;
            } else if (r2 >= 0 && r2 <= 1) {
                return r2;
            } 

            return 0.0;
        } 
        
        if (var == tot) {
            double a = rho, b = -rho * w, c = -var;
            double r1 = -1.0f, r2 = -1.0f;
            solve_quadratic(a, b, c, r1, r2);
            if (r1 >= 0 && r1 <= 1) {
                return r1;
            } else if (r2 >= 0 && r2 <= 1) {
                return r2;
            } 

            return 1.0;
        }

        // since var != 0 && var != tot, we need to solve a cubic eq over open
        // interval (0, 1)
        double a = -rho, b = rho + rho * w, c = tot - rho * w, d = -var;
        double r1 = -1.0f, r2 = -1.0f, r3 = -1.0f;
        solve_cubic(a, b, c, d, r1, r2, r3);
        if (r1 >= 0 && r1 <= 1) {
            return r1;
        } else if (r2 >= 0 && r2 <= 1) {
            return r2;
        } else if (r3 >= 0 && r3 <= 1) {
            return r3;
        }

        return -1.0;
    };


    SeparableADMM::Solver solver(
            compute_obj, 
            compute_minimizer, 
            clone_tree, 
            vertex_map, 
            variant_matrix, 
            total_matrix, 
            root, 
            20, 100, 1e-6 // ADMM parameters
    );

    auto start = std::chrono::high_resolution_clock::now();
    solver.solve();
    auto end = std::chrono::high_resolution_clock::now();
    double runtime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    return {runtime, solver.objective, solver.usages, solver.frequencies};
}

SolverResult log_binomial_fixed_solve(
    const std::vector<int>& vertex_map,
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

    std::vector<std::vector<float>> frequency_matrix;

    auto start = std::chrono::high_resolution_clock::now();
    double objective = 0;
    for (size_t i = 0; i < variant_matrix.size(); i++) {
        std::vector<int> ref_vector(variant_matrix[i].size(), 0);
        for (size_t j = 0; j < variant_matrix[i].size(); j++) {
            ref_vector[j] = total_matrix[i][j] - variant_matrix[i][j];
        }

        PiecewiseLinearSolver::SolverBinomial solver(K);
        solver.init(variant_matrix[i], ref_vector, link_list, root);
        objective += solver.solve(1e-4);

        std::vector<float> frequencies(variant_matrix[i].size(), 0);
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
    const std::vector<int>& vertex_map,
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

    std::vector<std::vector<float>> frequency_matrix;

    auto start = std::chrono::high_resolution_clock::now();
    float objective = 0;
    for (size_t i = 0; i < variant_matrix.size(); i++) {
        std::vector<int> ref_vector(variant_matrix[i].size(), 0);
        for (size_t j = 0; j < variant_matrix[i].size(); j++) {
            ref_vector[j] = total_matrix[i][j] - variant_matrix[i][j];
        }

        PiecewiseLinearSolver::SolverBinomial solver(K);
        solver.init(variant_matrix[i], ref_vector, link_list, root);
        objective += solver.solve_iteratively(0.75, 1e-4); // TODO: make these parameters configurable
        
        std::vector<float> frequencies(variant_matrix[i].size(), 0);
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

SolverResult log_beta_binomial_fixed_solve(
        const std::vector<int>& vertex_map,
        const std::vector<std::vector<int>>& variant_matrix,
        const std::vector<std::vector<int>>& total_matrix,
        const digraph<int>& clone_tree,
        size_t root,
        int K,
        double s
) {
    auto edges = clone_tree.edges();
    int n_clones = edges.size() + 1;
    std::vector<std::list<int>> link_list(n_clones);
    for (auto& [u, v] : edges) {
        link_list[clone_tree[u].data].push_back(clone_tree[v].data);
    }

    std::vector<std::vector<float>> frequency_matrix;

    auto start = std::chrono::high_resolution_clock::now();
    double objective = 0;
    for (size_t i = 0; i < variant_matrix.size(); i++) {
        std::vector<int> ref_vector(variant_matrix[i].size(), 0);
        for (size_t j = 0; j < variant_matrix[i].size(); j++) {
            ref_vector[j] = total_matrix[i][j] - variant_matrix[i][j];
        }

        PiecewiseLinearSolver::SolverBetaBinomial solver(K, s);
        solver.init(variant_matrix[i], ref_vector, link_list, root);
        objective += solver.solve(1e-4);

        std::vector<float> frequencies(variant_matrix[i].size(), 0);
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

SolverResult log_beta_binomial_solve(
        const std::vector<int>& vertex_map,
        const std::vector<std::vector<int>>& variant_matrix,
        const std::vector<std::vector<int>>& total_matrix,
        const digraph<int>& clone_tree,
        size_t root,
        int K,
        double s
) {
    auto edges = clone_tree.edges();
    int n_clones = edges.size() + 1;
    std::vector<std::list<int>> link_list(n_clones);
    for (auto& [u, v] : edges) {
        link_list[clone_tree[u].data].push_back(clone_tree[v].data);
    }

    std::vector<std::vector<float>> frequency_matrix;

    auto start = std::chrono::high_resolution_clock::now();
    float objective = 0;
    for (size_t i = 0; i < variant_matrix.size(); i++) {
        std::vector<int> ref_vector(variant_matrix[i].size(), 0);
        for (size_t j = 0; j < variant_matrix[i].size(); j++) {
            ref_vector[j] = total_matrix[i][j] - variant_matrix[i][j];
        }

        PiecewiseLinearSolver::SolverBetaBinomial solver(K, s);
        solver.init(variant_matrix[i], ref_vector, link_list, root);
        objective += solver.solve_iteratively(0.75, 1e-4, 1e-4); // TODO: make these parameters configurable

        std::vector<float> frequencies(variant_matrix[i].size(), 0);
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
