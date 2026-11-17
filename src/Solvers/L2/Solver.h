#ifndef QUADRATIC_SOLVER_H
#define QUADRATIC_SOLVER_H

#include <unordered_map>
#include <vector>
#include <algorithm>

#include "DiGraph.h"
#include "PiecewiseQuadraticF.h"

namespace L2Solver {
    /*
     * Solves min_f sum_i w_i (f_i - fhat_i)^2 subject to the sum condition,
     * by the tree-structured dual dynamic program.
     *
     * By default the root frequency is at most 1 (the root price alpha_0 is
     * clamped at 0). With set_fix_root_frequency(true), the root frequency is
     * exactly 1 and alpha_0 is free.
     *
     * T is the scalar type; Solver = BasicSolver<float> is the original solver.
     *
     * Possible speed-ups, not measured (the solver costs about the sum of the
     * subtree sizes per row, so each is a constant factor):
     * - Skip the objective when the caller does not need it (fastppm-CN does
     *   not): f0 needs one evaluation of F per node in update_representation.
     * - Keep the derivative intercepts of each node from the forward pass,
     *   instead of recomputing them in compute_argmin during the backtrack.
     * - Write new targets into the existing buffer in update_frequency_matrix,
     *   instead of reallocating it on every call.
     * - Reuse the index buffer of PiecewiseQuadraticF::sum, and merge with a
     *   heap at nodes with many children.
     * - Find pieces by binary search (breakpoints are sorted).
     */
    template <typename T>
    class BasicSolver {
        private:
            digraph<int> clone_tree;
            std::vector<int> vertex_map;
            std::vector<std::vector<T>> frequency_matrix;
            std::vector<std::vector<T>> weight_matrix;
            int root;
            bool fix_root_frequency = false;

            std::vector<PiecewiseQuadraticF<T>> fs;
            std::vector<PiecewiseQuadraticF<T>> gs;
            std::vector<T> cs;
            std::vector<std::vector<T>> alphas;

            std::vector<int> postorder;
            std::vector<int> preorder;
            std::vector<int> num_descendants;

            /* Backtracks to compute the dual
             * variables for sample j */
            void backtrack(T alpha_0, int j);
        public:
            T objective;
            std::vector<std::vector<T>> frequencies;

            BasicSolver() {}
            BasicSolver(
                const digraph<int> &clone_tree,
                const std::vector<int> &vertex_map,
                const std::vector<std::vector<T>> &frequency_matrix,
                const std::vector<std::vector<T>> &weight_matrix,
                int root
            ) : clone_tree(clone_tree), vertex_map(vertex_map), root(root) {
                this->frequency_matrix.reserve(frequency_matrix.size());
                for (size_t i = 0; i < frequency_matrix.size(); i++) {
                    this->frequency_matrix.push_back(std::vector<T>(frequency_matrix[i].begin(), frequency_matrix[i].end()));
                }

                for (size_t i = 0; i < frequency_matrix.size(); i++) {
                    this->alphas.push_back(std::vector<T>(frequency_matrix[i].size()));
                    this->frequencies.push_back(std::vector<T>(frequency_matrix[i].size()));
                }

                this->weight_matrix.reserve(weight_matrix.size());
                for (size_t i = 0; i < weight_matrix.size(); i++) {
                    this->weight_matrix.push_back(std::vector<T>(weight_matrix[i].begin(), weight_matrix[i].end()));
                }
            }

            void initialize();
            void solve();
            void set_fix_root_frequency(bool fix) { fix_root_frequency = fix; }

            void update_frequency_matrix(const std::vector<std::vector<T>> &frequency_matrix) {
                this->frequency_matrix.clear();
                this->frequency_matrix.reserve(frequency_matrix.size());
                for (size_t i = 0; i < frequency_matrix.size(); i++) {
                    this->frequency_matrix.push_back(std::vector<T>(frequency_matrix[i].begin(), frequency_matrix[i].end()));
                }
            }

            void update_weight_matrix(const std::vector<std::vector<T>> &weight_matrix) {
                this->weight_matrix.clear();
                this->weight_matrix.reserve(weight_matrix.size());
                for (size_t i = 0; i < weight_matrix.size(); i++) {
                    this->weight_matrix.push_back(std::vector<T>(weight_matrix[i].begin(), weight_matrix[i].end()));
                }
            }
    };

    using Solver = BasicSolver<float>;
};

#endif
