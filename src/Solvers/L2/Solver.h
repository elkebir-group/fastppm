#ifndef QUADRATIC_SOLVER_H
#define QUADRATIC_SOLVER_H

#include <unordered_map>
#include <vector>
#include <algorithm>

#include "DiGraph.h"
#include "PiecewiseQuadraticF.h"

namespace L2Solver {
    class Solver {
        private:
            digraph<int> clone_tree;
            std::vector<int> vertex_map;
            std::vector<std::vector<float>> frequency_matrix;
            std::vector<std::vector<float>> weight_matrix;
            int root;
        
            std::vector<PiecewiseQuadraticF> fs;
            std::vector<PiecewiseQuadraticF> gs;
            std::vector<float> cs;
            std::vector<std::vector<float>> alphas;

            std::vector<int> postorder;
            std::vector<int> preorder;
            std::vector<int> num_descendants;

            /* Backtracks to compute the dual 
             * variables for sample j */
            void backtrack(float alpha_0, int j);
        public:
            float objective;
            std::vector<std::vector<float>> frequencies;

            Solver() {}
            Solver(
                const digraph<int> &clone_tree, 
                const std::vector<int> &vertex_map, 
                const std::vector<std::vector<float>> &frequency_matrix, 
                const std::vector<std::vector<float>> &weight_matrix,
                int root
            ) : clone_tree(clone_tree), vertex_map(vertex_map), root(root) {
                this->frequency_matrix.reserve(frequency_matrix.size());
                for (size_t i = 0; i < frequency_matrix.size(); i++) {
                    this->frequency_matrix.push_back(std::vector<float>(frequency_matrix[i].begin(), frequency_matrix[i].end()));
                }

                for (size_t i = 0; i < frequency_matrix.size(); i++) {
                    this->alphas.push_back(std::vector<float>(frequency_matrix[i].size()));
                    this->frequencies.push_back(std::vector<float>(frequency_matrix[i].size()));
                }

                this->weight_matrix.reserve(weight_matrix.size());
                for (size_t i = 0; i < weight_matrix.size(); i++) {
                    this->weight_matrix.push_back(std::vector<float>(weight_matrix[i].begin(), weight_matrix[i].end()));
                }
            }

            void initialize();
            void solve();
            void update_frequency_matrix(const std::vector<std::vector<float>> &frequency_matrix) {
                this->frequency_matrix.clear();
                this->frequency_matrix.reserve(frequency_matrix.size());
                for (size_t i = 0; i < frequency_matrix.size(); i++) {
                    this->frequency_matrix.push_back(std::vector<float>(frequency_matrix[i].begin(), frequency_matrix[i].end()));
                }
            }
    };
};

#endif
