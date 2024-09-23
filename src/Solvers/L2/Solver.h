#ifndef QUADRATIC_SOLVER_H
#define QUADRATIC_SOLVER_H

#include <unordered_map>
#include <vector>
#include <algorithm>

#include "../../DiGraph.h"
#include "PiecewiseQuadraticF.h"

namespace L2Solver {
    class Solver {
        private:
            digraph<int> clone_tree;
            std::unordered_map<int, int> vertex_map;
            std::vector<std::vector<double>> frequency_matrix;
            std::vector<std::vector<double>> weight_matrix;
            int root;
        
            std::unordered_map<int, std::vector<PiecewiseQuadraticF>> fs;
            std::vector<std::vector<double>> alphas;

            /* Backtracks to compute the dual 
             * variables for sample j */
            void backtrack(double alpha_0, int j);
        public:
            double objective;
            std::vector<std::vector<double>> frequencies;

            Solver(
                digraph<int> clone_tree, 
                std::unordered_map<int, int> vertex_map, 
                std::vector<std::vector<double>> frequency_matrix, 
                std::vector<std::vector<double>> weight_matrix,
                int root
            ) : clone_tree(clone_tree), vertex_map(vertex_map), root(root) {
                this->frequency_matrix.reserve(frequency_matrix.size());
                for (size_t i = 0; i < frequency_matrix.size(); i++) {
                    this->frequency_matrix.push_back(std::vector<double>(frequency_matrix[i].begin(), frequency_matrix[i].end()));
                }

                for (size_t i = 0; i < frequency_matrix.size(); i++) {
                    this->alphas.push_back(std::vector<double>(frequency_matrix[i].size()));
                    this->frequencies.push_back(std::vector<double>(frequency_matrix[i].size()));
                }

                this->weight_matrix.reserve(weight_matrix.size());
                for (size_t i = 0; i < weight_matrix.size(); i++) {
                    this->weight_matrix.push_back(std::vector<double>(weight_matrix[i].begin(), weight_matrix[i].end()));
                }
            }

            void solve();
    };
};

#endif
