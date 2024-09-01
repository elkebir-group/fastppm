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
            std::vector<std::vector<double>> variant_reads;
            std::vector<std::vector<double>> total_reads;
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
                std::vector<std::vector<int>> variant_reads, 
                std::vector<std::vector<int>> total_reads, 
                int root
            ) : clone_tree(clone_tree), vertex_map(vertex_map), root(root) {
                this->variant_reads.reserve(variant_reads.size());
                for (size_t i = 0; i < variant_reads.size(); i++) {
                    this->variant_reads.push_back(std::vector<double>(variant_reads[i].begin(), variant_reads[i].end()));
                }

                this->total_reads.reserve(total_reads.size());
                for (size_t i = 0; i < total_reads.size(); i++) {
                    this->total_reads.push_back(std::vector<double>(total_reads[i].begin(), total_reads[i].end()));
                }

                for (size_t i = 0; i < total_reads.size(); i++) {
                    this->alphas.push_back(std::vector<double>(total_reads[i].size()));
                    this->frequencies.push_back(std::vector<double>(total_reads[i].size()));
                }
            }

            void solve();
    };
};

#endif
