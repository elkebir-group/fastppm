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
            std::vector<std::vector<float>> variant_reads;
            std::vector<std::vector<float>> total_reads;
            int root;
        
            std::unordered_map<int, std::vector<PiecewiseQuadraticF>> fs;
        public:
            float objective;

            Solver(
                digraph<int> clone_tree, 
                std::unordered_map<int, int> vertex_map, 
                std::vector<std::vector<int>> variant_reads, 
                std::vector<std::vector<int>> total_reads, 
                int root
            ) : clone_tree(clone_tree), vertex_map(vertex_map), root(root) {
                this->variant_reads.reserve(variant_reads.size());
                for (size_t i = 0; i < variant_reads.size(); i++) {
                    this->variant_reads.push_back(std::vector<float>(variant_reads[i].begin(), variant_reads[i].end()));
                }

                this->total_reads.reserve(total_reads.size());
                for (size_t i = 0; i < total_reads.size(); i++) {
                    this->total_reads.push_back(std::vector<float>(total_reads[i].begin(), total_reads[i].end()));
                }
            }

            void solve();
    };
};

#endif
