#ifndef LOG_BINOMIAL_ADMM_SOLVER_H
#define LOG_BINOMIAL_ADMM_SOLVER_H

#include <unordered_map>
#include <vector>
#include <algorithm>

#include "../../DiGraph.h"

namespace LogBinomialADMM {
    class Solver {
        private:
            digraph<int> clone_tree;
            std::unordered_map<int, int> vertex_map;
            std::vector<std::vector<int>> variant_reads;
            std::vector<std::vector<int>> total_reads;
            int root;

            double rho = 100.0;

            void frequency_update();
            void usage_update();
            void residual_update();
            void objective_update();
            void lagrangian_objective_update();
        public:
            double objective;
            double lagrangian_objective;

            std::vector<std::vector<double>> frequencies;
            std::vector<std::vector<double>> usages;
            std::vector<std::vector<double>> residuals;

            Solver(
                digraph<int> clone_tree, 
                std::unordered_map<int, int> vertex_map, 
                std::vector<std::vector<int>> variant_reads,
                std::vector<std::vector<int>> total_reads,
                int root
            ) : clone_tree(clone_tree), vertex_map(vertex_map), root(root) {
                this->variant_reads.reserve(variant_reads.size());
                for (size_t i = 0; i < variant_reads.size(); i++) {
                    this->variant_reads.push_back(std::vector<int>(variant_reads[i].begin(), variant_reads[i].end()));
                }

                for (size_t i = 0; i < variant_reads.size(); i++) {
                    this->total_reads.push_back(std::vector<int>(total_reads[i].begin(), total_reads[i].end()));
                }

                for (size_t i = 0; i < variant_reads.size(); i++) {
                    this->usages.push_back(std::vector<double>(variant_reads[i].size(), 0.0));
                    this->frequencies.push_back(std::vector<double>(variant_reads[i].size(), 0.0));
                    this->residuals.push_back(std::vector<double>(variant_reads[i].size(), 0.0));
                }
            }

            void solve();
    };
};

#endif
