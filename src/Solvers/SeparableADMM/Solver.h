#ifndef SEPARABLE_ADMM_SOLVER_H
#define SEPARABLE_ADMM_SOLVER_H

#include <unordered_map>
#include <vector>
#include <algorithm>
#include <functional>

#include "DiGraph.h"

namespace SeparableADMM {
    class Solver {
        private:
            digraph<int> clone_tree;
            std::vector<int> vertex_map;
            std::vector<std::vector<int>> variant_reads;
            std::vector<std::vector<int>> total_reads;
            int root;

            int num_admm_iterations;
            double rho;
            float frequency_clamp;
            bool verbose;

            std::function<double(double,int,int)> compute_obj;
            std::function<double(double,double,int,int)> compute_minimizer;

            L2Solver::Solver l2_solver;

            void frequency_update();
            void usage_update();
            void residual_update();
            void objective_update();
        public:
            double objective;

            std::vector<std::vector<float>> frequencies;
            std::vector<std::vector<float>> usages;
            std::vector<std::vector<float>> residuals;
            std::vector<std::vector<float>> buffer;

            Solver(
                std::function<double(double,int,int)> compute_obj,
                std::function<double(double,double,int,int)> compute_minimizer,
                digraph<int> clone_tree, 
                std::vector<int> vertex_map, 
                std::vector<std::vector<int>> variant_reads,
                std::vector<std::vector<int>> total_reads,
                int root,
                int num_admm_iterations,
                float rho,
                float frequency_clamp,
                bool verbose = false
            ) : clone_tree(clone_tree), 
                vertex_map(vertex_map), 
                root(root), 
                num_admm_iterations(num_admm_iterations),
                rho(rho),
                frequency_clamp(frequency_clamp),
                verbose(verbose),
                compute_obj(compute_obj),
                compute_minimizer(compute_minimizer)
            {
                this->variant_reads.reserve(variant_reads.size());
                for (size_t i = 0; i < variant_reads.size(); i++) {
                    this->variant_reads.push_back(std::vector<int>(variant_reads[i].begin(), variant_reads[i].end()));
                }

                for (size_t i = 0; i < variant_reads.size(); i++) {
                    this->total_reads.push_back(std::vector<int>(total_reads[i].begin(), total_reads[i].end()));
                }

                for (size_t i = 0; i < variant_reads.size(); i++) {
                    this->usages.push_back(std::vector<float>(variant_reads[i].size(), 0.0f));
                    this->frequencies.push_back(std::vector<float>(variant_reads[i].size(), 0.0f));
                    this->residuals.push_back(std::vector<float>(variant_reads[i].size(), 0.0f));
                    this->buffer.push_back(std::vector<float>(variant_reads[i].size(), 0.0f));
                }
            }

            void solve();
    };
};

#endif
