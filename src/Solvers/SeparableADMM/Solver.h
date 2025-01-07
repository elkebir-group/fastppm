#ifndef SEPARABLE_ADMM_SOLVER_H
#define SEPARABLE_ADMM_SOLVER_H

#include <unordered_map>
#include <vector>
#include <algorithm>

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
            int max_newton_iterations;
            double newton_tolerance;
            float frequency_clamp;
            double rho = 1000.0;

            std::function<double(double,int,int)> compute_obj;
            std::function<double(double,int,int)> compute_gradient;
            std::function<double(double,int,int)> compute_hessian;

            L2Solver::Solver l2_solver;

            void frequency_update();
            void usage_update(const std::vector<std::vector<float>>& weight_matrix);
            void residual_update();
            void objective_update();
        public:
            double objective;

            std::vector<std::vector<float>> frequencies;
            std::vector<std::vector<float>> usages;
            std::vector<std::vector<float>> residuals;

            Solver(
                std::function<double(double,int,int)> compute_obj,
                std::function<double(double,int,int)> compute_gradient,
                std::function<double(double,int,int)> compute_hessian,
                digraph<int> clone_tree, 
                std::vector<int> vertex_map, 
                std::vector<std::vector<int>> variant_reads,
                std::vector<std::vector<int>> total_reads,
                int root,
                int num_admm_iterations,
                int max_newton_iterations,
                double newton_tolerance,
                float frequency_clamp
            ) : clone_tree(clone_tree), 
                vertex_map(vertex_map), 
                root(root), 
                num_admm_iterations(num_admm_iterations),
                max_newton_iterations(max_newton_iterations),
                newton_tolerance(newton_tolerance),
                frequency_clamp(frequency_clamp),
                compute_obj(compute_obj),
                compute_gradient(compute_gradient),
                compute_hessian(compute_hessian)
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
                }
            }

            void solve();
    };
};

#endif
