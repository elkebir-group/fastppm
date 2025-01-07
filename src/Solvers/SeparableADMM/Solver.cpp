#include <spdlog/spdlog.h>

#include "../../CloneTree.h"
#include "../L2/Solver.h"
#include "Solver.h" 

namespace SeparableADMM {
    void Solver::solve() {
        std::vector<std::vector<float>> frequency_matrix;
        for (size_t i = 0; i < variant_reads.size(); i++) {
            std::vector<float> frequencies;
            for (size_t j = 0; j < variant_reads[i].size(); j++) {
                float freq = total_reads[i][j] == 0 ? 0 : static_cast<float>(variant_reads[i][j]) / total_reads[i][j];
                frequencies.push_back(freq);
            }
            frequency_matrix.push_back(frequencies);
        }

        std::vector<std::vector<float>> weight_matrix = frequency_matrix;
        for (size_t i = 0; i < weight_matrix.size(); i++) {
            for (size_t j = 0; j < weight_matrix[i].size(); j++) {
                weight_matrix[i][j] = 1.0f;
            }
        }

        // Step 1. Find decent initial solution
        l2_solver = L2Solver::Solver(clone_tree, vertex_map, frequency_matrix, weight_matrix, root);
        l2_solver.initialize();
        l2_solver.solve();

        frequencies = l2_solver.frequencies;
        usages = left_inverse(clone_tree, vertex_map, frequencies);
        objective_update();

        for (int i = 0; i < num_admm_iterations; i++) {
            std::cout << "Objective: " << objective << std::endl;

            frequency_update();          // ADMM Step 1
            usage_update(weight_matrix); // ADMM Step 2
            residual_update();           // ADMM Step 3

            objective_update(); 
        }

        frequencies = left_multiply(clone_tree, root, vertex_map, usages);
    }

    void Solver::frequency_update() {
        auto ws = left_multiply(clone_tree, root, vertex_map, usages);
        for (size_t i = 0; i < ws.size(); i++) {
            for (size_t j = 0; j < ws[i].size(); j++) {
                ws[i][j] = +ws[i][j] + (residuals[i][j] / rho);
            }
        }
        /* 
         * Solve 1D optimization problems using 
         * damped Newton's with backtracking line search. 
         * TODO: Investigate further.
         */
        for (size_t i = 0; i < frequencies.size(); i++) {
            for (size_t j = 0; j < frequencies[i].size(); j++) {
                double freq = 0.5;
                double current_obj = compute_obj(freq, variant_reads[i][j], total_reads[i][j]) + rho * (freq - ws[i][j]) * (freq - ws[i][j]);
                for (int k = 0; k < max_newton_iterations; k++) {
                    double current_grad = compute_gradient(freq, variant_reads[i][j], total_reads[i][j]) + rho * (freq-ws[i][j]);
                    double current_hess = compute_hessian(freq, variant_reads[i][j], total_reads[i][j]) + rho;

                    double step_size = 1.0;
                    double freq_p = freq - current_grad * step_size / current_hess;
                    double updated_obj = compute_obj(freq_p, variant_reads[i][j], total_reads[i][j]) + rho * (freq_p - ws[i][j]) * (freq_p - ws[i][j]);
                    while (updated_obj > current_obj) {
                        step_size *= 0.90;
                        freq_p = freq - current_grad * step_size / current_hess;
                        updated_obj = compute_obj(freq_p, variant_reads[i][j], total_reads[i][j]) + rho * (freq_p - ws[i][j]) * (freq_p - ws[i][j]);
                    }

                    freq = freq_p;

                    if (std::abs(current_obj - updated_obj) < newton_tolerance 
                            || std::abs(current_grad) < newton_tolerance
                            || freq > 1 - frequency_clamp 
                            || freq < frequency_clamp) {
                        break;
                    }
                }

                frequencies[i][j] = freq;
            }
        }
    }

    void Solver::usage_update(const std::vector<std::vector<float>>& weight_matrix) {
        std::vector<std::vector<float>> frequency_matrix = frequencies;
        for (size_t i = 0; i < frequency_matrix.size(); i++) {
            for (size_t j = 0; j < frequency_matrix[i].size(); j++) {
                frequency_matrix[i][j] = frequencies[i][j] - (residuals[i][j] / rho);
            }
        }
        l2_solver.update_frequency_matrix(frequency_matrix);
        l2_solver.solve();
        usages = left_inverse(clone_tree, vertex_map, l2_solver.frequencies);
    }

    void Solver::residual_update() {
        auto ws = left_multiply(clone_tree, root, vertex_map, usages);
            
        for (size_t i = 0; i < ws.size(); i++) {
            for (size_t j = 0; j < ws[i].size(); j++) {
                residuals[i][j] += rho * (ws[i][j] - frequencies[i][j]);
            }
        }
    }

    void Solver::objective_update() {
        auto frequency_matrix = left_multiply(clone_tree, root, vertex_map, usages);
        float obj = 0.0;
        for (size_t i = 0; i < variant_reads.size(); i++) {
            for (size_t j = 0; j < variant_reads[i].size(); j++) {
                if (frequency_matrix[i][j] < frequency_clamp) {
                    frequency_matrix[i][j] = frequency_clamp;
                } else if (frequency_matrix[i][j] > 1.0f - frequency_clamp) {
                    frequency_matrix[i][j] = 1.0f - frequency_clamp;
                }

                obj += compute_obj(frequency_matrix[i][j], variant_reads[i][j], total_reads[i][j]);
            }
        }

        objective = obj;
    }
}
