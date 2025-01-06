#include <spdlog/spdlog.h>

#include "../../CloneTree.h"
#include "../L2/Solver.h"
#include "Solver.h" 

namespace LogBinomialADMM {
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

        /* TODO:
         *  Take compute_gradient(.), compute_hessian(.), compute_obj(.) as an argument 
         *  to the Solver to support other loss functions.
         */
        const auto& compute_gradient = [](double freq, int var, int tot, double rho, double w) {
            if (var == 0) { 
                return ((tot - var) / (1 - freq)) + rho * (freq-w);
            } else if (var == tot) {
                return - var / freq + rho * (freq-w);
            } else {
                return ((tot - var) / (1 - freq)) - var / freq + rho * (freq-w);
            }
        };

        const auto& compute_hessian = [](double freq, int var, int tot, double rho, double w) {
            if (var == 0) {
                return rho + (tot - var) / ((1 - freq) * (1 - freq));
            } else if (var == tot) {
                return rho + var / (freq * freq);
            } else {
                return rho + (tot - var) / ((1 - freq) * (1 - freq)) + var / (freq * freq);
            }
        };

        const auto& compute_obj = [](double freq, int var, int tot, double rho, double w) {
            if (freq > 1.0 || freq < 0.0) return 1e9;

            if (var == 0) {
                return -(tot - var) * log(1 - freq) + rho * (freq - w) * (freq - w);
            } else if (var == tot) {
                return -var * log(freq) + rho * (freq - w) * (freq - w);
            } else {
                return -var * log(freq) - (tot - var) * log(1 - freq) + rho * (freq - w) * (freq - w);
            }
        };

        /* 
         * Solve 1D optimization problems using 
         * damped Newton's with backtracking line search. 
         * TODO: Investigate further.
         */
        for (size_t i = 0; i < frequencies.size(); i++) {
            for (size_t j = 0; j < frequencies[i].size(); j++) {
                double freq = 0.5;
                double current_obj = compute_obj(freq, variant_reads[i][j], total_reads[i][j], rho, ws[i][j]);
                for (int k = 0; k < max_newton_iterations; k++) {
                    double current_grad = compute_gradient(freq, variant_reads[i][j], total_reads[i][j], rho, ws[i][j]);
                    double current_hess = compute_hessian(freq, variant_reads[i][j], total_reads[i][j], rho, ws[i][j]);

                    double step_size = 1.0;
                    while (compute_obj(freq - current_grad * step_size / current_hess, variant_reads[i][j], total_reads[i][j], rho, ws[i][j]) > current_obj) {
                        step_size *= 0.90;
                    }

                    freq = freq - current_grad * step_size / current_hess;
                    double updated_obj = compute_obj(freq, variant_reads[i][j], total_reads[i][j], rho, ws[i][j]);

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
                float inc = 0.0;
                if (frequency_matrix[i][j] < frequency_clamp) {
                    frequency_matrix[i][j] = frequency_clamp;
                } else if (frequency_matrix[i][j] > 1.0f - frequency_clamp) {
                    frequency_matrix[i][j] = 1.0f - frequency_clamp;
                }

                if (variant_reads[i][j] == 0) {
                    inc = (total_reads[i][j] - variant_reads[i][j]) * log(1.0f - frequency_matrix[i][j]);
                } else if (variant_reads[i][j] == total_reads[i][j]) {
                    inc = variant_reads[i][j] * log(frequency_matrix[i][j]);
                } else {
                    inc = variant_reads[i][j] * log(frequency_matrix[i][j]) + (total_reads[i][j] - variant_reads[i][j]) * log(1.0f - frequency_matrix[i][j]);
                }

                obj += inc;
            }
        }

        objective = -obj;
    }
}
