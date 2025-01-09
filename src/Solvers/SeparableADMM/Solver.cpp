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
        // l2_solver.solve();

        // frequencies = l2_solver.frequencies;
        // usages = left_inverse(clone_tree, vertex_map, frequencies);
        // residual_update();
        // objective_update();

        for (int i = 0; i < num_admm_iterations; i++) {
            std::cout << "Objective: " << objective << std::endl;
            auto ws = left_multiply(clone_tree, root, vertex_map, usages);
            float primal_residual = 0.0;
            for (size_t i = 0; i < variant_reads.size(); i++) {
                for (size_t j = 0; j < variant_reads[i].size(); j++) {
                    primal_residual += std::abs(ws[i][j] - frequencies[i][j]);
                }
            }
            std::cout << "Primal Residual: " << primal_residual << std::endl;

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
                ws[i][j] += residuals[i][j];
            }
        }
        for (size_t i = 0; i < frequencies.size(); i++) {
            for (size_t j = 0; j < frequencies[i].size(); j++) {
                double freq = compute_minimizer(rho, ws[i][j], variant_reads[i][j], total_reads[i][j]);
                // std::cout << "Frequency: " << freq << std::endl;
                // std::cout << "RHO: " << rho << std::endl;
                // std::cout << "WS: " << ws[i][j] << std::endl;
                // std::cout << "Variant Reads: " << variant_reads[i][j] << std::endl;
                // std::cout << "Total Reads: " << total_reads[i][j] << std::endl;
                // std::cout << "-------------------" << std::endl;
                if (freq < frequency_clamp) {
                    freq = frequency_clamp;
                } else if (freq > 1.0f - frequency_clamp) {
                    freq = 1.0f - frequency_clamp;
                }

                frequencies[i][j] = freq;
            }
        }
    }

    void Solver::usage_update(const std::vector<std::vector<float>>& weight_matrix) {
        std::vector<std::vector<float>> frequency_matrix = frequencies;
        for (size_t i = 0; i < frequency_matrix.size(); i++) {
            for (size_t j = 0; j < frequency_matrix[i].size(); j++) {
                frequency_matrix[i][j] -= residuals[i][j];
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
                residuals[i][j] += ws[i][j] - frequencies[i][j];
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
