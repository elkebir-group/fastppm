#include <spdlog/spdlog.h>

#include "../../CloneTree.h"
#include "../L2/Solver.h"
#include "Solver.h" 

namespace SeparableADMM {
    void Solver::solve() {
        for (size_t i = 0; i < variant_reads.size(); i++) {
            for (size_t j = 0; j < variant_reads[i].size(); j++) {
                float freq = total_reads[i][j] == 0 ? 0 : static_cast<float>(variant_reads[i][j]) / static_cast<float>(total_reads[i][j]);
                frequencies[i][j] = freq;
                buffer[i][j] = freq;
            }
        }

        for (size_t i = 0; i < buffer.size(); i++) {
            for (size_t j = 0; j < buffer[i].size(); j++) {
                if (total_reads[i][j] == 0) {
                    buffer[i][j] = 0.0f;
                } else if (variant_reads[i][j] == 0 || variant_reads[i][j] == total_reads[i][j]) {
                    buffer[i][j] = total_reads[i][j];
                } else {
                    buffer[i][j] = ((float) variant_reads[i][j]) / (frequencies[i][j] * frequencies[i][j]);
                    buffer[i][j] += ((float) total_reads[i][j] - variant_reads[i][j]) / ((1.0f - frequencies[i][j]) * (1.0f - frequencies[i][j]));
                }
            }
        }

        // Step 1. Find decent initial solution
        l2_solver = L2Solver::Solver(clone_tree, vertex_map, frequencies, buffer, root);
        l2_solver.initialize();
        l2_solver.solve();

        frequencies = l2_solver.frequencies;
        usages = left_inverse(clone_tree, vertex_map, frequencies);
        objective_update();

        for (size_t i = 0; i < buffer.size(); i++) {
            for (size_t j = 0; j < buffer[i].size(); j++) {
                buffer[i][j] = 1.0f;
            }
        }

        l2_solver.update_weight_matrix(buffer);

        for (int i = 0; i < num_admm_iterations; i++) {
            left_multiply(clone_tree, root, vertex_map, usages, buffer);
            float primal_residual = 0.0;
            for (size_t i = 0; i < variant_reads.size(); i++) {
                for (size_t j = 0; j < variant_reads[i].size(); j++) {
                    primal_residual += std::abs(buffer[i][j] - frequencies[i][j]);
                }
            }

            primal_residual = primal_residual / (variant_reads.size() * sqrt(variant_reads[0].size()));

            if (verbose) {
                spdlog::info("ADMM Iteration: {}, Objective: {}, Normalized Primal Residual: {}", i, objective, primal_residual);
            }

            frequency_update(); // ADMM Step 1
            usage_update();     // ADMM Step 2
            residual_update();  // ADMM Step 3

            objective_update(); 
        }

        frequencies = left_multiply(clone_tree, root, vertex_map, usages);
    }

    void Solver::frequency_update() {
        auto &ws = buffer;
        left_multiply(clone_tree, root, vertex_map, usages, buffer);
        for (size_t i = 0; i < ws.size(); i++) {
            for (size_t j = 0; j < ws[i].size(); j++) {
                ws[i][j] += residuals[i][j];
            }
        }
        for (size_t i = 0; i < frequencies.size(); i++) {
            for (size_t j = 0; j < frequencies[i].size(); j++) {
                double freq = compute_minimizer(rho, ws[i][j], variant_reads[i][j], total_reads[i][j]);
                if (freq < frequency_clamp) {
                    freq = frequency_clamp;
                } else if (freq > 1.0f - frequency_clamp) {
                    freq = 1.0f - frequency_clamp;
                }

                frequencies[i][j] = freq;
            }
        }
    }

    void Solver::usage_update() {
        auto &frequency_matrix = buffer;
        for (size_t i = 0; i < frequency_matrix.size(); i++) {
            for (size_t j = 0; j < frequency_matrix[i].size(); j++) {
                frequency_matrix[i][j] = frequencies[i][j] - residuals[i][j];
            }
        }
        l2_solver.update_frequency_matrix(frequency_matrix);
        l2_solver.solve();
        left_inverse(clone_tree, vertex_map, l2_solver.frequencies, usages);
    }

    void Solver::residual_update() {
        auto &ws = buffer;
        left_multiply(clone_tree, root, vertex_map, usages, ws);
        for (size_t i = 0; i < ws.size(); i++) {
            for (size_t j = 0; j < ws[i].size(); j++) {
                residuals[i][j] += ws[i][j] - frequencies[i][j];
            }
        }
    }

    void Solver::objective_update() {
        auto &frequency_matrix = buffer;
        left_multiply(clone_tree, root, vertex_map, usages, buffer);
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
