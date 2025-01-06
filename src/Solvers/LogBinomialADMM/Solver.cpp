#include <spdlog/spdlog.h>

#include "../../CloneTree.h"
#include "../L2/Solver.h"
#include "Solver.h" 
#include "Poly34.h"

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
        L2Solver::Solver l2_solver(clone_tree, vertex_map, frequency_matrix, weight_matrix, root);
        l2_solver.solve();

        frequencies = std::move(l2_solver.frequencies);
        usages = left_inverse(clone_tree, vertex_map, frequencies);
        objective_update();

        for (int i = 0; i < 200; i++) {
            spdlog::info("Objective: {}", objective);

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

        double x[3];
        for (size_t i = 0; i < frequencies.size(); i++) {
            for (size_t j = 0; j < frequencies[i].size(); j++) {
                // solves cubic a + (b * x) + (c * x^2) + (d * x^3) = 0
                float a = -variant_reads[i][j];
                float b = total_reads[i][j] - rho * ws[i][j];
                float c = rho + rho*ws[i][j];
                float d = -rho;

                // solve the above cubic by making leading coefficient 1.0
                int k = SolveP3(x, c / d, b / d, a / d);

                bool found = false;
                for (int l = 0; l < k; l++) {
                    if ((x[l] >= 0.0f - 1e-4f) && (x[l] <= 1.0f + 1e-4f)) {
                        if (x[l] < 0.0f) {
                            x[l] = 0.0f;
                        } else if (x[l] > 1.0f) {
                            x[l] = 1.0f;
                        }

                        frequencies[i][j] = x[l];
                        found = true;
                        break;
                    }
                }

                if (!found) {
                    std::cout << "a + (b * x) + (c * x^2) + (d * x^3) = " << a << " + " << b << " * x + " << c << " * x^2 + " << d << " * x^3" << std::endl;
                    std::cout << "var/total: " << ((float) variant_reads[i][j]) / ((float) total_reads[i][j]) << std::endl;
                    std::cout << "ws[i][j]: " << ws[i][j] << std::endl;
                    std::cout << "frequencies[i][j]: " << frequencies[i][j] << std::endl;
                    for (size_t l = 0; l < k; l++) {
                        std::cout << "x[" << l << "]: " << x[l] << std::endl;
                    }
                    throw std::runtime_error("No solution found for cubic equation");
                }
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
        L2Solver::Solver l2_solver(clone_tree, vertex_map, frequency_matrix, weight_matrix, root);
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
                if (frequency_matrix[i][j] < 1e-5) {
                    frequency_matrix[i][j] = 1e-5;
                } else if (frequency_matrix[i][j] > 1 - 1e-5) {
                    frequency_matrix[i][j] = 1 - 1e-5;
                }

                /*
                if (frequency_matrix[i][j] == 0 && variant_reads[i][j] > 0) {
                    std::cout << "frequency_matrix[i][j]: " << frequency_matrix[i][j] << " variant_reads[i][j]: " << variant_reads[i][j] << " total_reads[i][j]: " << total_reads[i][j] << std::endl;
                    throw std::runtime_error("frequency_matrix[i][j] == 0 && variant_reads[i][j] > 0");
                }

                if (frequency_matrix[i][j] == 1 && variant_reads[i][j] < total_reads[i][j]) {
                    std::cout << "frequency_matrix[i][j]: " << frequency_matrix[i][j] << " variant_reads[i][j]: " << variant_reads[i][j] << " total_reads[i][j]: " << total_reads[i][j] << std::endl;
                    throw std::runtime_error("frequency_matrix[i][j] == 1 && variant_reads[i][j] < total_reads[i][j]");
                }

                std::cout << "frequency_matrix[i][j]: " << frequency_matrix[i][j] << " variant_reads[i][j]: " << variant_reads[i][j] << " total_reads[i][j]: " << total_reads[i][j] << std::endl;
                */

                if (variant_reads[i][j] == 0) {
                    inc = (total_reads[i][j] - variant_reads[i][j]) * log(1 - frequency_matrix[i][j]);
                } else if (variant_reads[i][j] == total_reads[i][j]) {
                    inc = variant_reads[i][j] * log(frequency_matrix[i][j]);
                } else {
                    inc = variant_reads[i][j] * log(frequency_matrix[i][j]) + (total_reads[i][j] - variant_reads[i][j]) * log(1 - frequency_matrix[i][j]);
                }

                obj += inc;
            }
        }

        objective = -obj;
    }
}
