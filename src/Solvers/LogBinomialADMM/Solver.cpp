#include <spdlog/spdlog.h>

#include "../../CloneTree.h"
#include "../L2/Solver.h"
#include "Solver.h" 
#include "Poly34.h"

namespace LogBinomialADMM {
    void Solver::solve() {
        std::vector<std::vector<double>> frequency_matrix;
        for (size_t i = 0; i < variant_reads.size(); i++) {
            std::vector<double> frequencies;
            for (size_t j = 0; j < variant_reads[i].size(); j++) {
                double freq = total_reads[i][j] == 0 ? 0 : static_cast<double>(variant_reads[i][j]) / total_reads[i][j];
                frequencies.push_back(freq);
            }
            frequency_matrix.push_back(frequencies);
        }

        // Step 1. Find decent initial solution
        L2Solver::Solver l2_solver(clone_tree, vertex_map, frequency_matrix, root);
        l2_solver.solve();

        frequencies = std::move(l2_solver.frequencies);
        usages = left_inverse(clone_tree, vertex_map, frequencies);
        lagrangian_objective_update();
        objective_update();

        for (int i = 0; i < 100; i++) {
            //spdlog::info("Objective: {}", objective);
            spdlog::info("Lagrangian Objective: {}", lagrangian_objective);

            frequency_update(); // ADMM Step 1
            usage_update();     // ADMM Step 2
            residual_update();  // ADMM Step 3

            lagrangian_objective_update();
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
                    if (x[l] >= 0 - 1e-6 && x[l] <= 1 + 1e-6) {
                        if (x[l] < 0) {
                            x[l] = 0;
                        } else if (x[l] > 1) {
                            x[l] = 1;
                        }

                        frequencies[i][j] = x[l];
                        found = true;
                        break;
                    }
                }

                if (!found) {
                    std::cout << "var/total: " << ((float) variant_reads[i][j]) / ((float) total_reads[i][j]) << std::endl;
                    std::cout << "ws[i][j]: " << ws[i][j] << std::endl;
                    std::cout << "frequencies[i][j]: " << frequencies[i][j] << std::endl;
                    for (size_t l = 0; l < k; l++) {
                        std::cout << "x[" << l << "]: " << x[l] << std::endl;
                    }
                    throw std::runtime_error("No solution found for cubic equation");
                }

                // std::cout << "var/total: " << ((float) variant_reads[i][j]) / ((float) total_reads[i][j]) << std::endl;
                // std::cout << "ws[i][j]: " << ws[i][j] << std::endl;
                // std::cout << "frequencies[i][j]: " << frequencies[i][j] << std::endl;
            }
        }
    }

    void Solver::usage_update() {
        std::vector<std::vector<double>> frequency_matrix = frequencies;
        for (size_t i = 0; i < frequency_matrix.size(); i++) {
            for (size_t j = 0; j < frequency_matrix[i].size(); j++) {
                frequency_matrix[i][j] = frequencies[i][j] - (residuals[i][j] / rho);
            }
        }
        L2Solver::Solver l2_solver(clone_tree, vertex_map, frequency_matrix, root);
        l2_solver.solve();
        usages = left_inverse(clone_tree, vertex_map, l2_solver.frequencies);

        for (size_t i = 0; i < l2_solver.frequencies.size(); i++) {
            for (size_t j = 0; j < l2_solver.frequencies[i].size(); j++) {
                // std::cout << "l2_solver.frequencies[i][j]: " << l2_solver.frequencies[i][j] << std::endl;
            }
        }
    }

    void Solver::residual_update() {
        auto ws = left_multiply(clone_tree, root, vertex_map, usages);
            
        for (size_t i = 0; i < ws.size(); i++) {
            for (size_t j = 0; j < ws[i].size(); j++) {
                residuals[i][j] += rho * (ws[i][j] - frequencies[i][j]);
                // std::cout << "residuals[i][j]: " << residuals[i][j] << std::endl;
            }
        }
    }

    void Solver::objective_update() {
        auto frequency_matrix = left_multiply(clone_tree, root, vertex_map, usages);
        float obj = 0.0;
        for (size_t i = 0; i < variant_reads.size(); i++) {
            for (size_t j = 0; j < variant_reads[i].size(); j++) {
                double inc = 0.0;

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

        objective = obj;
    }

    void Solver::lagrangian_objective_update() {
        auto frequency_matrix = left_multiply(clone_tree, root, vertex_map, usages);

        float obj = 0.0;
        for (size_t i = 0; i < variant_reads.size(); i++) {
            for (size_t j = 0; j < variant_reads[i].size(); j++) {
                if (variant_reads[i][j] == 0) {
                    obj += (total_reads[i][j] - variant_reads[i][j]) * log(1 - frequencies[i][j]);
                } else if (variant_reads[i][j] == total_reads[i][j]) {
                    obj += variant_reads[i][j] * log(frequencies[i][j]);
                } else {
                    obj += variant_reads[i][j] * log(frequencies[i][j]) + (total_reads[i][j] - variant_reads[i][j]) * log(1 - frequencies[i][j]);
                }

                obj *= -1;

                obj += (rho / 2) * pow(frequencies[i][j] - frequency_matrix[i][j], 2);
                obj += residuals[i][j] * (frequency_matrix[i][j] - frequencies[i][j]);
            }
        }

        lagrangian_objective = obj;
    }
}
