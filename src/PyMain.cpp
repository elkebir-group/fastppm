#include "Lib.h"
#include "DiGraph.h"
#include "CloneTree.h"
#include "Solvers/L2/Solver.h"
#include "Solvers/CopyNumberADMM/Solver.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

py::dict regress_frequencies(
    const std::vector<std::vector<int>>& adjacency_list,
    const std::vector<std::vector<float>> &frequency_matrix, 
    const std::optional<std::vector<std::vector<float>>>& weight_matrix,
    const std::string loss_function
) {
    std::vector<int> vertex_map(adjacency_list.size());
    digraph<int> clone_tree;
    for (size_t i = 0; i < adjacency_list.size(); i++) {
        vertex_map[i] = clone_tree.add_vertex(i);
    }

    for (size_t p = 0; p < adjacency_list.size(); p++) {
        for (size_t j = 0; j < adjacency_list[p].size(); j++) {
            size_t c = adjacency_list[p][j];
            clone_tree.add_edge(vertex_map[p], vertex_map[c]);
        }
    }

    int root_value;
    bool seen_root = false;
    for (auto n : clone_tree.nodes()) {
        if (clone_tree.in_degree(n) == 0) {
            root_value = clone_tree[n].data;
            seen_root = true;
            break;
        }
    }

    if (!seen_root) {
        throw std::invalid_argument("no root node found in the adjacency list");
    }

    std::vector<std::vector<float>> weight_matrix_value(frequency_matrix.size(), std::vector<float>(frequency_matrix[0].size(), 1.0));
    if (weight_matrix.has_value()) {
        weight_matrix_value = weight_matrix.value();
    }

    SolverResult r;
    if (loss_function == "l2") {
        L2Solver::Solver solver(clone_tree, vertex_map, frequency_matrix, weight_matrix_value, root_value);
        solver.initialize();
        solver.solve();
        auto usage_matrix = left_inverse(clone_tree, vertex_map, solver.frequencies);
        r = {-1.0, solver.objective, usage_matrix, solver.frequencies};
    } else {
        throw std::invalid_argument("loss_function must be one of ['l2']");
    }

    py::dict result;
    result["objective"] = r.objective;
    if (r.usage_matrix.has_value()) {
        result["usage_matrix"] = r.usage_matrix.value();
    }

    if (r.frequency_matrix.has_value()) {
        result["frequency_matrix"] = r.frequency_matrix.value();
    }

    return result;
}

py::list regress_read_counts_multiple_trees(
        const std::vector<std::vector<std::vector<int>>>& adjacency_lists,
        const std::vector<std::vector<int>> &variant_matrix,
        const std::vector<std::vector<int>> &total_matrix,
        const std::optional<std::vector<std::vector<float>>>& weight_matrix,
        const std::string loss_function,
        int nr_segments,
        float precision
) {
    if (variant_matrix.size() != total_matrix.size()) {
        throw std::invalid_argument("variant_matrix and total_matrix must have the same number of rows");
    }

    if (variant_matrix.size() == 0) {
        throw std::invalid_argument("variant_matrix and total_matrix must be non-empty");
    }

    if (variant_matrix[0].size() != total_matrix[0].size()) {
        throw std::invalid_argument("variant_matrix and total_matrix must have the same number of columns");
    }

    // Initialize result list
    py::list result;

    // Iterate over each tree in the list
    for (const auto& adjacency_list : adjacency_lists) {
        std::vector<int> vertex_map(adjacency_list.size());
        digraph<int> clone_tree;
        for (size_t i = 0; i < adjacency_list.size(); i++) {
            vertex_map[i] = clone_tree.add_vertex(i);
        }

        for (size_t p = 0; p < adjacency_list.size(); p++) {
            for (size_t j = 0; j < adjacency_list[p].size(); j++) {
                size_t c = adjacency_list[p][j];
                clone_tree.add_edge(vertex_map[p], vertex_map[c]);
            }
        }

        int root_value;
        bool seen_root = false;
        for (auto n: clone_tree.nodes()) {
            if (clone_tree.in_degree(n) == 0) {
                root_value = clone_tree[n].data;
                seen_root = true;
                break;
            }
        }

        if (!seen_root) {
            throw std::invalid_argument("no root node found in the adjacency list");
        }

        std::vector<std::vector<float>> weight_matrix_value(variant_matrix.size(),
                                                            std::vector<float>(variant_matrix[0].size(), 1.0));
        if (weight_matrix.has_value()) {
            weight_matrix_value = weight_matrix.value();
        }

        SolverResult r;
        if (loss_function == "l2") {
            r = l2_solve(vertex_map, variant_matrix, total_matrix, weight_matrix_value, clone_tree, root_value);
        } else if (loss_function == "binomial") {
            r = log_binomial_admm_solve(vertex_map, variant_matrix, total_matrix, clone_tree, root_value);
        } else if (loss_function == "binomial_pla") {
            r = log_binomial_fixed_solve(vertex_map, variant_matrix, total_matrix, clone_tree, root_value, nr_segments);
        } else if (loss_function == "binomial_ppla") {
            r = log_binomial_solve(vertex_map, variant_matrix, total_matrix, clone_tree, root_value, nr_segments);
        } else if (loss_function == "beta_binomial_pla") {
            r = log_beta_binomial_solve(vertex_map, variant_matrix, total_matrix, clone_tree, root_value, nr_segments,
                                        precision);
        } else if (loss_function == "beta_binomial_ppla") {
            r = log_beta_binomial_fixed_solve(vertex_map, variant_matrix, total_matrix, clone_tree, root_value,
                                              nr_segments, precision);
        } else if (loss_function == "l1") {
            throw std::invalid_argument("loss_function 'l1' is not yet implemented");
        } else {
            throw std::invalid_argument(
                    "loss_function must be one of ['l1', 'l2', 'binomial', 'binomial_pla', 'binomial_ppla', 'beta_binomial_pla', 'beta_binomial_ppla']");
        }

        // Store the result for this tree in a dictionary
        py::dict tree_result;
        tree_result["objective"] = r.objective;
        if (r.usage_matrix.has_value()) {
            tree_result["usage_matrix"] = r.usage_matrix.value();
        }
        if (r.frequency_matrix.has_value()) {
            tree_result["frequency_matrix"] = r.frequency_matrix.value();
        }

        // Append the tree result to the result list
        result.append(tree_result);
    }

    return result;

}

py::dict regress_read_counts(
    const std::vector<std::vector<int>>& adjacency_list,
    const std::vector<std::vector<int>> &variant_matrix, 
    const std::vector<std::vector<int>> &total_matrix, 
    const std::optional<std::vector<std::vector<float>>>& weight_matrix,
    const std::string loss_function,
    int nr_segments,
    float precision
) {
    if (variant_matrix.size() != total_matrix.size()) {
        throw std::invalid_argument("variant_matrix and total_matrix must have the same number of rows");
    }

    if (variant_matrix.size() == 0) {
        throw std::invalid_argument("variant_matrix and total_matrix must be non-empty");
    }

    if (variant_matrix[0].size() != total_matrix[0].size()) {
        throw std::invalid_argument("variant_matrix and total_matrix must have the same number of columns");
    }

    std::vector<int> vertex_map(adjacency_list.size());
    digraph<int> clone_tree;
    for (size_t i = 0; i < adjacency_list.size(); i++) {
        vertex_map[i] = clone_tree.add_vertex(i);
    }

    for (size_t p = 0; p < adjacency_list.size(); p++) {
        for (size_t j = 0; j < adjacency_list[p].size(); j++) {
            size_t c = adjacency_list[p][j];
            clone_tree.add_edge(vertex_map[p], vertex_map[c]);
        }
    }

    int root_value;
    bool seen_root = false;
    for (auto n : clone_tree.nodes()) {
        if (clone_tree.in_degree(n) == 0) {
            root_value = clone_tree[n].data;
            seen_root = true;
            break;
        }
    }

    if (!seen_root) {
        throw std::invalid_argument("no root node found in the adjacency list");
    }

    std::vector<std::vector<float>> weight_matrix_value(variant_matrix.size(), std::vector<float>(variant_matrix[0].size(), 1.0));
    if (weight_matrix.has_value()) {
        weight_matrix_value = weight_matrix.value();
    }

    SolverResult r;
    if (loss_function == "l2") {
        r = l2_solve(vertex_map, variant_matrix, total_matrix, weight_matrix_value, clone_tree, root_value);
    } else if (loss_function == "binomial") {
        r = log_binomial_admm_solve(vertex_map, variant_matrix, total_matrix, clone_tree, root_value);
    } else if (loss_function == "binomial_pla") {
        r = log_binomial_fixed_solve(vertex_map, variant_matrix, total_matrix, clone_tree, root_value, nr_segments);
    } else if (loss_function == "binomial_ppla") {
        r = log_binomial_solve(vertex_map, variant_matrix, total_matrix, clone_tree, root_value, nr_segments);
    } else if (loss_function == "beta_binomial_pla") {
        r = log_beta_binomial_solve(vertex_map, variant_matrix, total_matrix, clone_tree, root_value, nr_segments, precision);
    } else if (loss_function == "beta_binomial_ppla") {
        r = log_beta_binomial_fixed_solve(vertex_map, variant_matrix, total_matrix, clone_tree, root_value, nr_segments, precision);
    } else if (loss_function == "l1") {
        throw std::invalid_argument("loss_function 'l1' is not yet implemented");
    } else {
        throw std::invalid_argument("loss_function must be one of ['l1', 'l2', 'binomial', 'binomial_pla', 'binomial_ppla', 'beta_binomial_pla', 'beta_binomial_ppla']");
    }

    py::dict result;
    result["objective"] = r.objective;
    if (r.usage_matrix.has_value()) {
        result["usage_matrix"] = r.usage_matrix.value();
    }

    if (r.frequency_matrix.has_value()) {
        result["frequency_matrix"] = r.frequency_matrix.value();
    }

    return result;

}

py::dict regress_copy_number(
    const std::vector<int>& parent,
    const std::vector<std::vector<std::pair<int, int>>>& states,
    const py::array_t<int, py::array::c_style | py::array::forcecast>& multiplicity,
    const py::array_t<int, py::array::c_style | py::array::forcecast>& segment,
    const py::array_t<int, py::array::c_style | py::array::forcecast>& variant,
    const py::array_t<int, py::array::c_style | py::array::forcecast>& total,
    const std::vector<std::vector<std::tuple<int, int, double>>>& state_fractions,
    const std::optional<std::vector<double>>& mean_copy_number,
    const std::optional<std::vector<double>>& segment_weight,
    double rho,
    double kappa,
    double eps,
    int max_iterations,
    int check_every,
    int polish_steps,
    double alpha,
    bool verbose,
    bool warm_prox,
    int anderson_memory,
    double anderson_reg,
    const std::optional<std::vector<double>>& initial_usages
) {
    // SNV inputs are NumPy arrays (lists are converted): one copy instead of a conversion per element
    if (multiplicity.ndim() != 2 || multiplicity.shape(1) != static_cast<py::ssize_t>(parent.size()))
        throw std::invalid_argument("multiplicity must be an array of shape (number of SNVs, number of clones)");
    const py::ssize_t n = multiplicity.shape(0), k = multiplicity.shape(1);
    if (segment.size() != n || variant.size() != n || total.size() != n)
        throw std::invalid_argument("segment, variant and total must have one entry per SNV");
    std::vector<std::vector<int>> mult(n);
    const int* m = multiplicity.data();
    for (py::ssize_t i = 0; i < n; ++i) mult[i].assign(m + i * k, m + (i + 1) * k);
    auto vec = [](const py::array_t<int, py::array::c_style | py::array::forcecast>& a) {
        return std::vector<int>(a.data(), a.data() + a.size());
    };
    CopyNumberADMM::Instance instance{parent, states, mult, vec(segment), vec(variant), vec(total), state_fractions,
                                      mean_copy_number.value_or(std::vector<double>()),
                                      segment_weight.value_or(std::vector<double>())};
    CopyNumberADMM::Options options{rho, kappa, eps, max_iterations, check_every, polish_steps, alpha, verbose, anderson_memory, anderson_reg, warm_prox, {}};
    if (initial_usages) options.initial_usages = *initial_usages;
    CopyNumberADMM::Result r = CopyNumberADMM::solve(instance, options);

    py::dict result;
    result["usages"] = r.usages;
    result["frequencies"] = r.frequencies;
    result["loss"] = r.loss;
    result["lower_bound"] = r.lower_bound;
    result["gap"] = r.loss - r.lower_bound;
    result["iterations"] = r.iterations;
    result["kappa"] = r.kappa;
    result["anderson_accepted"] = r.anderson_accepted;
    result["anderson_rejected"] = r.anderson_rejected;
    result["converged"] = r.converged;
    result["runtime"] = r.runtime_ms;
    return result;
}

PYBIND11_MODULE(fastppm, m) {
    m.doc() = "fastppm: fast perfect phylogeny mixture";

    std::string version = std::to_string(FASTPPM_VERSION_MAJOR) + "." +
                          std::to_string(FASTPPM_VERSION_MINOR) + "." +
                          std::to_string(FASTPPM_VERSION_PATCH);

    m.attr("__version__") = version;

    m.def(
        "regress_counts", &regress_read_counts, "Regresses read counts against the perfect phylogeny mixture model",
        py::arg("adjacency_list"), 
        py::arg("variant_matrix"), 
        py::arg("total_matrix"), 
        py::arg("weight_matrix") = std::nullopt,
        py::arg("loss_function") = "binomial",
        py::arg("nr_segments") = 10,
        py::arg("precision") = 2.0
    );

    m.def(
        "regress_counts_multiple_trees", &regress_read_counts_multiple_trees, "Regresses read counts against the perfect phylogeny mixture model for multiple trees",
        py::arg("adjacency_lists"),
        py::arg("variant_matrix"),
        py::arg("total_matrix"),
        py::arg("weight_matrix") = std::nullopt,
        py::arg("loss_function") = "binomial",
        py::arg("nr_segments") = 10,
        py::arg("precision") = 2.0
    );

    m.def(
        "regress_frequencies", &regress_frequencies, "Regresses frequencies against the perfect phylogeny mixture model",
        py::arg("adjacency_list"), 
        py::arg("frequency_matrix"), 
        py::arg("weight_matrix") = std::nullopt,
        py::arg("loss_function") = "l2"
    );

    m.def(
        "regress_copy_number", &regress_copy_number,
        "Copy-number-aware clonal regression (fastppm-CN): binomial loss on SNVs, squared loss on state fractions",
        py::arg("parent"),
        py::arg("states"),
        py::arg("multiplicity"),
        py::arg("segment"),
        py::arg("variant"),
        py::arg("total"),
        py::arg("state_fractions"),
        py::arg("mean_copy_number") = std::nullopt,
        py::arg("segment_weight") = std::nullopt,
        py::arg("rho") = 100.0,
        py::arg("kappa") = -1.0,
        py::arg("eps") = 1.0,
        py::arg("max_iterations") = 100000,
        py::arg("check_every") = 10,
        py::arg("polish_steps") = 1,
        py::arg("alpha") = 1.6,
        py::arg("verbose") = false,
        py::arg("warm_prox") = true,
        py::arg("anderson_memory") = 10,
        py::arg("anderson_reg") = 1e-10,
        py::arg("initial_usages") = std::nullopt
    );
}
