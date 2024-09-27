#include "Main.h"
#include "DiGraph.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

py::dict regress(
    const std::vector<std::vector<int>>& adjacency_list,
    const std::vector<std::vector<int>> &variant_matrix, 
    const std::vector<std::vector<int>> &total_matrix, 
    const std::optional<std::vector<std::vector<double>>>& weight_matrix,
    std::optional<int> root,
    const std::string loss_function
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

    std::unordered_map<int, int> vertex_map;
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
    if (root.has_value()) {
        root_value = root.value();
    } else {
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
    }

    std::vector<std::vector<double>> weight_matrix_value(variant_matrix.size(), std::vector<double>(variant_matrix[0].size(), 1.0));
    if (weight_matrix.has_value()) {
        weight_matrix_value = weight_matrix.value();
    }

    SolverResult r;
    if (loss_function == "l2") {
        r = l2_solve(vertex_map, variant_matrix, total_matrix, weight_matrix_value, clone_tree, root_value);
    } else if (loss_function == "binomial") {
        r = log_binomial_solve(vertex_map, variant_matrix, total_matrix, clone_tree, root_value, 10);
    } else if (loss_function == "l1") {
        throw std::invalid_argument("loss_function 'l1' is not yet implemented");
    } else {
        throw std::invalid_argument("loss_function must be one of ['l1', 'l2', 'binomial']");
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

PYBIND11_MODULE(fastppm, m) {
    m.doc() = "fastppm: fast perfect phylogeny mixture";
    m.def(
        "regress", &regress, "Regress the perfect phylogeny mixture model",
        py::arg("adjacency_list"), 
        py::arg("variant_matrix"), 
        py::arg("total_matrix"), 
        py::arg("weight_matrix") = std::nullopt,
        py::arg("root") = std::nullopt,
        py::arg("loss_function") = "l2"
    );
}