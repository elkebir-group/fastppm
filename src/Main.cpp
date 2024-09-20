#include <vector>
#include <fstream>
#include <iostream>
#include <chrono>

#include <argparse/argparse.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/fmt/ostr.h>
#include <nlohmann/json.hpp>

#include "CloneTree.h"
#include "DiGraph.h"
#include "Solvers/LogBinomialPiecewise/Solver.h"
#include "Solvers/L2/Solver.h"
#include "Solvers/LogBinomialADMM/Solver.h"

#define FASTPPM_VERSION_MAJOR 1
#define FASTPPM_VERSION_MINOR 0

/* 
 * Represents the solver output.
 * 
 * Fields:
 *  - runtime: The time taken to solve the optimization problem in milliseconds.
 *  - objective: The value of the objective function after optimization.
 *  - usage_matrix: The usage matrix after optimization. This is an optional field.
 *  - frequency_matrix: The frequency matrix after optimization. This is an optional field.
 */
struct SolverResult {
    double runtime;
    double objective;
    std::optional<std::vector<std::vector<double>>> usage_matrix;
    std::optional<std::vector<std::vector<double>>> frequency_matrix;
};

/* 
 * Solves the optimization problem using the L2 loss function.
 */
SolverResult l2_solve(
    const std::unordered_map<int, int>& vertex_map,
    const std::vector<std::vector<int>>& variant_matrix,
    const std::vector<std::vector<int>>& total_matrix,
    const digraph<int>& clone_tree,
    size_t root
) {
    std::vector<std::vector<double>> frequency_matrix;
    for (size_t i = 0; i < variant_matrix.size(); i++) {
        std::vector<double> frequencies;
        for (size_t j = 0; j < variant_matrix[i].size(); j++) {
            double freq = total_matrix[i][j] == 0 ? 0 : static_cast<double>(variant_matrix[i][j]) / total_matrix[i][j];
            frequencies.push_back(freq);
        }
        frequency_matrix.push_back(frequencies);
    }

    L2Solver::Solver solver(clone_tree, vertex_map, frequency_matrix, root);

    auto start = std::chrono::high_resolution_clock::now();
    solver.solve();
    auto end = std::chrono::high_resolution_clock::now();
    double runtime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    auto usage_matrix = left_inverse(clone_tree, vertex_map, solver.frequencies);
    return {runtime, solver.objective, usage_matrix, solver.frequencies};
}

SolverResult log_binomial_admm_solve(
    const std::unordered_map<int, int>& vertex_map,
    const std::vector<std::vector<int>>& variant_matrix,
    const std::vector<std::vector<int>>& total_matrix,
    const digraph<int>& clone_tree,
    size_t root
) {
    LogBinomialADMM::Solver solver(clone_tree, vertex_map, variant_matrix, total_matrix, root);

    auto start = std::chrono::high_resolution_clock::now();
    solver.solve();
    auto end = std::chrono::high_resolution_clock::now();
    double runtime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    return {runtime, solver.objective, solver.usages, solver.frequencies};
}

/* 
 * Solves the optimization problem using the binomial loss function
 * with Yuanyuan's progressive piecewise linear solver.
 */
SolverResult log_binomial_solve(
    const std::unordered_map<int, int>& vertex_map,
    const std::vector<std::vector<int>>& variant_matrix,
    const std::vector<std::vector<int>>& total_matrix,
    const digraph<int>& clone_tree,
    size_t root,
    int K
) {
    auto edges = clone_tree.edges();
    int n_clones = edges.size() + 1;
    std::vector<std::list<int>> link_list(n_clones);
    for (auto& [u, v] : edges) {
        link_list[clone_tree[u].data].push_back(clone_tree[v].data);
    }

    std::vector<std::vector<double>> frequency_matrix;

    auto start = std::chrono::high_resolution_clock::now();
    double objective = 0;
    for (size_t i = 0; i < variant_matrix.size(); i++) {
        std::vector<int> ref_vector(variant_matrix[i].size(), 0);
        for (size_t j = 0; j < variant_matrix[i].size(); j++) {
            ref_vector[j] = total_matrix[i][j] - variant_matrix[i][j];
        }

        LogBinomialPiecewiseLinearSolver::Solver solver(K);
        solver.init(variant_matrix[i], ref_vector, link_list, root);
        objective += solver.main(0.75, 1e-6); // TODO: make these parameters configurable
        
        std::vector<double> frequencies(variant_matrix[i].size(), 0);
        for (size_t j = 0; j < variant_matrix[i].size(); j++) {
            frequencies[j] = solver.F[j];
        }

        frequency_matrix.push_back(frequencies);
    }
    auto end = std::chrono::high_resolution_clock::now();
    double runtime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    
    auto usage_matrix = left_inverse(clone_tree, vertex_map, frequency_matrix);
    return {runtime, objective, usage_matrix, frequency_matrix};
}

/* 
 * This function parses a matrix from a text file and returns it as a 2D vector.
 *
 * Input: 
 *  - filename: The name of the file that contains the matrix.
 *              Each line in the file should correspond to a row in the matrix.
 *              Within a line, numbers should be space-separated and represent the values 
 *              in the columns of that row.
 * 
 * Returns a 2D vector (vector of vectors) that represents the matrix. Each inner 
 * vector represents a row of the matrix. Throws a std::runtime_error if the file can't be opened
 * or if the file does not represent a valid matrix.
*/
template <typename T>
std::vector<std::vector<T>> parse_matrix(const std::string& filename) {
    std::ifstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Failed to open the file.");
    }

    std::vector<std::vector<T>> matrix;
    std::string line;
    size_t num_cols = 0;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::vector<T> row;
        T value;

        while (iss >> value) {
            row.push_back(value);
        }

        if (matrix.empty()) {
            num_cols = row.size();
        } else if (row.size() != num_cols) {
            throw std::runtime_error("The file does not represent a matrix. Number of columns is not consistent across rows.");
        }

        matrix.push_back(row);
    }

    file.close();
    return matrix;
}

int main(int argc, char ** argv) {
    // Set up logging with spdlog
    auto console_logger = spdlog::stdout_color_mt("fastppm");
    auto error_logger = spdlog::stderr_color_mt("error");
    spdlog::set_default_logger(console_logger);

    // Parse command line arguments
    argparse::ArgumentParser program(
        "fastppm",
        std::to_string(FASTPPM_VERSION_MAJOR) + "." + std::to_string(FASTPPM_VERSION_MAJOR),
        argparse::default_arguments::help
    );

    program.add_argument("-v", "--variant")
        .help("Path to the variant read matrix file")
        .required();

    program.add_argument("--version")
        .action([&](const auto & /*unused*/) {
            std::cout << "fastppm version " << FASTPPM_VERSION_MAJOR << "." << FASTPPM_VERSION_MINOR << std::endl;
            std::exit(0);
        })
        .default_value(false)
        .help("prints version information and exits")
        .implicit_value(true)
        .nargs(0);

    program.add_argument("-d", "--total")
        .help("Path to the total read matrix file")
        .required();

    program.add_argument("-t", "--tree")
        .help("Path to the tree file")
        .required();

    program.add_argument("-o", "--output")
        .help("Path to the output file")
        .required();

    program.add_argument("-r", "--root")
        .help("Root node of the tree")
        .default_value(0)
        .scan<'d', int>();

    program.add_argument("-l", "--loss")
        .help("Loss function L_i(.) to use for optimization")
        .default_value("binomial")
        .choices("l1", "l2", "binomial", "binomial_admm");

    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        std::exit(1);
    }

    auto variant_matrix = parse_matrix<int>(program.get<std::string>("-v"));
    auto total_matrix = parse_matrix<int>(program.get<std::string>("-d"));

    if (variant_matrix.size() != total_matrix.size()) {
        error_logger->error("The variant and total read matrices have different number of rows.");
        std::exit(1);
    }

    if (variant_matrix.empty()) {
        error_logger->error("The variant read matrix is empty.");
        std::exit(1);
    }

    if (total_matrix.empty()) {
        error_logger->error("The total read matrix is empty.");
        std::exit(1);
    }

    for (size_t i = 0; i < variant_matrix.size(); i++) {
        if (variant_matrix[i].size() != total_matrix[i].size()) {
            error_logger->error("The variant and total read matrices have different number of columns.");
            std::exit(1);
        }
    }

    // vertex_map : takes the vertex ID in the adjacency list to the vertex ID in the digraph
    auto [clone_tree, vertex_map] = parse_adjacency_list(program.get<std::string>("tree"));
    size_t root = program.get<int>("-r");

    SolverResult result;
    if (program.get<std::string>("-l") == "binomial") {
        result = log_binomial_solve(vertex_map, variant_matrix, total_matrix, clone_tree, root, 10);
    } else if (program.get<std::string>("-l") == "l2") {
        result = l2_solve(vertex_map, variant_matrix, total_matrix, clone_tree, root);
    } else if (program.get<std::string>("-l") == "binomial_admm") {
        result = log_binomial_admm_solve(vertex_map, variant_matrix, total_matrix, clone_tree, root);
    } else {
        error_logger->error("The loss function specified is not supported.");
        std::exit(1);
    }

    // Write the output to a JSON file
    nlohmann::json output;
    output["objective"] = result.objective;
    output["runtime"] = result.runtime;

    if (result.usage_matrix.has_value()) {
        output["usage_matrix"] = result.usage_matrix.value();
    }

    if (result.frequency_matrix.has_value()) {
        output["frequency_matrix"] = result.frequency_matrix.value();
    }

    std::ofstream output_file(program.get<std::string>("-o"));
    output_file << output.dump();

    return 0;
}
