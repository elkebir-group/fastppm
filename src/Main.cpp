#include "Lib.h"

#include <vector>
#include <argparse/argparse.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/fmt/ostr.h>
#include <yyjson.h>

#define FASTPPM_VERSION_MAJOR 1
#define FASTPPM_VERSION_MINOR 0

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

    program.add_argument("-w", "--weights")
        .help("Path to the weights matrix file")
        .default_value("");

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

    program.add_argument("-f", "--format")
        .help("Output format, either 'concise' or 'verbose'")
        .default_value("concise")
        .choices("concise", "verbose");

    program.add_argument("-l", "--loss")
        .help("Loss function L_i(.) to use for optimization")
        .default_value("l2")
        .choices("l1", "l2", "binomial", "binomial_K");

    program.add_argument("-K", "--segments")
        .help("Number of segments, only used when loss function is 'binomial' or 'binomial_K'")
        .default_value(10)
        .scan<'d', int>();

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

    std::vector<std::vector<double>> weights(variant_matrix.size(), std::vector<double>(variant_matrix[0].size(), 1.0));
    if (!program.get<std::string>("-w").empty()) {
        weights = parse_matrix<double>(program.get<std::string>("-w"));
        if (weights.size() != variant_matrix.size()) {
            error_logger->error("The weights matrix has different number of rows than the variant read matrix.");
            std::exit(1);
        }

        for (size_t i = 0; i < weights.size(); i++) {
            if (weights[i].size() != variant_matrix[i].size()) {
                error_logger->error("The weights matrix has different number of columns than the variant read matrix.");
                std::exit(1);
            }
        }

        for (size_t i = 0; i < weights.size(); i++) {
            for (size_t j = 0; j < weights[i].size(); j++) {
                if (weights[i][j] < 0) {
                    error_logger->error("The weights matrix contains negative values.");
                    std::exit(1);
                }
            }
        }
    }

    // vertex_map : takes the vertex ID in the adjacency list to the vertex ID in the digraph
    auto [clone_tree, vertex_map] = parse_adjacency_list(program.get<std::string>("tree"));
    size_t tree_root = program.get<int>("-r");
    size_t nr_segments = program.get<int>("-K");

    SolverResult result;
    if (program.get<std::string>("-l") == "binomial") {
        result = log_binomial_solve(vertex_map, variant_matrix, total_matrix, clone_tree, tree_root, nr_segments);
    } else if (program.get<std::string>("-l") == "binomial_K") {
        result = log_binomial_fixed_solve(vertex_map, variant_matrix, total_matrix, clone_tree, tree_root, nr_segments);
    } else if (program.get<std::string>("-l") == "l2") {
            result = l2_solve(vertex_map, variant_matrix, total_matrix, weights, clone_tree, tree_root);
    } else {
        error_logger->error("The loss function specified is not yet supported.");
        std::exit(1);
    }

    // Write the output to a JSON file
    std::ofstream output_file(program.get<std::string>("-o"));

    yyjson_mut_doc *doc = yyjson_mut_doc_new(NULL);
    yyjson_mut_val *root = yyjson_mut_obj(doc);
    yyjson_mut_doc_set_root(doc, root);

    yyjson_mut_obj_add_real(doc, root, "objective", result.objective);
    yyjson_mut_obj_add_real(doc, root, "runtime", result.runtime);

    float rounding_factor = 1e5;

    bool verbose = program.get<std::string>("-f") == "verbose";

    if (verbose && result.usage_matrix.has_value()) {
        auto &usage_matrix = result.usage_matrix.value();
        for (size_t i = 0; i < usage_matrix.size(); i++) {
            for (size_t j = 0; j < usage_matrix[i].size(); j++) {
                usage_matrix[i][j] = std::round(usage_matrix[i][j] * rounding_factor) / rounding_factor;

                if (usage_matrix[i][j] == 0 && std::signbit(usage_matrix[i][j])) {
                    usage_matrix[i][j] = 0;
                }
            }
        }

        yyjson_mut_val *usage_matrix_json = yyjson_mut_arr(doc);
        for (size_t i = 0; i < usage_matrix.size(); i++) {
            yyjson_mut_val *row_json = yyjson_mut_arr(doc);
            for (size_t j = 0; j < usage_matrix[i].size(); j++) {
                yyjson_mut_arr_add_real(doc, row_json, usage_matrix[i][j]);
            }
            yyjson_mut_arr_add_val(usage_matrix_json, row_json);
        }
        yyjson_mut_obj_add_val(doc, root, "usage_matrix", usage_matrix_json);
    }

    if (verbose && result.frequency_matrix.has_value()) {
        auto &frequency_matrix = result.frequency_matrix.value();
        for (size_t i = 0; i < frequency_matrix.size(); i++) {
            for (size_t j = 0; j < frequency_matrix[i].size(); j++) {
                frequency_matrix[i][j] = std::round(frequency_matrix[i][j] * rounding_factor) / rounding_factor;

                if (frequency_matrix[i][j] == 0 && std::signbit(frequency_matrix[i][j])) {
                    frequency_matrix[i][j] = 0;
                }
            }
        }

        yyjson_mut_val *frequency_matrix_json = yyjson_mut_arr(doc);
        for (size_t i = 0; i < frequency_matrix.size(); i++) {
            yyjson_mut_val *row_json = yyjson_mut_arr(doc);
            for (size_t j = 0; j < frequency_matrix[i].size(); j++) {
                yyjson_mut_arr_add_real(doc, row_json, frequency_matrix[i][j]);
            }
            yyjson_mut_arr_add_val(frequency_matrix_json, row_json);
        }
        yyjson_mut_obj_add_val(doc, root, "frequency_matrix", frequency_matrix_json);
    }

    size_t len;
    char *json_str = yyjson_mut_write(doc, 0, &len);

    if (json_str) {
        output_file.write(json_str, len);
        free(json_str);
    } else {
        std::cerr << "Failed to serialize JSON document." << std::endl;
    }

    yyjson_mut_doc_free(doc);
    return 0;
}
