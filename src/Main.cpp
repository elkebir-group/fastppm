#include "Lib.h"
#include "fast_float/fast_float.h"

#include <cstdlib>
#include <vector>
#include <argparse/argparse.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/fmt/ostr.h>
#include <yyjson.h>

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
    std::ifstream file(filename, std::ios::binary | std::ios::ate);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open the file.");
    }
    
    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);

    std::string buffer(size, ' ');
    if (!file.read(&buffer[0], size)) {
        throw std::runtime_error("Failed to read the file into memory.");
    }

    file.close();

    std::vector<std::vector<T>> matrix;
    std::vector<T> row;
    const char* ptr = buffer.c_str();
    const char* end = ptr + size;
    size_t num_cols = 0;

    while (ptr < end) {
        while (ptr < end && std::isspace(*ptr)) {
            if (*ptr == '\n' || *ptr == '\r') {
                if (!row.empty()) {
                    if (matrix.empty()) {
                        num_cols = row.size();
                    } else if (row.size() != num_cols) {
                        throw std::runtime_error("Inconsistent row sizes in matrix file.");
                    }
                    matrix.push_back(std::move(row));
                    row.clear();
                }
            }
            ++ptr;
        }

        if (ptr >= end) break;

        T value;
        auto res = fast_float::from_chars(ptr, end, value);
        row.push_back(value);
        ptr = res.ptr;
    }

    if (!row.empty()) {
        if (matrix.empty()) {
            num_cols = row.size();
        } else if (row.size() != num_cols) {
            throw std::runtime_error("Inconsistent row sizes in matrix file.");
        }
        matrix.push_back(std::move(row));
    }

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
        std::to_string(FASTPPM_VERSION_MAJOR)
            + "." + std::to_string(FASTPPM_VERSION_MAJOR)
            + "." + std::to_string(FASTPPM_VERSION_PATCH),
        argparse::default_arguments::help
    );

    program.add_argument("-v", "--variant")
        .help("Path to the variant read matrix file")
        .required();

    program.add_argument("--version")
        .action([&](const auto & /*unused*/) {
            std::cout << "fastppm version " << FASTPPM_VERSION_MAJOR
                << "." << FASTPPM_VERSION_MINOR
                << "." << FASTPPM_VERSION_PATCH
                << std::endl;
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

    program.add_argument("-f", "--format")
        .help("Output format, either 'concise' or 'verbose'")
        .default_value("concise")
        .choices("concise", "verbose");

    program.add_argument("-l", "--loss")
        .help("Loss function L_i(.) to use for optimization")
        .default_value("l2")
        .choices("l2", "binomial", "binomial_pla", "binomial_ppla", "beta_binomial_pla", "beta_binomial_ppla");

    program.add_argument("-s", "--precision")
        .help("Precision parameter, only used when loss function is 'beta_binomial*'")
        .default_value(10.)
        .scan<'g', double>();

    program.add_argument("-K", "--segments")
        .help("Number of segments, only used when loss function is '*_pla' or '*_ppla'")
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

    std::vector<std::vector<float>> weights(variant_matrix.size(), std::vector<float>(variant_matrix[0].size(), 1.0));
    if (!program.get<std::string>("-w").empty()) {
        weights = parse_matrix<float>(program.get<std::string>("-w"));
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
    int ncols = variant_matrix[0].size();
    auto [clone_tree, vertex_map] = parse_adjacency_list(program.get<std::string>("tree"));
    std::vector<int> vertex_map_vector(vertex_map.size());
    for (const auto& [key, value] : vertex_map) {
        if (key >= ncols) {
            error_logger->error("The nodes must range from 0, 1, ..., {} as {} is the number of columns in the input matrix.", ncols - 1, ncols);
            error_logger->error("The node {} is out of bounds.", key);
            std::exit(1);
        }

        vertex_map_vector[key] = value;
    }


    size_t tree_root;
    bool seen_root = false;
    for (auto n : clone_tree.nodes()) {
        if (clone_tree.in_degree(n) == 0) {
            if (seen_root) {
                throw std::invalid_argument("multiple root nodes found in the adjacency list");
            }

            tree_root = clone_tree[n].data;
            seen_root = true;
        }
    }

    if (!seen_root) {
        throw std::invalid_argument("no root node found in the adjacency list");
    }

    size_t nr_segments = program.get<int>("-K");
    double precision = program.get<double>("-s");

    SolverResult result;
    if (program.get<std::string>("-l") == "binomial") {
        result = log_binomial_admm_solve(vertex_map_vector, variant_matrix, total_matrix, clone_tree, tree_root);
    } else if (program.get<std::string>("-l") == "binomial_pla") {
        result = log_binomial_fixed_solve(vertex_map_vector, variant_matrix, total_matrix, clone_tree, tree_root, nr_segments);
    } else if (program.get<std::string>("-l") == "binomial_ppla") {
        result = log_binomial_solve(vertex_map_vector, variant_matrix, total_matrix, clone_tree, tree_root, nr_segments);
    } else if (program.get<std::string>("-l") == "beta_binomial_pla") {
        result = log_beta_binomial_solve(vertex_map_vector, variant_matrix, total_matrix, clone_tree, tree_root, nr_segments, precision);
    } else if (program.get<std::string>("-l") == "beta_binomial_ppla") {
        result = log_beta_binomial_fixed_solve(vertex_map_vector, variant_matrix, total_matrix, clone_tree, tree_root, nr_segments, precision);
    } else if (program.get<std::string>("-l") == "l2") {
        result = l2_solve(vertex_map_vector, variant_matrix, total_matrix, weights, clone_tree, tree_root);
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
        error_logger->error("Failed to serialize JSON document.");
    }

    yyjson_mut_doc_free(doc);
    return 0;
}
