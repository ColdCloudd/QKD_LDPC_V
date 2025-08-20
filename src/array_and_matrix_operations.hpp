#pragma once 
#include <set>
#include <random>
#include <vector>
#include <numeric>
#include <fstream>
#include <algorithm>
#include <filesystem>

#include <fmt/core.h>
#include <fmt/color.h>
#include <fmt/ranges.h>
#include <XoshiroCpp.hpp>

#include "config.hpp"

namespace fs = std::filesystem;

// Structure for storing information about second-order neighbors.
struct second_order_neighbors 
{
    int bit_node_idx;               // Index of bit node.
    std::set<int> neighbors;        // Set of second-order neighbors of this node.
};

// Structure for code rate modulation (https://arxiv.org/abs/1007.1616).
struct H_matrix_params
{
    // Fraction of punctured (p) and shortened (s) symbols (δ = π + σ).
    double delta{};

    // Efficiency of the reconciliation (f_EC).
    double efficiency{};

    // Fraction of punctured symbols (π).
    double punctured_fraction{};

    // Fraction of shortened symbols (σ).
    double shortened_fraction{};

    // Modulated code rate (R).
    double adapted_code_rate{};

    // Positions of punctured (p) bits.
    std::vector<int> punctured_bits{};

    // Positions of shortened (s) bits.
    std::vector<int> shortened_bits{};

    // A set of bits that must be removed.
    // Four cases:
    // 1. Code rate adaptation is disabled and privacy maintenance is disabled - array is empty.
    // 2. Code rate adaptation is disabled and privacy maintenance is enabled - array contains m bits positions (m - number of check nodes).
    // 3. Code rate adaptation is enabled and privacy maintenance is disabled - array contains only positions of punctured and shortened bits (p + s).
    // 4. Code rate adaptation is enabled and privacy maintenance is enabled - array contains s + m bits positions (p < m).
    std::vector<int> bits_to_remove{}; 
};

// Structure for storing parity-check matrix in compressed form.
struct H_matrix
{
    // A two-dimensional array that establishes a correspondence between the positions of bit nodes
    // and the positions of check nodes that control them in parity check matrix.
    std::vector<std::vector<int>> bit_nodes{};

    // A two-dimensional array that establishes a correspondence between the positions of check nodes
    // and the positions of bit nodes that controlled by them in parity check matrix.
    std::vector<std::vector<int>> check_nodes{};  

    // Puncturable bits (symbol nodes) indexes obtained using untainted puncturing algorithm
    // (https://arxiv.org/abs/1103.6149). Contains the maximum possible number of bits that 
    // can be punctured for the given matrix.
    std::vector<int> punctured_bits_untainted{};

    // Matrix type (regular or irregular).
    bool is_regular{};                  
};

std::vector<std::vector<int>> get_bit_nodes(
    const std::vector<std::vector<int8_t>> &matrix, 
    const std::vector<int> &bit_nodes_weight
);

std::vector<std::vector<int>> get_check_nodes(
    const std::vector<std::vector<int8_t>> &matrix, 
    const std::vector<int> &check_nodes_weight
);

std::vector<std::vector<int>> get_bit_nodes_from_check_nodes(const std::vector<std::vector<int>> &check_nodes);

double get_max_llr(const std::vector<std::vector<double>> &matrix);

bool arrays_equal(
    const std::vector<int> &array1,
    const std::vector<int> &array2
);

int find_available_index(
    const std::vector<int> &array1, 
    const std::vector<int> &array2
);

std::vector<int> get_bits_positions_to_remove(const H_matrix &matrix);

std::vector<int> get_bits_positions_to_remove_rate_adapt(
    const H_matrix &matrix, 
    const H_matrix_params &mat_params
);

void remove_bits(
    const std::vector<int> &bits_to_remove, 
    const std::vector<int> &array1, 
    const std::vector<int> &array2,
    std::vector<int> &array1_out,
    std::vector<int> &array2_out
);

H_matrix read_sparse_matrix_alist(const fs::path &matrix_path);

H_matrix read_sparse_matrix_1(const fs::path &matrix_path);

H_matrix read_sparse_matrix_2(const fs::path &matrix_path);

H_matrix read_uncompressed_matrix(const fs::path &matrix_path);

void fill_random_bits(
    XoshiroCpp::Xoshiro256PlusPlus &prng,
    std::vector<int> &bit_array
);

double inject_errors(
    XoshiroCpp::Xoshiro256PlusPlus &prng, 
    const std::vector<int> &bit_array, 
    double QBER, 
    std::vector<int> &bit_array_with_errors_out
);

void calculate_syndrome(
    const std::vector<int> &bit_array,
    const H_matrix &matrix, 
    std::vector<int> &syndrome_out
);
                        
void threshold_matrix(
    std::vector<std::vector<double>> &matrix, 
    const double &msg_threshold
);

std::vector<second_order_neighbors> get_second_order_neighbors(const H_matrix &matrix);

std::vector<int> select_punctured_bits_untainted(
    XoshiroCpp::Xoshiro256PlusPlus &prng,
    const H_matrix &matrix
);

std::vector<int> get_punctured_bits_untainted(
    const fs::path &matrix_path,
    XoshiroCpp::Xoshiro256PlusPlus &prng,
    const H_matrix &matrix
);

H_matrix_params adapt_code_rate(
    XoshiroCpp::Xoshiro256PlusPlus &prng, 
    const H_matrix &matrix,
    double QBER,
    double delta,
    double efficiency
);
