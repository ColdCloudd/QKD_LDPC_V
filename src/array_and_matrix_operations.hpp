#pragma once 
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

namespace fs = std::filesystem;

struct H_matrix
{
    std::vector<std::vector<int>> bit_nodes{};          // A two-dimensional array that establishes a correspondence between the positions of bit nodes and the positions of check nodes that control them in parity check matrix.
    std::vector<std::vector<int>> check_nodes{};        // A two-dimensional array that establishes a correspondence between the positions of check nodes and the positions of bit nodes that controlled by them in parity check matrix.
    bool is_regular{};                                  // Matrix type.
};

std::vector<std::vector<int>> get_bit_nodes(const std::vector<std::vector<int>> &matrix, 
                                            const std::vector<int> &bit_nodes_weight);
std::vector<std::vector<int>> get_check_nodes(const std::vector<std::vector<int>> &matrix, 
                                              const std::vector<int> &check_nodes_weight);
double get_max_llr(const std::vector<std::vector<double>> &matrix);
bool arrays_equal(const std::vector<int> &array1,
                  const std::vector<int> &array2);
H_matrix read_sparse_alist_matrix(const fs::path &matrix_path);
H_matrix read_dense_matrix(const fs::path &matrix_path);
void fill_random_bits(XoshiroCpp::Xoshiro256PlusPlus &prng,
                      std::vector<int> &bit_array);
double introduce_errors(XoshiroCpp::Xoshiro256PlusPlus &prng, 
                        const std::vector<int> &bit_array, 
                        double QBER, 
                        std::vector<int> &bit_array_with_errors_out);
void calculate_syndrome(const std::vector<int> &bit_array,
                        const H_matrix &matrix, 
                        std::vector<int> &syndrome_out);
void threshold_matrix(std::vector<std::vector<double>> &matrix, 
                      const double &msg_threshold);















