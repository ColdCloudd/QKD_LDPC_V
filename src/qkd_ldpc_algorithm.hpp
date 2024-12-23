#pragma once

#include <fmt/core.h>
#include <fmt/color.h>
#include <fmt/ranges.h>

#include "config.hpp"
#include "utils.hpp"
#include "array_and_matrix_operations.hpp"

namespace fs = std::filesystem;

// Result of SPA or MSA algorithm
struct decoding_result
{
    size_t iterations_num{};
    bool syndromes_match{};
};

struct LDPC_result
{
    decoding_result decoding_res{};
    bool keys_match{};
};

decoding_result sum_product_decoding(const std::vector<double> &bit_array_llr,
                               const H_matrix &matrix,
                               const std::vector<int> &syndrome,
                               const size_t &max_num_iterations, 
                               const double &msg_threshold, 
                               std::vector<int> &bit_array_out);

decoding_result min_sum_decoding(const std::vector<double> &bit_array_llr,
                           const H_matrix &matrix,
                           const std::vector<int> &syndrome,          
                           const size_t &max_num_iterations, 
                           const double &msg_threshold, 
                           std::vector<int> &bit_array_out);

LDPC_result QKD_LDPC(const std::vector<int> &alice_bit_array, 
                     const std::vector<int> &bob_bit_array, 
                     const double &QBER, 
                     const H_matrix &matrix);