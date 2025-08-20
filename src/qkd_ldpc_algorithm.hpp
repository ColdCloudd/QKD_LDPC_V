#pragma once

#include <fmt/core.h>
#include <fmt/color.h>
#include <fmt/ranges.h>

#include "config.hpp"
#include "utils.hpp"
#include "array_and_matrix_operations.hpp"

namespace fs = std::filesystem;

const double ALMOST_ZERO = 1e-4;

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

decoding_result sum_product_decoding(
    const std::vector<double> &bit_array_llr,
    const H_matrix &matrix,
    const std::vector<int> &syndrome,          
    const size_t &max_num_iterations, 
    const double &msg_threshold, 
    std::vector<int> &bit_array_out
);

double tanh_lin_approx(double x);

double atanh_lin_approx(double x);

decoding_result sum_product_linear_approx_decoding(
    const std::vector<double> &bit_array_llr,
    const H_matrix &matrix,
    const std::vector<int> &syndrome,          
    const size_t &max_num_iterations, 
    const double &msg_threshold, 
    std::vector<int> &bit_array_out
);

decoding_result min_sum_normalized_decoding(
    const std::vector<double> &bit_array_llr,
    const H_matrix &matrix, 
    const std::vector<int> &syndrome,              
    const size_t &max_num_iterations,
    const double &alpha,   
    const double &msg_threshold,      
    std::vector<int> &bit_array_out
);

decoding_result min_sum_offset_decoding(
    const std::vector<double> &bit_array_llr,
    const H_matrix &matrix, 
    const std::vector<int> &syndrome,              
    const size_t &max_num_iterations,
    const double &beta,   
    const double &msg_threshold,      
    std::vector<int> &bit_array_out
);

decoding_result adaptive_min_sum_normalized_decoding(
    const std::vector<double> &bit_array_llr,
    const H_matrix &matrix, 
    const std::vector<int> &syndrome,              
    const size_t &max_num_iterations,
    const double &alpha, 
    const double &nu, 
    const double &msg_threshold,      
    std::vector<int> &bit_array_out
);

decoding_result adaptive_min_sum_offset_decoding(
    const std::vector<double> &bit_array_llr,
    const H_matrix &matrix, 
    const std::vector<int> &syndrome,              
    const size_t &max_num_iterations,
    const double &beta, 
    const double &sigma, 
    const double &msg_threshold,      
    std::vector<int> &bit_array_out
);

LDPC_result QKD_LDPC(
    const H_matrix &matrix,
    const std::vector<int> &alice_bit_array, 
    const std::vector<int> &bob_bit_array, 
    const double &QBER,
    const decoding_scaling_factors &scaling_factors = {},
    const H_matrix_params &matrix_params = {}
);

LDPC_result QKD_LDPC_RATE_ADAPT(
    const H_matrix &matrix,
    const std::vector<int> &alice_bit_array, 
    const std::vector<int> &bob_bit_array, 
    const double &QBER,
    const decoding_scaling_factors &scaling_factors,
    const H_matrix_params &matrix_params,
    XoshiroCpp::Xoshiro256PlusPlus &prng
);