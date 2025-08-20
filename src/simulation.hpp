#pragma once
#include <cmath>
#include <string>
#include <vector>
#include <chrono>
#include <limits>
#include <locale>
#include <optional>
#include <algorithm>

#include <fmt/core.h>
#include <fmt/color.h>
#include <fmt/ranges.h>
#include <XoshiroCpp.hpp>
#include <BS_thread_pool.hpp>
#include <indicators/progress_bar.hpp>
#include <indicators/cursor_control.hpp>

#include "config.hpp"
#include "qkd_ldpc_algorithm.hpp"

struct sim_combination 
{
    double QBER;                
    H_matrix_params matrix_params;                  // Matrix parameters after modulation.
    decoding_scaling_factors scaling_factors;       // Scaling factors for NMSA, OMSA, ANMSA and AOMSA algorithms. 
};

struct sim_input
{
    H_matrix matrix{};
    fs::path matrix_path{};   
    std::vector<sim_combination> combinations;
};

struct trial_result
{
    LDPC_result ldpc_res{};
    double accurate_QBER{};
    std::chrono::microseconds runtime{};
};

struct sim_result
{
    size_t sim_number{};
    std::string matrix_filename{};
    bool is_regular{};                              // Matrix type.
    size_t num_bit_nodes{};                         // Number of bit nodes, which is defined as the number of columns in the parity check matrix.
    size_t num_check_nodes{};                       // Number of check nodes, which is defined as the number of rows in the parity check matrix.
    double delta{};                                 // Code rate modulation parameter. Fraction of punctured (p) and shortened (s) symbols (δ = π + σ).
    double efficiency{};                            // Code rate modulation parameter. Efficiency of the reconciliation (f_EC).
    double punctured_fraction{};                    // Fraction of punctured symbols (π) in modulated code.
    double shortened_fraction{};                    // Fraction of shortened symbols (σ) in modulated code.
    double adapted_code_rate{};                     // Modulated code rate (R).
    double accurate_QBER{};                         // An accurate QBER that corresponds to the number of errors in the key.
    decoding_scaling_factors scaling_factors{};     // Scaling factor(s) applied in NMSA, OMSA (or ANMSA, AOMSA) decoding algorithm.
    size_t iter_success_dec_alg_max{};              // The maximum number of iterations of the decoding algorithm in which Alice's syndrome matched Bob's syndrome (i.e. successful).
    size_t iter_success_dec_alg_min{};              // The minimum number of iterations of the decoding algorithm.
    double iter_success_dec_alg_mean{};             // The mean number of iterations of the decoding algorithm. 
    double iter_success_dec_alg_std_dev{};          // The standard deviation of iterations of the decoding algorithm. 
    double ratio_trials_success_dec_alg{};          // Success rate of the decoding algorithm. Success when Bob's syndrome matches Alice's.
    double ratio_trials_success_ldpc{};             // Success rate of the QKD LDPC error reconciliation. Success when Bob and Alice's keys match.
    size_t throughput_max{};                        // Throughput is measured as the ratio of the number of bits remaining after protocol execution to the protocol run time (bits/s).
    size_t throughput_min{};
    size_t throughput_mean{};
    size_t throughput_std_dev{};
};

fs::path write_file(
    const std::vector<sim_result> &data,
    std::string sim_duration,
    fs::path directory
);

std::vector<double> get_rate_based_QBER_range(
    const double code_rate,
    const std::vector<R_QBER_map> &R_QBER_maps
);

std::vector<double> get_scaling_factor_range_values(const scaling_factor_range &scaling_factor_range);

double get_rate_based_scaling_factor_value(
    const double code_rate,
    const std::vector<R_scaling_factor_map> &R_scaling_factor_maps
);

std::vector<sim_input> prepare_sim_inputs(const std::vector<fs::path> &matrix_paths);

trial_result run_trial(
    const H_matrix &matrix, 
    double QBER, 
    size_t seed,
    const H_matrix_params &matrix_params = {},
    const decoding_scaling_factors &scaling_factors = {}
);

void process_trials_results(
    const std::vector<trial_result> &trial_results, 
    const H_matrix &matrix,
    const H_matrix_params &matrix_params,
    sim_result &result
);

std::vector<sim_result> QKD_LDPC_batch_simulation(const std::vector<sim_input> &sim_in);
