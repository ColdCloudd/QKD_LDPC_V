#pragma once
#include <vector>
#include <fstream>
#include <filesystem>

#include <fmt/core.h>
#include <fmt/color.h>
#include <fmt/ranges.h>
#include <nlohmann/json.hpp>

using json = nlohmann::json;
namespace fs = std::filesystem;

// Structure that stores code rate values that correspond to a range of QBER values from QBER_begin to QBER_end in QBER_step increments.
struct R_QBER_params
{
    double code_rate{};
    double QBER_begin{};
    double QBER_end{};
    double QBER_step{};
};

struct config_data
{
    // Number of threads for parallelizing runs.
    size_t THREADS_NUMBER{};

    // Number of runs with one combination.
    size_t TRIALS_NUMBER{};

    // Seed of simulation.
    size_t SIMULATION_SEED{};

    // If true, the user can select the parity check matrix and run tests with it. Otherwise, the BATCH MODE works, 
    // which reads all parity check matrices from the specified directory and runs tests with them and writes the results to a file.
    bool INTERACTIVE_MODE{};

    // Enables privacy maintenance after protocol execution.
    bool ENABLE_PRIVACY_MAINTENANCE{}; 

    // Measurement of protocol throughput (T). As the ratio of the number of bits remaining after protocol execution to the protocol run time (bits/s). 
    // It is recommended to perform experiments in single-threaded mode.
    bool ENABLE_THROUGHPUT_MEASUREMENT{};

    // Take RTT into account when calculating protocol throughput.
    bool CONSIDER_RTT{};

    // RTT (Round-Trip Time) in milliseconds.
    size_t RTT{};

    // Use MSA (Min-Sum Algorithm) decoding algorithm instead of SPA (Sum-Product Algorithm).
    bool USE_MIN_SUM_DECODING_ALG{};

    // The maximum number of iterations of the decoding algorithm.
    // If the maximum number of iterations is reached, error reconciliation in the key is considered unsuccessful.
    size_t DECODING_ALG_MAX_ITERATIONS{};

    // Use dense matrices (folder dense_matrices) instead of sparse matrices (folder alist_sparse_matrices).
    // Three options:
    // 0. Dense matrices (folder dense_matrices).
    // 1. Sparse matrices in .alist format (folder alist_sparse_matrices).
    // 2. Sparse matrices in format specified below (folder sparse_matrices).
    // The first line contains the block length, N. The second line defines the number of parity-checks, M.
    // The third line defines the number of columns of the compressed parity-check matrix. 
    // The following M lines are then the compressed parity-check matrix. Each of the M rows contains the 
    // indices (1 ... N) of 1's in the compressed row of parity-check matrix. If not all column entries are used, 
    // the column is filled up with 0's. Program for matrix generation: https://www.inference.org.uk/mackay/PEG_ECC.html
    size_t MATRIX_FORMAT{};

    // Output intermediate results of LDPC operation to the console.
    bool TRACE_QKD_LDPC{};

    // Output intermediate results of the decoding algorithm to the console.
    bool TRACE_DECODING_ALG{};

    // Console output of maximum log likelihood ratios (LLR) values of the message during the decoding algorithm.
    bool TRACE_DECODING_ALG_LLR{};

    // Enables limitation on the maximum LLR value of the message, which cannot exceed the set threshold. 
    // If the set threshold is exceeded, the value is set equal to DECODING_ALG_MSG_LLR_THRESHOLD.
    bool ENABLE_DECODING_ALG_MSG_LLR_THRESHOLD{};

    // The maximum LLR value a message can have.
    double DECODING_ALG_MSG_LLR_THRESHOLD{};

    // Code rate and QBER correspondence set.
    std::vector<R_QBER_params> R_QBER_PARAMETERS{};
};

extern config_data CFG;

config_data get_config_data(fs::path config_path);