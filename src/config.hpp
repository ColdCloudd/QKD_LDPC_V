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

// Alpha (or beta, nu, sigma) range specified using 'begin', 'end' and 'step'. If 'begin'=='end', only one value is used and 'step' is not taken into account.
struct scaling_factor_range
{
    double begin{};
    double end{};
    double step{};
};

// Structure that stores code rate value that correspond to alpha (or beta, nu, sigma) value.
struct R_scaling_factor_map
{
    double code_rate{};
    double scaling_factor{};
};

// Parameters for parity-check matrices when using NMSA, OMSA, ANMSA and AOMSA algorithms.
struct algorithm_params 
{
    // Common parameters for NMSA, OMSA, ANMSA and AOMSA algorithms.
    struct 
    {
        bool use_range{};                           // For all matrices, scaling factor values will be generated based on the range specified by 'range'.
        scaling_factor_range range{};               // Range for all matrices, specified using 'begin', 'end' and 'step'.
        std::vector<R_scaling_factor_map> maps{};   // Code rate and scaling factor correspondence set.
    } primary;
    
    // Additional parameters for ANMSA and AOMSA algorithms.
    struct 
    {
        bool use_range{};
        scaling_factor_range range{};
        std::vector<R_scaling_factor_map> maps{};
    } secondary;
};

// Structure for storing factors passed to NMSA, OMSA, ANMSA and AOMSA algorithms.
struct decoding_scaling_factors
{
    double primary{};
    double secondary{};
};

// Structure that stores code rate value that correspond to a range of QBER values from 'QBER_begin' to 'QBER_end' in 'QBER_step' increments.
// If 'QBER_begin'=='QBER_end', only one value is used and 'step' is not taken into account.
struct R_QBER_map
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

    // Four options:
    // 0.   SPA (Sum-Product Algorithm).
    // 1.   SPA with linear approximation of tanh and atanh functions.
    // 2.   NMSA (Normalized Min-Sum Algorithm) with α-factor.
    // 3.   OMSA (Offset Min-Sum Algorithm) with β-factor.
    // 4.   ANMSA (Adaptive Normalized Min-Sum Algorithm) with α-factor (0 < α =< 1) and 𝜈-factor (0 < 𝜈 =< 1).
    // 5.   AOMSA (Adaptive Offset Min-Sum Algorithm) with β-factor (0 < β =< 1) and 𝜍-factor (0 < 𝜍).
    size_t DECODING_ALGORITHM{};

    // Parameters specifying scaling factors for parity-check matrices when using NMSA, OMSA, ANMSA and AOMSA algorithms.
    algorithm_params DECODING_ALG_PARAMS{};

    // The maximum number of iterations of the decoding algorithm.
    // If the maximum number of iterations is reached, error reconciliation in the key is considered unsuccessful.
    size_t DECODING_ALG_MAX_ITERATIONS{};

    // Four options:
    // 0.   Dense matrices (folder dense_matrices).
    // 1.   Sparse matrices in '.alist' format (folder sparse_matrices_alist).
    // 2.   Sparse matrices in format specified below (folder sparse_matrices_1).
    // The first line contains the block length, N. The second line defines the number of parity-checks, M.
    // The third line defines the number of columns of the compressed parity-check matrix. 
    // The following M lines are then the compressed parity-check matrix. Each of the M rows contains the 
    // indices (1 ... N) of 1's in the compressed row of parity-check matrix. If not all column entries are used, 
    // the column is filled up with 0's. Program for matrix generation: https://www.inference.org.uk/mackay/PEG_ECC.html.
    // 3.   Sparse matrices in format specified below (folder sparse_matrices_2).
    // Read sparse matrix from file in format:
    // The first line contains two numbers: the first is the block length (N) and the second is the number of parity-checks (M).
    // The following M lines are then the compressed parity-check matrix. Each of the M rows contains the 
    // indices (0 ... N-1) of 1's in the compressed row of parity-check matrix. 
    // The next N lines contains the indices (0 ... M-1) of 1's in the compressed column of parity-check matrix.
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
    std::vector<R_QBER_map> R_QBER_MAPS{};
};

extern config_data CFG;
const double EPSILON = 1e-6;

inline constexpr size_t SPA = 0, SPA_APPROX = 1, NMSA = 2, OMSA = 3, ANMSA = 4, AOMSA = 5;
inline constexpr size_t DENSE_MAT = 0, SPARSE_MAT_ALIST = 1, SPARSE_MAT_1 = 2, SPARSE_MAT_2 = 3;

scaling_factor_range parse_scaling_factor_range(const json& scaling_factor_range);

std::vector<R_scaling_factor_map> parse_scaling_factor_maps(const json& scaling_factor_maps, const std::string& key);

config_data parse_config_data(fs::path config_path);