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
struct decoding_algorithm_params 
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
    double primary{};       // Scaling factor applied in NMSA (α), OMSA (β), ANMSA (α) and AOMSA (β) decoding algorithms.
    double secondary{};     // Scaling factor applied in ANMSA (𝜈) and AOMSA (𝜍) decoding algorithms.
};

// Structure that stores code rate value that correspond to a range of QBER values from 'QBER_begin' to 'QBER_end' in 'QBER_step' increments.
// If 'QBER_begin'=='QBER_end', only one value is used and 'step' is not taken into account.
struct R_QBER_range
{
    double code_rate{};
    double QBER_begin{};
    double QBER_end{};
    double QBER_step{};
};

// Structure that stores code rate value that correspond to a range of δ/f_EC values from 'delta_begin/efficiency_begin' 
// to 'delta_end/efficiency_end' in 'delta_step/efficiency_step' increments.
// If 'delta_begin/efficiency_begin'=='delta_end/efficiency_end', only one value is used and 'step' is not taken into account.
// https://arxiv.org/pdf/1007.1616
struct R_adaptation_parameters_range
{
    double code_rate{};
    double delta_begin{};
    double delta_end{};
    double delta_step{};
    double efficiency_begin{};
    double efficiency_end{};
    double efficiency_step{};
};

// Ranges of values for δ and f_EC for a specific code rate, calculated based on the structure 'R_adaptation_parameters_range'.
struct R_adaptation_parameters_values
{
    std::vector<double> delta{};
    std::vector<double> efficiency{};
};

// Structure for 'R_QBER_adaptation_parameters_map' that stores δ and f_EC values that corresponds to QBER value.
struct QBER_adaptation_parameters
{
    double QBER{};
    double delta{};
    double efficiency{};
};

// Structure that stores QBER with code rate adaptation parameters that corresponds to code rate(R) value from config file.
struct R_QBER_adaptation_parameters_map
{
    double code_rate{};
    QBER_adaptation_parameters QBER_adapt_params{};
};

struct config_data
{
    // Number of threads for parallelizing runs.
    size_t THREADS_NUMBER{};

    // Number of runs with one combination.
    size_t TRIALS_NUMBER{};

    // Seed of simulation.
    size_t SIMULATION_SEED{};

    // Enables privacy maintenance after protocol execution.
    bool ENABLE_PRIVACY_MAINTENANCE{}; 

    // Measurement of protocol throughput (T). As the ratio of the number of bits remaining after protocol execution to the protocol run time (bits/s). 
    // It is recommended to perform experiments in single-threaded mode.
    bool ENABLE_THROUGHPUT_MEASUREMENT{};

    // Take RTT into account when calculating protocol throughput.
    bool CONSIDER_RTT{};

    // RTT (Round-Trip Time) in milliseconds.
    double RTT{};

    // Options:
    // 0.   SPA (Sum-Product Algorithm).
    // 1.   SPA with linear approximation of tanh and atanh functions.
    // 2.   NMSA (Normalized Min-Sum Algorithm) with α-factor.
    // 3.   OMSA (Offset Min-Sum Algorithm) with β-factor.
    // 4.   ANMSA (Adaptive Normalized Min-Sum Algorithm) with α-factor (0 < α =< 1) and 𝜈-factor (0 < 𝜈 =< 1).
    // 5.   AOMSA (Adaptive Offset Min-Sum Algorithm) with β-factor (0 < β =< 1) and 𝜍-factor (0 < 𝜍).
    size_t DECODING_ALGORITHM{};

    // Parameters specifying scaling factors for parity-check matrices when using NMSA, OMSA, ANMSA and AOMSA algorithms.
    decoding_algorithm_params DECODING_ALG_PARAMS{};

    // The maximum number of iterations of the decoding algorithm.
    size_t DECODING_ALG_MAX_ITERATIONS{};

    // Four options:
    //    0) Uncompressed matrices (folder matrices_uncompressed).
    //    1) Matrices in 'alist' format (folder matrices_alist).
    //        About 'alist' format: https://rptu.de/channel-codes/matrix-file-formats.
    //    2) Matrices in format specified below (folder matrices_1).
    //        The first line contains the block length, N. The second line defines the
    //        number of parity-checks, M. The third line defines the number of columns
    //        of the compressed parity-check matrix. The following M lines are then the
    //        compressed parity-check matrix. Each of the M rows contains the indices
    //        (1 ... N) of 1's in the compressed row of parity-check matrix. If not all
    //        column entries are used, the column is filled up with 0's. Program for
    //        matrix generation: https://www.inference.org.uk/mackay/PEG_ECC.html.
    //    3) Matrices in format specified below (folder matrices_2).
    //        The first line contains two numbers: the first is the block length (N)
    //        and the second is the number of parity-checks (M). The following M 
    //        lines are then the compressed parity-check matrix. Each of the M rows
    //        contains the indices (0 ... N-1) of 1's in the compressed row of
    //        parity-check matrix. The next N lines contains the indices (0 ... M-1)
    //        of 1's in the compressed column of parity-check matrix.
    size_t MATRIX_FORMAT{};

    // Output intermediate results of LDPC operation to the console.
    bool TRACE_QKD_LDPC{};

    // Output intermediate results of the decoding algorithm to the console.
    bool TRACE_DECODING_ALG{};

    // Console output of maximum log likelihood ratios (LLR) values of the message during the decoding algorithm.
    bool TRACE_DECODING_ALG_LLR{};

    // Enables limitation on the maximum LLR value of the message.
    // If the set threshold is exceeded, the value is set equal to DECODING_ALG_MSG_LLR_THRESHOLD.
    bool ENABLE_DECODING_ALG_MSG_LLR_THRESHOLD{};

    // The LLR value to which the message value is limited.
    double DECODING_ALG_MSG_LLR_THRESHOLD{};

    // Code rate and QBER correspondence set.
    std::vector<R_QBER_range> R_QBER_RANGES{};

    // Enables code rate modulation of pre-built codes by applying puncturing and shortening (https://arxiv.org/abs/1007.1616).
    bool ENABLE_CODE_RATE_ADAPTATION{};

    // Enables the use of the untainted puncturing when determining the positions of punctured bits (https://arxiv.org/pdf/1103.6149).
    bool ENABLE_UNTAINTED_PUNCTURING{};

    // If 'true', then 'R_ADAPT_PARAMS_RANGES' is used, otherwise 'R_QBER_ADAPT_PARAMS_MAPS'.
    bool USE_ADAPTATION_PARAMETERS_RANGES{};

    // Code rate correspondence set of δ and f_EC ranges for rate modulation.
    std::vector<R_adaptation_parameters_range> R_ADAPT_PARAMS_RANGES{};

    // Code rate correspondence set of QBER (instead of QBER values from 'R_QBER_RANGES') with δ and f_EC maps for rate modulation. 
    std::vector<R_QBER_adaptation_parameters_map> R_QBER_ADAPT_PARAMS_MAPS{};
};

extern config_data CFG;
const double EPSILON = 1e-6;

inline constexpr size_t DEC_SPA = 0, DEC_SPA_APPROX = 1, DEC_NMSA = 2, DEC_OMSA = 3, DEC_ANMSA = 4, DEC_AOMSA = 5;
inline constexpr size_t MAT_SPARSE_UNCOMPRESSED = 0, MAT_SPARSE_ALIST = 1, MAT_SPARSE_1 = 2, MAT_SPARSE_2 = 3;

scaling_factor_range parse_scaling_factor_range(const json& scaling_factor_range);

std::vector<R_scaling_factor_map> parse_scaling_factor_maps(
    const json& scaling_factor_maps,
    const std::string& key
);

void print_config_info(
    config_data cfg, 
    std::string cfg_name,
    size_t cfg_number
);

config_data parse_config_data(fs::path config_path);