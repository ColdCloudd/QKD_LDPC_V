#include "config.hpp"
#include "simulation.hpp"

namespace fs = std::filesystem;

#ifdef USE_CURRENT_DIR
const fs::path CONFIG_DIR_PATH = std::filesystem::current_path() / "configs";
const fs::path UNCOMPRESSED_MATRIX_DIR_PATH = std::filesystem::current_path() / "matrices" / "uncompressed_matrices";
const fs::path SPARSE_MATRIX_ALIST_DIR_PATH = std::filesystem::current_path() / "matrices" / "sparse_matrices_alist";
const fs::path SPARSE_MATRIX_1_DIR_PATH = std::filesystem::current_path() / "matrices" / "sparse_matrices_1";
const fs::path SPARSE_MATRIX_2_DIR_PATH = std::filesystem::current_path() / "matrices" / "sparse_matrices_2";
const fs::path RESULTS_DIR_PATH = std::filesystem::current_path() / "results";
#else
const fs::path CONFIG_DIR_PATH = fs::path(SOURCE_DIR) / "configs";
const fs::path UNCOMPRESSED_MATRIX_DIR_PATH = fs::path(SOURCE_DIR) / "matrices" / "uncompressed_matrices";
const fs::path SPARSE_MATRIX_ALIST_DIR_PATH = fs::path(SOURCE_DIR) / "matrices" / "sparse_matrices_alist";
const fs::path SPARSE_MATRIX_1_DIR_PATH = fs::path(SOURCE_DIR) / "matrices" / "sparse_matrices_1";
const fs::path SPARSE_MATRIX_2_DIR_PATH = fs::path(SOURCE_DIR) / "matrices" / "sparse_matrices_2";
const fs::path RESULTS_DIR_PATH = fs::path(SOURCE_DIR) / "results";
#endif

config_data CFG;

int main(int argc, char* argv[])
{
    try
    {
        for (int i = 1; i < argc; ++i) 
        {
            std::string arg = argv[i];
            if (arg == "-help" || arg == "--help") 
            {
                std::string help = 
                "-----------------------------CONFIGURATION PARAMETERS-----------------------------\n"
                
                "threads_number - Number of threads for parallelizing runs.\n\n"

                "trials_number - Number of runs with one combination of parameters.\n\n"

                "use_config_simulation_seed - If true, the seed specified in configuration is used,\n"
                "   otherwise a random seed is used.\n\n"

                "simulation_seed - Initial value of the PRNG.\n\n"

                "enable_privacy_maintenance - Enables privacy maintenance after protocol execution.\n\n"

                "enable_throughput_measurement - Measurement of protocol throughput (T, bits/s).\n" 
                "   As the ratio of the number of bits remaining after protocol execution to the\n"
                "   protocol run time. Recommended to perform experiments in single-threaded mode.\n\n"

                "throughput_measurement_parameters - Parameters of throughput (T) measurement.\n"
                "   consider_RTT - Include RTT value in throughput calculations.\n"
                "   RTT - Round Trip Time in milliseconds (ms).\n\n"

                "decoding_algorithm - LDPC decoding algorithms. Options:\n"
                "   0) SPA (Sum-Product Algorithm).\n"
                "   1) SPA with linear approximation of tanh and atanh functions.\n"
                "   2) NMSA (Normalized Min-Sum Algorithm) with α-factor. More info:\n"
                "      http://ieeexplore.ieee.org/document/1495850.\n"
                "   3) OMSA (Offset Min-Sum Algorithm) with β-factor. More info:\n"
                "      http://ieeexplore.ieee.org/document/1495850.\n"
                "   4) ANMSA (Adaptive Normalized Min-Sum Algorithm) with α-factor (0 < α =< 1) and\n"
                "      𝜈-factor (0 < 𝜈 =< 1). More info: https://ieeexplore.ieee.org/document/5545625.\n"
                "   5) AOMSA (Adaptive Offset Min-Sum Algorithm) with β-factor (0 < β =< 1) and\n" 
                "      𝜍-factor (0 < 𝜍). More info: https://ieeexplore.ieee.org/document/5545625.\n\n"

                "min_sum_normalized(offset)_parameters - Parameters specifying scaling factors\n"
                "   for parity-check matrices (decoding algorithms NMSA and OMSA).\n"
                "   use_alpha(beta)_range - For all matrices, scaling factor values will be\n"
                "       generated based on the range specified by 'alpha(beta)_range'.\n"
                "   alpha(beta)_range - Range for all matrices, specified using 'begin', 'end'\n"
                "       and 'step'.\n"
                "   code_rate_alpha(beta)_maps - Code rate(R) of parity-check matrix and scaling\n"
                "       factor (α or β) correspondence set.\n\n"

                "adaptive_min_sum_normalized(offset)_parameters - Parameters specifying scaling \n"
                "   factors for parity-check matrices (decoding algorithms ANMSA and AOMSA).\n"
                "   Similar to 'min_sum_normalized(offset)_parameters', but with one additional\n"
                "   scaling factor (𝜈 or 𝜍) specified.\n\n"

                "decoding_algorithm_max_iterations - The maximum number of iterations of the\n"
                "   decoding algorithm.\n\n"

                "matrix_format - Matrix representation format in read files.Options:\n"
                "   0) Uncompressed matrices (folder uncompressed_matrices).\n"
                "   1) Sparse matrices in 'alist' format (folder sparse_matrices_alist).\n" 
                "       About 'alist' format: https://rptu.de/channel-codes/matrix-file-formats.\n"
                "   2) Sparse matrices in format specified below (folder sparse_matrices_1).\n"
                "       The first line contains the block length, N. The second line defines the\n"
                "       number of parity-checks, M. The third line defines the number of columns\n"
                "       of the compressed parity-check matrix. The following M lines are then the\n"
                "       compressed parity-check matrix. Each of the M rows contains the indices\n"
                "       (1 ... N) of 1's in the compressed row of parity-check matrix. If not all\n"
                "       column entries are used, the column is filled up with 0's. Program for\n"
                "       matrix generation: https://www.inference.org.uk/mackay/PEG_ECC.html.\n"
                "   3) Sparse matrices in format specified below (folder sparse_matrices_2).\n"
                "       The first line contains two numbers: the first is the block length (N)\n"
                "       and the second is the number of parity-checks (M). The following M \n"
                "       lines are then the compressed parity-check matrix. Each of the M rows\n"
                "       contains the indices (0 ... N-1) of 1's in the compressed row of\n"
                "       parity-check matrix. The next N lines contains the indices (0 ... M-1)\n"
                "       of 1's in the compressed column of parity-check matrix.\n\n"

                "trace_qkd_ldpc - Output to the console of input and output data of the decoding\n"
                "   algorithm.\n\n"

                "trace_decoding_algorithm - Output to the console intermediate results of the\n"
                "   decoding algorithm.\n\n"

                "trace_decoding_algorithm_llr - Output to the console of maximum log likelihood\n"
                "   ratios (LLR) values of the message in Tanner graph during the decoding algorithm.\n\n"

                "enable_decoding_algorithm_msg_llr_threshold - Enables limitation on the maximum LLR\n"
                "   value of the message. If the set threshold is exceeded, the value is set equal\n"
                "   to 'decoding_algorithm_msg_llr_threshold'.\n\n"

                "decoding_algorithm_msg_llr_threshold - The LLR value to which the message value is\n"
                "   limited.\n\n"

                "code_rate_QBER_maps - Code rate(R) and QBER correspondence set. The 'code_rate'\n"
                "   value of parity-check matrix correspond to a range of QBER values from\n"
                "   'begin' to 'end' in 'step' increments. If 'begin' == 'end', only one value\n"
                "   will be used, and 'step' will not be taken into account.\n\n"
                
                "enable_code_rate_adaptation - Enables code rate(R) modulation using puncturing\n"
                "   and shortening techniques. More info: https://arxiv.org/abs/1007.1616.\n\n"

                "enable_untainted_puncturing - Enables (when code rate(R) modulation enabled) the use\n"
                "   of bits obtained using the algorithm for puncturing untainted nodes as punctured\n"
                "   symbols. Otherwise, random positions of punctured nodes are used. More info:\n"
                "   https://ieeexplore.ieee.org/document/6290312.\n\n"

                "code_rate_adaptation_parameters_maps - Code rate(R) and delta(δ) with reconciliation\n" 
                "   efficiency(f_EC) correspondence set. The 'code_rate' value of parity-check matrix\n"
                "   correspond to a range of δ and f_EC values from 'begin' to 'end' in 'step'\n"
                "   increments. If 'begin' == 'end', only one value will be used, and 'step' will\n"
                "   not be taken into account.\n"

                "----------------------------------------------------------------------------------\n\n";
                
                fmt::print(fg(fmt::color::honey_dew), "{}", help);
                fmt::print(fg(fmt::color::honey_dew), "Press Enter to exit...");
                std::cin.get();
                return EXIT_SUCCESS;
            }
        }

        std::vector<fs::path> config_paths = get_file_paths_in_directory(CONFIG_DIR_PATH, ".json");
        
        for (size_t i = 0; i < config_paths.size(); ++i)
        {
            CFG = parse_config_data(config_paths[i]);
            print_config_info(CFG, config_paths[i].filename().string(), i + 1);
            fs::path matrix_dir_path;
            if (CFG.MATRIX_FORMAT == MAT_UNCOMPRESSED)
                matrix_dir_path = UNCOMPRESSED_MATRIX_DIR_PATH;
            else if (CFG.MATRIX_FORMAT == MAT_SPARSE_ALIST)
                matrix_dir_path = SPARSE_MATRIX_ALIST_DIR_PATH;
            else if (CFG.MATRIX_FORMAT == MAT_SPARSE_1)
                matrix_dir_path = SPARSE_MATRIX_1_DIR_PATH;
            else if (CFG.MATRIX_FORMAT == MAT_SPARSE_2)
                matrix_dir_path = SPARSE_MATRIX_2_DIR_PATH;
            
            std::vector<fs::path> matrix_paths = get_file_paths_in_directory(matrix_dir_path, ".mtrx");
            std::vector<sim_input> sim_inputs = prepare_sim_inputs(matrix_paths);

            auto sim_start = std::chrono::high_resolution_clock::now();
            std::vector<sim_result> sim_results = QKD_LDPC_batch_simulation(sim_inputs);
            auto sim_end = std::chrono::high_resolution_clock::now();

            auto hours = std::chrono::duration_cast<std::chrono::hours>(sim_end - sim_start);
            auto minutes = std::chrono::duration_cast<std::chrono::minutes>(sim_end - sim_start - hours);
            auto seconds = std::chrono::duration_cast<std::chrono::seconds>(sim_end - sim_start - hours - minutes);
            std::string sim_duration = fmt::format("{:02}h-{:02}m-{:02}s", hours.count(), minutes.count(), seconds.count());
            
            fs::path result_file_path = write_file(sim_results, sim_duration, RESULTS_DIR_PATH);
            fmt::print(fg(fmt::color::lawn_green), "The results are written to the file: {}{}\n\n", 
                (RESULTS_DIR_PATH / "").string(), fmt::styled(result_file_path.filename().string(), fg(fmt::color::crimson)));
            
        }   
    }
    catch (const std::exception &e)
    {
        fmt::print(stderr, fg(fmt::color::red), "ERROR: {}\n", e.what());
        fmt::print(fg(fmt::color::red), "Press Enter to exit...");
        std::cin.get();
        return EXIT_FAILURE;
    }

    fmt::print(fg(fmt::color::lawn_green), "Simulations successfully completed! Press Enter to exit...");
    std::cin.get();

    return EXIT_SUCCESS;
}
