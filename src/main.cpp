﻿#include "config.hpp"
#include "simulation.hpp"

namespace fs = std::filesystem;

#ifdef USE_CURRENT_DIR
const fs::path CONFIG_PATH = std::filesystem::current_path() / "config.json";
const fs::path DENSE_MATRIX_DIR_PATH = std::filesystem::current_path() / "dense_matrices";
const fs::path SPARSE_MATRIX_ALIST_DIR_PATH = std::filesystem::current_path() / "sparse_matrices_alist";
const fs::path SPARSE_MATRIX_1_DIR_PATH = std::filesystem::current_path() / "sparse_matrices_1";
const fs::path SPARSE_MATRIX_2_DIR_PATH = std::filesystem::current_path() / "sparse_matrices_2";
const fs::path RESULTS_DIR_PATH = std::filesystem::current_path() / "results";
#else
const fs::path CONFIG_PATH = fs::path(SOURCE_DIR) / "config.json";
const fs::path DENSE_MATRIX_DIR_PATH = fs::path(SOURCE_DIR) / "dense_matrices";
const fs::path SPARSE_MATRIX_ALIST_DIR_PATH = fs::path(SOURCE_DIR) / "sparse_matrices_alist";
const fs::path SPARSE_MATRIX_1_DIR_PATH = fs::path(SOURCE_DIR) / "sparse_matrices_1";
const fs::path SPARSE_MATRIX_2_DIR_PATH = fs::path(SOURCE_DIR) / "sparse_matrices_2";
const fs::path RESULTS_DIR_PATH = fs::path(SOURCE_DIR) / "results";
#endif

config_data CFG;

int main()
{
    try
    {
        CFG = parse_config_data(CONFIG_PATH);
        fs::path matrix_dir_path;
        if (CFG.MATRIX_FORMAT == DENSE_MAT)
            matrix_dir_path = DENSE_MATRIX_DIR_PATH;
        else if (CFG.MATRIX_FORMAT == SPARSE_MAT_ALIST)
            matrix_dir_path = SPARSE_MATRIX_ALIST_DIR_PATH;
        else if (CFG.MATRIX_FORMAT == SPARSE_MAT_1)
            matrix_dir_path = SPARSE_MATRIX_1_DIR_PATH;
        else if (CFG.MATRIX_FORMAT == SPARSE_MAT_2)
            matrix_dir_path = SPARSE_MATRIX_2_DIR_PATH;
        
        if (CFG.INTERACTIVE_MODE)
        {
            fmt::print(fg(fmt::color::purple), "INTERACTIVE MODE\n");
            QKD_LDPC_interactive_simulation(matrix_dir_path);
        }
        else
        {
            fmt::print(fg(fmt::color::purple), "BATCH MODE\n");
            
            std::vector<fs::path> matrix_paths = get_file_paths_in_directory(matrix_dir_path);
            if (matrix_paths.empty())
            {
                throw std::runtime_error("Matrix folder is empty: " + matrix_dir_path.string());
            }
            
            std::vector<sim_input> sim_inputs = prepare_sim_inputs(matrix_paths);
            std::vector<sim_result> sim_results = QKD_LDPC_batch_simulation(sim_inputs);

            fmt::print(fg(fmt::color::green), "The results will be written to the directory: {}\n", RESULTS_DIR_PATH.string());
            write_file(sim_results, RESULTS_DIR_PATH);
        }
    }
    catch (const std::exception &e)
    {
        fmt::print(stderr, fg(fmt::color::red), "ERROR: {}\n", e.what());
        fmt::print(fg(fmt::color::red), "Press Enter to exit...");
        std::cin.get();
        return EXIT_FAILURE;
    }

    fmt::print(fg(fmt::color::green), "Press Enter to exit...");
    std::cin.get();

    return EXIT_SUCCESS;
}
