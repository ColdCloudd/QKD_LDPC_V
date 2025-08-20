#include <filesystem>

#include "config.hpp"
#include "qkd_ldpc_algorithm.hpp"

namespace fs = std::filesystem;

const fs::path EXAMPLE_DENSE_MATRIX_PATH = fs::path(SOURCE_DIR) / "dense_matrices/(N=6,K=2,M=4,R=0.34).txt";

config_data CFG;

int main()
{
    try
    {
        CFG.DECODING_ALG_MAX_ITERATIONS = 100;

        CFG.DECODING_ALGORITHM = 0;               // Use Sum-Product
        CFG.ENABLE_DECODING_ALG_MSG_LLR_THRESHOLD = true;
        CFG.DECODING_ALG_MSG_LLR_THRESHOLD = 100.;

        CFG.TRACE_QKD_LDPC = true;
        CFG.TRACE_DECODING_ALG = true;
        CFG.TRACE_DECODING_ALG_LLR = true;
        
        H_matrix matrix = read_uncompressed_matrix(EXAMPLE_DENSE_MATRIX_PATH);

        // Page 33, example 2.5: https://www.researchgate.net/publication/228977165_Introducing_Low-Density_Parity-Check_Codes
        std::vector<int> alice_bit_array{0, 0, 1, 0, 1, 1};
        std::vector<int> bob_bit_array{1, 0, 1, 0, 1, 1};

        double QBER = 0.2;
        QKD_LDPC(matrix, alice_bit_array, bob_bit_array, QBER);
    }
    catch (const std::exception &e)
    {
        fmt::print(stderr, fg(fmt::color::red), "ERROR: {}\n", e.what());
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}