#include "qkd_ldpc_algorithm.hpp"

decoding_result sum_product_decoding(const std::vector<double> &bit_array_llr,
                                     const H_matrix &matrix,
                                     const std::vector<int> &syndrome,          
                                     const size_t &max_num_iterations, 
                                     const double &msg_threshold, 
                                     std::vector<int> &bit_array_out)
{
    double max_llr_c2b = 0.;
    double max_llr_b2c = 0.;
    double max_llr = 0.;

    auto& bit_nodes = matrix.bit_nodes;
    auto& check_nodes = matrix.check_nodes;
    size_t num_bit_nodes = bit_nodes.size();
    size_t num_check_nodes = check_nodes.size();

    std::vector<std::vector<double>> bit_to_check_msg(num_check_nodes);
    for (size_t i = 0; i < num_check_nodes; ++i)
    {
        bit_to_check_msg[i].resize(check_nodes[i].size());
        for (size_t j = 0; j < check_nodes[i].size(); ++j)
        {
            bit_to_check_msg[i][j] = bit_array_llr[check_nodes[i][j]]; // Initialization
        }
    }

    std::vector<std::vector<double>> check_to_bit_msg(num_bit_nodes);
    for (size_t i = 0; i < num_bit_nodes; ++i)
    {
        check_to_bit_msg[i].resize(bit_nodes[i].size());
    }

    double prod;
    double row_prod;
    int curr_bit_pos;

    std::vector<int> check_pos_idx(num_bit_nodes);
    std::vector<double> total_bit_llr(num_bit_nodes);
    std::vector<int> decision_syndrome(num_check_nodes);
    std::vector<int> bit_pos_idx(num_check_nodes);

    double sum;
    double col_sum;
    int curr_check_pos;

    size_t curr_iteration = 0;
    while (curr_iteration != max_num_iterations)
    {
        if (CFG.TRACE_DECODING_ALG)
        {
            fmt::print(fg(fmt::color::blue), "\n\nIteration: {}\n", curr_iteration + 1);
        }

        // Compute extrinsic messages from check nodes to bit nodes (Step 1: Check messages)
        for (size_t i = 0; i < num_check_nodes; ++i)
        {
            for (size_t j = 0; j < check_nodes[i].size(); ++j)
            {
                bit_to_check_msg[i][j] = tanh(bit_to_check_msg[i][j] / 2.);
            }
        }

        std::fill(check_pos_idx.begin(), check_pos_idx.end(), 0);
        for (size_t j = 0; j < num_check_nodes; ++j)
        {
            row_prod = (syndrome[j]) ? -1. : 1.;
            for (size_t i = 0; i < check_nodes[j].size(); ++i)
            {
                row_prod *= bit_to_check_msg[j][i];
            }

            for (size_t i = 0; i < check_nodes[j].size(); ++i)
            {
                prod = row_prod / bit_to_check_msg[j][i];
                curr_bit_pos = matrix.check_nodes[j][i];
                check_to_bit_msg[curr_bit_pos][check_pos_idx[curr_bit_pos]] = 2. * atanh(prod);
                check_pos_idx[curr_bit_pos]++;
            }
        }

        if (CFG.ENABLE_DECODING_ALG_MSG_LLR_THRESHOLD)
        {
            threshold_matrix(check_to_bit_msg, msg_threshold);
        }
        if (CFG.TRACE_DECODING_ALG)
        {
            fmt::print(fg(fmt::color::blue), "\nE:\n");
            print_matrix(check_to_bit_msg);
        }

        for (size_t i = 0; i < num_bit_nodes; ++i)
        {
            total_bit_llr[i] = std::accumulate(check_to_bit_msg[i].begin(), check_to_bit_msg[i].end(), bit_array_llr[i]);
            if (total_bit_llr[i] <= 0)
            {
                bit_array_out[i] = 1;
            }
            else
            {
                bit_array_out[i] = 0;
            }
        }

        if (CFG.TRACE_DECODING_ALG)
        {
            fmt::print(fg(fmt::color::blue), "\nL:\n");
            print_array(total_bit_llr);
            fmt::print(fg(fmt::color::blue), "\n\nz:\n");
            print_array(bit_array_out);
        }

        calculate_syndrome(bit_array_out, matrix, decision_syndrome);

        if (CFG.TRACE_DECODING_ALG)
        {
            fmt::print(fg(fmt::color::blue), "\n\ns:\n");
            print_array(decision_syndrome);
        }

        if (arrays_equal(syndrome, decision_syndrome))
        {
            if (CFG.TRACE_DECODING_ALG_LLR)
            {
                fmt::print(fg(fmt::color::blue), "\n\nMAX_LLR = {}\n", max_llr);
            }
            return {curr_iteration + 1, true};
        }

        std::fill(bit_pos_idx.begin(), bit_pos_idx.end(), 0);
        for (size_t i = 0; i < num_bit_nodes; ++i)
        {
            col_sum = total_bit_llr[i];
            for (size_t j = 0; j < bit_nodes[i].size(); ++j)
            {
                sum = col_sum - check_to_bit_msg[i][j];
                curr_check_pos = matrix.bit_nodes[i][j];
                bit_to_check_msg[curr_check_pos][bit_pos_idx[curr_check_pos]] = sum;
                bit_pos_idx[curr_check_pos]++;
            }
        }

        if (CFG.ENABLE_DECODING_ALG_MSG_LLR_THRESHOLD)
        {
            threshold_matrix(bit_to_check_msg, msg_threshold);
        }
        if (CFG.TRACE_DECODING_ALG)
        {
            fmt::print(fg(fmt::color::blue), "\n\nM:\n");
            print_matrix(bit_to_check_msg);
        }
        if (CFG.TRACE_DECODING_ALG_LLR)
        {
            max_llr_c2b = get_max_llr(check_to_bit_msg);
            max_llr_b2c = get_max_llr(bit_to_check_msg);
            max_llr = std::max({max_llr, max_llr_c2b, max_llr_b2c});
        }

        curr_iteration++;
    }

    if (CFG.TRACE_DECODING_ALG_LLR)
    {
        fmt::print(fg(fmt::color::blue), "\n\nMAX_LLR = {}\n", max_llr);
    }

    return {max_num_iterations, false};
}

double get_product_sign(const std::vector<double>& array)
{
    int negative_count = 0;
    for (size_t i = 0; i < array.size(); ++i)
    {
        if (array[i] < 0.)
        {
            ++negative_count;
        }
    }
    return (negative_count % 2 == 0) ? 1. : -1.;
}

double get_min_abs_except(const std::vector<double>& array, const size_t &k)
{
    double min_abs = std::numeric_limits<double>::max();
    double abs = 0;
    for (size_t i = 0; i < array.size(); ++i) 
    {
        if (i == k) continue; 

        abs = std::abs(array[i]);
        if (abs < min_abs) 
        {
            min_abs = abs;
        }
    }

    return min_abs;
}

decoding_result min_sum_normalized_decoding(const std::vector<double> &bit_array_llr,
                                            const H_matrix &matrix, 
                                            const std::vector<int> &syndrome,              
                                            const size_t &max_num_iterations,
                                            const double &alpha,   
                                            const double &msg_threshold,      
                                            std::vector<int> &bit_array_out)
{
    double max_llr_c2b = 0.;
    double max_llr_b2c = 0.;
    double max_llr = 0.;

    auto& bit_nodes = matrix.bit_nodes;
    auto& check_nodes = matrix.check_nodes;
    size_t num_bit_nodes = bit_nodes.size();
    size_t num_check_nodes = check_nodes.size();

    std::vector<std::vector<double>> bit_to_check_msg(num_check_nodes);
    for (size_t i = 0; i < num_check_nodes; ++i)
    {
        bit_to_check_msg[i].resize(check_nodes[i].size());
        for (size_t j = 0; j < check_nodes[i].size(); ++j)
        {
            bit_to_check_msg[i][j] = bit_array_llr[check_nodes[i][j]]; // Initialization
        }
    }

    std::vector<std::vector<double>> check_to_bit_msg(num_bit_nodes);
    for (size_t i = 0; i < num_bit_nodes; ++i)
    {
        check_to_bit_msg[i].resize(bit_nodes[i].size());
    }

    double prod;
    double sign_prod;
    int curr_bit_pos;

    std::vector<int> check_pos_idx(num_bit_nodes);
    std::vector<double> total_bit_llr(num_bit_nodes);
    std::vector<int> decision_syndrome(num_check_nodes);
    std::vector<int> bit_pos_idx(num_check_nodes);

    double sum;
    double col_sum;
    int curr_check_pos;

    size_t curr_iteration = 0;
    while (curr_iteration != max_num_iterations)
    {
        if (CFG.TRACE_DECODING_ALG)
        {
            fmt::print(fg(fmt::color::blue), "\n\nIteration: {}\n", curr_iteration + 1);
        }

        // Compute extrinsic messages from check nodes to bit nodes (Step 1: Check messages)
        std::fill(check_pos_idx.begin(), check_pos_idx.end(), 0);
        for (size_t j = 0; j < num_check_nodes; ++j)
        {
            sign_prod = (syndrome[j]) ? -1. : 1.;
            sign_prod *= get_product_sign(bit_to_check_msg[j]);     // Product of all signs in a row

            for (size_t i = 0; i < check_nodes[j].size(); ++i)
            {
                
                prod = sign_prod * ((bit_to_check_msg[j][i] > 0) ? 1. : -1.);       // Exclusion of the i-th sign
                curr_bit_pos = matrix.check_nodes[j][i];
                check_to_bit_msg[curr_bit_pos][check_pos_idx[curr_bit_pos]] = alpha * prod * get_min_abs_except(bit_to_check_msg[j], i);
                check_pos_idx[curr_bit_pos]++;
            }
        }

        if (CFG.ENABLE_DECODING_ALG_MSG_LLR_THRESHOLD)
        {
            threshold_matrix(check_to_bit_msg, msg_threshold);
        }
        if (CFG.TRACE_DECODING_ALG)
        {
            fmt::print(fg(fmt::color::blue), "\nE:\n");
            print_matrix(check_to_bit_msg);
        }

        for (size_t i = 0; i < num_bit_nodes; ++i)
        {
            total_bit_llr[i] = std::accumulate(check_to_bit_msg[i].begin(), check_to_bit_msg[i].end(), bit_array_llr[i]);
            if (total_bit_llr[i] <= 0)
            {
                bit_array_out[i] = 1;
            }
            else
            {
                bit_array_out[i] = 0;
            }
        }

        if (CFG.TRACE_DECODING_ALG)
        {
            fmt::print(fg(fmt::color::blue), "\nL:\n");
            print_array(total_bit_llr);
            fmt::print(fg(fmt::color::blue), "\n\nz:\n");
            print_array(bit_array_out);
        }

        calculate_syndrome(bit_array_out, matrix, decision_syndrome);

        if (CFG.TRACE_DECODING_ALG)
        {
            fmt::print(fg(fmt::color::blue), "\n\ns:\n");
            print_array(decision_syndrome);
        }

        if (arrays_equal(syndrome, decision_syndrome))
        {
            if (CFG.TRACE_DECODING_ALG_LLR)
            {
                fmt::print(fg(fmt::color::blue), "\n\nMAX_LLR = {}\n", max_llr);
            }
            return {curr_iteration + 1, true};
        }

        std::fill(bit_pos_idx.begin(), bit_pos_idx.end(), 0);
        for (size_t i = 0; i < num_bit_nodes; ++i)
        {
            col_sum = total_bit_llr[i];
            for (size_t j = 0; j < bit_nodes[i].size(); ++j)
            {
                sum = col_sum - check_to_bit_msg[i][j];
                curr_check_pos = matrix.bit_nodes[i][j];
                bit_to_check_msg[curr_check_pos][bit_pos_idx[curr_check_pos]] = sum;
                bit_pos_idx[curr_check_pos]++;
            }
        }

        if (CFG.ENABLE_DECODING_ALG_MSG_LLR_THRESHOLD)
        {
            threshold_matrix(bit_to_check_msg, msg_threshold);
        }
        if (CFG.TRACE_DECODING_ALG)
        {
            fmt::print(fg(fmt::color::blue), "\n\nM:\n");
            print_matrix(bit_to_check_msg);
        }
        if (CFG.TRACE_DECODING_ALG_LLR)
        {
            max_llr_c2b = get_max_llr(check_to_bit_msg);
            max_llr_b2c = get_max_llr(bit_to_check_msg);
            max_llr = std::max({max_llr, max_llr_c2b, max_llr_b2c});
        }

        curr_iteration++;
    }

    if (CFG.TRACE_DECODING_ALG_LLR)
    {
        fmt::print(fg(fmt::color::blue), "\n\nMAX_LLR = {}\n", max_llr);
    }

    return {max_num_iterations, false};
}

LDPC_result QKD_LDPC(const H_matrix &matrix,
                     const std::vector<int> &alice_bit_array, 
                     const std::vector<int> &bob_bit_array, 
                     const double &QBER,
                     const double &alpha)
{
    size_t num_bit_nodes = matrix.bit_nodes.size();
    size_t num_check_nodes = matrix.check_nodes.size();

    double log_p = log((1. - QBER) / QBER);         
    std::vector<double> apriori_llr(num_bit_nodes);

    for (size_t i = 0; i < num_bit_nodes; ++i)
    {
        apriori_llr[i] = (bob_bit_array[i] ? -log_p : log_p);
    }

    if (CFG.TRACE_QKD_LDPC)
    {
        fmt::print(fg(fmt::color::blue), "\nr:\n");
        print_array(apriori_llr);
    }

    std::vector<int> alice_syndrome(num_check_nodes);
    calculate_syndrome(alice_bit_array, matrix, alice_syndrome);

    if (CFG.TRACE_QKD_LDPC)
    {
        fmt::print(fg(fmt::color::blue), "\n\nAlice syndrome:\n");
        print_array(alice_syndrome);
    }

    std::vector<int> bob_solution(num_bit_nodes);
    LDPC_result ldpc_res;
    if (CFG.USE_MIN_SUM_NORMALIZED_ALG)
    {
        ldpc_res.decoding_res = min_sum_normalized_decoding(apriori_llr, matrix, alice_syndrome, CFG.DECODING_ALG_MAX_ITERATIONS, 
                                                 alpha, CFG.DECODING_ALG_MSG_LLR_THRESHOLD, bob_solution);
    }
    else
    {
        ldpc_res.decoding_res = sum_product_decoding(apriori_llr, matrix, alice_syndrome, CFG.DECODING_ALG_MAX_ITERATIONS,
                                                     CFG.DECODING_ALG_MSG_LLR_THRESHOLD, bob_solution);
    }

    if (CFG.TRACE_QKD_LDPC)
    {
        fmt::print(fg(fmt::color::blue), "\nBob corrected bit array:\n");
        print_array(bob_solution);
    }

    ldpc_res.keys_match = arrays_equal(alice_bit_array, bob_solution);

    if (CFG.ENABLE_PRIVACY_MAINTENANCE)
    {
        std::vector<int> alice_bit_array_pm;
        std::vector<int> bob_bit_array_pm;
        privacy_maintenance(matrix, alice_bit_array, bob_solution, alice_bit_array_pm, bob_bit_array_pm);
        if (CFG.TRACE_QKD_LDPC)
        {
            fmt::print(fg(fmt::color::blue), "\nAlice bit array after privacy maintenance:\n");
            print_array(alice_bit_array_pm);
            fmt::print(fg(fmt::color::blue), "\nBob bit array after privacy maintenance:\n");
            print_array(bob_bit_array_pm);
        }
    }    

    if (CFG.TRACE_QKD_LDPC)
    {
        fmt::print(fg(fmt::color::blue), "\n\nIterations performed: {}\n", ldpc_res.decoding_res.iterations_num);
        fmt::print(fg(fmt::color::blue), "Syndromes are match: {}\n", ((ldpc_res.decoding_res.syndromes_match) ? "YES" : "NO"));
        fmt::print(fg(fmt::color::blue), "Keys are match: {}\n", ((ldpc_res.keys_match) ? "YES" : "NO"));
    }

    return ldpc_res;
}
