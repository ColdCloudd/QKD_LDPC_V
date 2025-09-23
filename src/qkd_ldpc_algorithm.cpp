#include "qkd_ldpc_algorithm.hpp"

decoding_result sum_product_decoding(
    const std::vector<double> &bit_array_llr,
    const H_matrix &matrix,
    const std::vector<int> &syndrome,          
    const size_t &max_num_iterations, 
    const double &msg_threshold, 
    std::vector<int> &bit_array_out
)
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
        // Compute extrinsic messages from check nodes to bit nodes (Step 1: Check messages)
        std::fill(check_pos_idx.begin(), check_pos_idx.end(), 0);
        for (size_t j = 0; j < num_check_nodes; ++j)
        {
            row_prod = (syndrome[j]) ? -1. : 1.;
            for (size_t i = 0; i < check_nodes[j].size(); ++i)
            {
                bit_to_check_msg[j][i] = tanh(bit_to_check_msg[j][i] / 2.);
                row_prod *= bit_to_check_msg[j][i];
            }

            for (size_t i = 0; i < check_nodes[j].size(); ++i)
            {
                prod = row_prod / bit_to_check_msg[j][i];
                curr_bit_pos = matrix.check_nodes[j][i];
                check_to_bit_msg[curr_bit_pos][check_pos_idx[curr_bit_pos]] = 2. * atanh(prod);
                ++check_pos_idx[curr_bit_pos];
            }
        }

        if (CFG.ENABLE_DECODING_ALG_MSG_LLR_THRESHOLD)
            threshold_matrix(check_to_bit_msg, msg_threshold);

        for (size_t i = 0; i < num_bit_nodes; ++i)
        {
            total_bit_llr[i] = std::accumulate(check_to_bit_msg[i].begin(), check_to_bit_msg[i].end(), bit_array_llr[i]);

            if (total_bit_llr[i] <= 0)
                bit_array_out[i] = 1;
            else
                bit_array_out[i] = 0;
        }

        calculate_syndrome(bit_array_out, matrix, decision_syndrome);

        if (CFG.TRACE_DECODING_ALG)
        {
            fmt::print(fg(fmt::color::blue), "\n\nIteration: {}\n", curr_iteration + 1);
            fmt::print(fg(fmt::color::blue), "\nE:\n");
            print_matrix(check_to_bit_msg);
            fmt::print(fg(fmt::color::blue), "\nL:\n");
            print_array(total_bit_llr);
            fmt::print(fg(fmt::color::blue), "\nz:\n");
            print_array(bit_array_out);
            fmt::print(fg(fmt::color::blue), "\ns:\n");
            print_array(decision_syndrome);
        }

        if (arrays_equal(syndrome, decision_syndrome))
        {
            if (CFG.TRACE_DECODING_ALG_LLR)
                fmt::print(fg(fmt::color::blue), "\nMAX_LLR = {}\n", max_llr);

            return {.iterations_num = (curr_iteration + 1), .syndromes_match = true};
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
                ++bit_pos_idx[curr_check_pos];
            }
        }

        if (CFG.ENABLE_DECODING_ALG_MSG_LLR_THRESHOLD)
            threshold_matrix(bit_to_check_msg, msg_threshold);

        if (CFG.TRACE_DECODING_ALG)
        {
            fmt::print(fg(fmt::color::blue), "\nM:\n");
            print_matrix(bit_to_check_msg);
        }
        if (CFG.TRACE_DECODING_ALG_LLR)
        {
            max_llr_c2b = get_max_llr(check_to_bit_msg);
            max_llr_b2c = get_max_llr(bit_to_check_msg);
            max_llr = std::max({max_llr, max_llr_c2b, max_llr_b2c});
        }

        ++curr_iteration;
    }

    if (CFG.TRACE_DECODING_ALG_LLR)
        fmt::print(fg(fmt::color::blue), "\nMAX_LLR = {}\n", max_llr);

    return {.iterations_num = max_num_iterations, .syndromes_match = false};
}

double tanh_lin_approx(double x)
{
    double abs_x = std::abs(x);
    double result;
    if (abs_x < 0.5) result = 0.9242 * abs_x;
    else if (abs_x < 0.9) result = 0.6355 * abs_x + 0.1444;
    else if (abs_x < 1.2) result = 0.3912 * abs_x + 0.3642;
    else if (abs_x < 1.75) result = 0.1958 * abs_x + 0.5986;
    else if (abs_x < 2.5) result = 0.0603 * abs_x + 0.8358;
    else if (abs_x < 3.5) result = 0.0115 * abs_x + 0.9577;
    else if (abs_x < 8) result = 0.0004 * abs_x + 0.9967;
    else result = 1;

    return (x < 0.) ? -result : result;
}

double atanh_lin_approx(double x)
{
    double abs_x = std::abs(x);
    double result;
    if (abs_x < 0.7) result = 1.196 * abs_x - 0.0323;
    else if (abs_x < 0.9) result = 2.9187 * abs_x - 1.214;
    else if (abs_x < 0.999) result = 10.8717 * abs_x - 8.3717;
    else result = 2510.9 * abs_x - 2505.9;

    return (x < 0.) ? -result : result;
}

decoding_result sum_product_linear_approx_decoding(
    const std::vector<double> &bit_array_llr,
    const H_matrix &matrix,
    const std::vector<int> &syndrome,          
    const size_t &max_num_iterations, 
    const double &msg_threshold, 
    std::vector<int> &bit_array_out
)
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
        // Compute extrinsic messages from check nodes to bit nodes (Step 1: Check messages)
        std::fill(check_pos_idx.begin(), check_pos_idx.end(), 0);
        for (size_t j = 0; j < num_check_nodes; ++j)
        {
            row_prod = (syndrome[j]) ? -1. : 1.;
            for (size_t i = 0; i < check_nodes[j].size(); ++i)
            {
                bit_to_check_msg[j][i] = tanh_lin_approx(bit_to_check_msg[j][i] / 2.);
                row_prod *= bit_to_check_msg[j][i];
            }

            for (size_t i = 0; i < check_nodes[j].size(); ++i)
            {
                prod = row_prod / bit_to_check_msg[j][i];
                curr_bit_pos = matrix.check_nodes[j][i];
                check_to_bit_msg[curr_bit_pos][check_pos_idx[curr_bit_pos]] = 2. * atanh_lin_approx(prod);
                ++check_pos_idx[curr_bit_pos];
            }
        }

        if (CFG.ENABLE_DECODING_ALG_MSG_LLR_THRESHOLD)
            threshold_matrix(check_to_bit_msg, msg_threshold);

        for (size_t i = 0; i < num_bit_nodes; ++i)
        {
            total_bit_llr[i] = std::accumulate(check_to_bit_msg[i].begin(), check_to_bit_msg[i].end(), bit_array_llr[i]);

            if (total_bit_llr[i] <= 0)
                bit_array_out[i] = 1;
            else
                bit_array_out[i] = 0;
        }

        calculate_syndrome(bit_array_out, matrix, decision_syndrome);

        if (CFG.TRACE_DECODING_ALG)
        {
            fmt::print(fg(fmt::color::blue), "\n\nIteration: {}\n", curr_iteration + 1);
            fmt::print(fg(fmt::color::blue), "\nE:\n");
            print_matrix(check_to_bit_msg);
            fmt::print(fg(fmt::color::blue), "\nL:\n");
            print_array(total_bit_llr);
            fmt::print(fg(fmt::color::blue), "\nz:\n");
            print_array(bit_array_out);
            fmt::print(fg(fmt::color::blue), "\ns:\n");
            print_array(decision_syndrome);
        }

        if (arrays_equal(syndrome, decision_syndrome))
        {
            if (CFG.TRACE_DECODING_ALG_LLR)
                fmt::print(fg(fmt::color::blue), "\nMAX_LLR = {}\n", max_llr);

            return {.iterations_num = (curr_iteration + 1), .syndromes_match = true};
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
                ++bit_pos_idx[curr_check_pos];
            }
        }

        if (CFG.ENABLE_DECODING_ALG_MSG_LLR_THRESHOLD)
            threshold_matrix(bit_to_check_msg, msg_threshold);

        if (CFG.TRACE_DECODING_ALG)
        {
            fmt::print(fg(fmt::color::blue), "\nM:\n");
            print_matrix(bit_to_check_msg);
        }
        if (CFG.TRACE_DECODING_ALG_LLR)
        {
            max_llr_c2b = get_max_llr(check_to_bit_msg);
            max_llr_b2c = get_max_llr(bit_to_check_msg);
            max_llr = std::max({max_llr, max_llr_c2b, max_llr_b2c});
        }

        ++curr_iteration;
    }

    if (CFG.TRACE_DECODING_ALG_LLR)
        fmt::print(fg(fmt::color::blue), "\nMAX_LLR = {}\n", max_llr);

    return {.iterations_num = max_num_iterations, .syndromes_match = false};
}

decoding_result min_sum_normalized_decoding(
    const std::vector<double> &bit_array_llr,
    const H_matrix &matrix, 
    const std::vector<int> &syndrome,              
    const size_t &max_num_iterations,
    const double &alpha,   
    const double &msg_threshold,      
    std::vector<int> &bit_array_out
)
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
    double bit_to_check_msg_min_1;
    double bit_to_check_msg_min_2;
    double current_abs;
    int negative_count;
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
        // Compute extrinsic messages from check nodes to bit nodes (Step 1: Check messages)
        std::fill(check_pos_idx.begin(), check_pos_idx.end(), 0);
        for (size_t j = 0; j < num_check_nodes; ++j)
        {
            sign_prod = (syndrome[j]) ? -1. : 1.;
            negative_count = 0;
            bit_to_check_msg_min_1 = std::numeric_limits<double>::max();      
            bit_to_check_msg_min_2 = std::numeric_limits<double>::max();  

            for (size_t k = 0; k < bit_to_check_msg[j].size(); ++k)         // Searching for two minimum values
            {
                if (bit_to_check_msg[j][k] < 0) ++negative_count;           // Count negative numbers to determine the sign of the product
                
                current_abs = std::abs(bit_to_check_msg[j][k]);
                if (current_abs < bit_to_check_msg_min_1) 
                {
                    // The current number becomes the new minimum number (modulo)
                    bit_to_check_msg_min_2 = bit_to_check_msg_min_1;
                    bit_to_check_msg_min_1 = current_abs;
                } 
                else if (current_abs < bit_to_check_msg_min_2) 
                {
                    // Update second minimum modulo number
                    bit_to_check_msg_min_2 = current_abs;
                }
            }
            sign_prod *= (negative_count % 2 == 0) ? 1. : -1.;        // Product of all signs in a row
            
            for (size_t i = 0; i < check_nodes[j].size(); ++i)
            {
                prod = sign_prod * ((bit_to_check_msg[j][i] > 0) ? 1. : -1.);       // Exclusion of the i-th sign
                curr_bit_pos = matrix.check_nodes[j][i];
                // If the i-th modulo value is equal to the found minimum value, then we use the second minimum value
                check_to_bit_msg[curr_bit_pos][check_pos_idx[curr_bit_pos]] = alpha * prod * 
                ((std::abs(bit_to_check_msg[j][i]) == bit_to_check_msg_min_1) ? bit_to_check_msg_min_2 : bit_to_check_msg_min_1); 
                ++check_pos_idx[curr_bit_pos];
            }
        }

        if (CFG.ENABLE_DECODING_ALG_MSG_LLR_THRESHOLD)
            threshold_matrix(check_to_bit_msg, msg_threshold);

        for (size_t i = 0; i < num_bit_nodes; ++i)
        {
            total_bit_llr[i] = std::accumulate(check_to_bit_msg[i].begin(), check_to_bit_msg[i].end(), bit_array_llr[i]);
            
            if (total_bit_llr[i] <= 0)
                bit_array_out[i] = 1;
            else
                bit_array_out[i] = 0;
        }

        calculate_syndrome(bit_array_out, matrix, decision_syndrome);

        if (CFG.TRACE_DECODING_ALG)
        {
            fmt::print(fg(fmt::color::blue), "\n\nIteration: {}\n", curr_iteration + 1);
            fmt::print(fg(fmt::color::blue), "\nE:\n");
            print_matrix(check_to_bit_msg);
            fmt::print(fg(fmt::color::blue), "\nL:\n");
            print_array(total_bit_llr);
            fmt::print(fg(fmt::color::blue), "\nz:\n");
            print_array(bit_array_out);
            fmt::print(fg(fmt::color::blue), "\ns:\n");
            print_array(decision_syndrome);
        }

        if (arrays_equal(syndrome, decision_syndrome))
        {
            if (CFG.TRACE_DECODING_ALG_LLR)
                fmt::print(fg(fmt::color::blue), "\nMAX_LLR = {}\n", max_llr);

            return {.iterations_num = (curr_iteration + 1), .syndromes_match = true};
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
                ++bit_pos_idx[curr_check_pos];
            }
        }

        if (CFG.ENABLE_DECODING_ALG_MSG_LLR_THRESHOLD)
            threshold_matrix(bit_to_check_msg, msg_threshold);

        if (CFG.TRACE_DECODING_ALG)
        {
            fmt::print(fg(fmt::color::blue), "\nM:\n");
            print_matrix(bit_to_check_msg);
        }
        if (CFG.TRACE_DECODING_ALG_LLR)
        {
            max_llr_c2b = get_max_llr(check_to_bit_msg);
            max_llr_b2c = get_max_llr(bit_to_check_msg);
            max_llr = std::max({max_llr, max_llr_c2b, max_llr_b2c});
        }

        ++curr_iteration;
    }

    if (CFG.TRACE_DECODING_ALG_LLR)
        fmt::print(fg(fmt::color::blue), "\nMAX_LLR = {}\n", max_llr);

    return {.iterations_num = max_num_iterations, .syndromes_match = false};
}

decoding_result min_sum_offset_decoding(
    const std::vector<double> &bit_array_llr,
    const H_matrix &matrix, 
    const std::vector<int> &syndrome,              
    const size_t &max_num_iterations,
    const double &beta,   
    const double &msg_threshold,      
    std::vector<int> &bit_array_out
)
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

    double diff;
    double prod;
    double sign_prod;
    double bit_to_check_msg_min_1;
    double bit_to_check_msg_min_2;
    double current_abs;
    int negative_count;
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
        // Compute extrinsic messages from check nodes to bit nodes (Step 1: Check messages)
        std::fill(check_pos_idx.begin(), check_pos_idx.end(), 0);
        for (size_t j = 0; j < num_check_nodes; ++j)
        {
            sign_prod = (syndrome[j]) ? -1. : 1.;
            negative_count = 0;
            bit_to_check_msg_min_1 = std::numeric_limits<double>::max();      
            bit_to_check_msg_min_2 = std::numeric_limits<double>::max();  

            for (size_t k = 0; k < bit_to_check_msg[j].size(); ++k)         // Searching for two minimum values
            {
                if (bit_to_check_msg[j][k] < 0) ++negative_count;           // Count negative numbers to determine the sign of the product
                
                current_abs = std::abs(bit_to_check_msg[j][k]);
                if (current_abs < bit_to_check_msg_min_1) 
                {
                    // The current number becomes the new minimum number (modulo)
                    bit_to_check_msg_min_2 = bit_to_check_msg_min_1;
                    bit_to_check_msg_min_1 = current_abs;
                } 
                else if (current_abs < bit_to_check_msg_min_2) 
                {
                    // Update second minimum modulo number
                    bit_to_check_msg_min_2 = current_abs;
                }
            }
            sign_prod *= (negative_count % 2 == 0) ? 1. : -1.;        // Product of all signs in a row

            for (size_t i = 0; i < check_nodes[j].size(); ++i)
            {
                prod = sign_prod * ((bit_to_check_msg[j][i] > 0) ? 1. : -1.);       // Exclusion of the i-th sign
                curr_bit_pos = matrix.check_nodes[j][i];
                // If the i-th modulo value is equal to the found minimum value, then we use the second minimum value
                diff = ((std::abs(bit_to_check_msg[j][i]) == bit_to_check_msg_min_1) ? bit_to_check_msg_min_2 : bit_to_check_msg_min_1) - beta;
                check_to_bit_msg[curr_bit_pos][check_pos_idx[curr_bit_pos]] = prod * ((diff < 0.) ? 0. : diff);        // If the absolute value of the selected constant 'beta' exceeds the minimum absolute value 
                ++check_pos_idx[curr_bit_pos];                                                                         // of the message 'bit_to_check_msg', then 'check_to_bit_msg' converts to 0
            }
        }

        if (CFG.ENABLE_DECODING_ALG_MSG_LLR_THRESHOLD)
            threshold_matrix(check_to_bit_msg, msg_threshold);

        for (size_t i = 0; i < num_bit_nodes; ++i)
        {
            total_bit_llr[i] = std::accumulate(check_to_bit_msg[i].begin(), check_to_bit_msg[i].end(), bit_array_llr[i]);

            if (total_bit_llr[i] <= 0)
                bit_array_out[i] = 1;
            else
                bit_array_out[i] = 0;
        }

        calculate_syndrome(bit_array_out, matrix, decision_syndrome);

        if (CFG.TRACE_DECODING_ALG)
        {
            fmt::print(fg(fmt::color::blue), "\n\nIteration: {}\n", curr_iteration + 1);
            fmt::print(fg(fmt::color::blue), "\nE:\n");
            print_matrix(check_to_bit_msg);
            fmt::print(fg(fmt::color::blue), "\nL:\n");
            print_array(total_bit_llr);
            fmt::print(fg(fmt::color::blue), "\nz:\n");
            print_array(bit_array_out);
            fmt::print(fg(fmt::color::blue), "\ns:\n");
            print_array(decision_syndrome);
        }

        if (arrays_equal(syndrome, decision_syndrome))
        {
            if (CFG.TRACE_DECODING_ALG_LLR)
                fmt::print(fg(fmt::color::blue), "\nMAX_LLR = {}\n", max_llr);

            return {.iterations_num = (curr_iteration + 1), .syndromes_match = true};
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
                ++bit_pos_idx[curr_check_pos];
            }
        }

        if (CFG.ENABLE_DECODING_ALG_MSG_LLR_THRESHOLD)
            threshold_matrix(bit_to_check_msg, msg_threshold);

        if (CFG.TRACE_DECODING_ALG)
        {
            fmt::print(fg(fmt::color::blue), "\nM:\n");
            print_matrix(bit_to_check_msg);
        }
        if (CFG.TRACE_DECODING_ALG_LLR)
        {
            max_llr_c2b = get_max_llr(check_to_bit_msg);
            max_llr_b2c = get_max_llr(bit_to_check_msg);
            max_llr = std::max({max_llr, max_llr_c2b, max_llr_b2c});
        }

        ++curr_iteration;
    }

    if (CFG.TRACE_DECODING_ALG_LLR)
        fmt::print(fg(fmt::color::blue), "\nMAX_LLR = {}\n", max_llr);

    return {.iterations_num = max_num_iterations, .syndromes_match = false};
}                                

decoding_result adaptive_min_sum_normalized_decoding(
    const std::vector<double> &bit_array_llr,
    const H_matrix &matrix, 
    const std::vector<int> &syndrome,              
    const size_t &max_num_iterations,
    const double &alpha, 
    const double &nu, 
    const double &msg_threshold,      
    std::vector<int> &bit_array_out
)
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
            bit_to_check_msg[i][j] = bit_array_llr[check_nodes[i][j]];      // Initialization
        }
    }

    std::vector<std::vector<double>> check_to_bit_msg(num_bit_nodes);
    for (size_t i = 0; i < num_bit_nodes; ++i)
    {
        if (bit_array_llr[i] <= 0)
            bit_array_out[i] = 1;
        else
            bit_array_out[i] = 0;

        check_to_bit_msg[i].resize(bit_nodes[i].size());
    }

    double prod;
    double sign_prod;
    double bit_to_check_msg_min_1;
    double bit_to_check_msg_min_2;
    double current_abs;
    double scaling_factor;
    bool syndromes_equal;
    int negative_count;
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
        // Compute extrinsic messages from check nodes to bit nodes (Step 1: Check messages)
        syndromes_equal = true;
        std::fill(check_pos_idx.begin(), check_pos_idx.end(), 0);
        std::fill(decision_syndrome.begin(), decision_syndrome.end(), 0);
        for (size_t j = 0; j < num_check_nodes; ++j)
        {
            sign_prod = (syndrome[j]) ? -1. : 1.;
            negative_count = 0;
            bit_to_check_msg_min_1 = std::numeric_limits<double>::max();      
            bit_to_check_msg_min_2 = std::numeric_limits<double>::max();  

            for (size_t k = 0; k < bit_to_check_msg[j].size(); ++k)         // Searching for two minimum values
            {
                if (bit_to_check_msg[j][k] < 0) ++negative_count;           // Count negative numbers to determine the sign of the product
                
                current_abs = std::abs(bit_to_check_msg[j][k]);
                if (current_abs < bit_to_check_msg_min_1) 
                {
                    // The current number becomes the new minimum number (modulo)
                    bit_to_check_msg_min_2 = bit_to_check_msg_min_1;
                    bit_to_check_msg_min_1 = current_abs;
                } 
                else if (current_abs < bit_to_check_msg_min_2) 
                {
                    // Update second minimum modulo number
                    bit_to_check_msg_min_2 = current_abs;
                }
            }
            sign_prod *= (negative_count % 2 == 0) ? 1. : -1.;        // Product of all signs in a row

            for (size_t k = 0; k < check_nodes[j].size(); ++k)      
            {
                decision_syndrome[j] ^= bit_array_out[check_nodes[j][k]];       // Calculation of j-check value
            }
            if (decision_syndrome[j] != syndrome[j])        // Messages from this node are considered unreliable and are suppressed
            {
                scaling_factor = nu;
                syndromes_equal = false;
            }
            else
            {
                scaling_factor = alpha;
            }
            
            for (size_t i = 0; i < check_nodes[j].size(); ++i)
            {
                prod = sign_prod * ((bit_to_check_msg[j][i] > 0) ? 1. : -1.);       // Exclusion of the i-th sign
                curr_bit_pos = matrix.check_nodes[j][i];
                // If the i-th modulo value is equal to the found minimum value, then we use the second minimum value
                check_to_bit_msg[curr_bit_pos][check_pos_idx[curr_bit_pos]] = scaling_factor * prod * 
                ((std::abs(bit_to_check_msg[j][i]) == bit_to_check_msg_min_1) ? bit_to_check_msg_min_2 : bit_to_check_msg_min_1); 
                ++check_pos_idx[curr_bit_pos];
            }
        }
        
        if (syndromes_equal)
        {
            if (CFG.TRACE_DECODING_ALG_LLR)
                fmt::print(fg(fmt::color::blue), "\nMAX_LLR = {}\n", max_llr);

            return {.iterations_num = (curr_iteration + 1), .syndromes_match = true};
        }

        if (CFG.ENABLE_DECODING_ALG_MSG_LLR_THRESHOLD)
            threshold_matrix(check_to_bit_msg, msg_threshold);

        for (size_t i = 0; i < num_bit_nodes; ++i)
        {
            total_bit_llr[i] = std::accumulate(check_to_bit_msg[i].begin(), check_to_bit_msg[i].end(), bit_array_llr[i]);

            if (total_bit_llr[i] <= 0)
                bit_array_out[i] = 1;
            else
                bit_array_out[i] = 0;
        }

        if (CFG.TRACE_DECODING_ALG)
        {
            fmt::print(fg(fmt::color::blue), "\n\nIteration: {}\n", curr_iteration + 1);
            fmt::print(fg(fmt::color::blue), "\nE:\n");
            print_matrix(check_to_bit_msg);
            fmt::print(fg(fmt::color::blue), "\nL:\n");
            print_array(total_bit_llr);
            fmt::print(fg(fmt::color::blue), "\nz:\n");
            print_array(bit_array_out);
            fmt::print(fg(fmt::color::blue), "\ns:\n");
            print_array(decision_syndrome);
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
                ++bit_pos_idx[curr_check_pos];
            }
        }

        if (CFG.ENABLE_DECODING_ALG_MSG_LLR_THRESHOLD)
            threshold_matrix(bit_to_check_msg, msg_threshold);

        if (CFG.TRACE_DECODING_ALG)
        {
            fmt::print(fg(fmt::color::blue), "\nM:\n");
            print_matrix(bit_to_check_msg);
        }
        if (CFG.TRACE_DECODING_ALG_LLR)
        {
            max_llr_c2b = get_max_llr(check_to_bit_msg);
            max_llr_b2c = get_max_llr(bit_to_check_msg);
            max_llr = std::max({max_llr, max_llr_c2b, max_llr_b2c});
        }

        ++curr_iteration;
    }

    if (CFG.TRACE_DECODING_ALG_LLR)
        fmt::print(fg(fmt::color::blue), "\nMAX_LLR = {}\n", max_llr);

    return {.iterations_num = max_num_iterations, .syndromes_match = false};
}

decoding_result adaptive_min_sum_offset_decoding(
    const std::vector<double> &bit_array_llr,
    const H_matrix &matrix, 
    const std::vector<int> &syndrome,              
    const size_t &max_num_iterations,
    const double &beta, 
    const double &sigma, 
    const double &msg_threshold,      
    std::vector<int> &bit_array_out
)
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
        if (bit_array_llr[i] <= 0)
            bit_array_out[i] = 1;
        else
            bit_array_out[i] = 0;

        check_to_bit_msg[i].resize(bit_nodes[i].size());
    }

    double diff;
    double prod;
    double sign_prod;
    double bit_to_check_msg_min_1;
    double bit_to_check_msg_min_2;
    double current_abs;
    double scaling_factor;
    bool syndromes_equal;
    int negative_count;
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
        // Compute extrinsic messages from check nodes to bit nodes (Step 1: Check messages)
        syndromes_equal = true;
        std::fill(check_pos_idx.begin(), check_pos_idx.end(), 0);
        std::fill(decision_syndrome.begin(), decision_syndrome.end(), 0);
        for (size_t j = 0; j < num_check_nodes; ++j)
        {
            sign_prod = (syndrome[j]) ? -1. : 1.;
            negative_count = 0;
            bit_to_check_msg_min_1 = std::numeric_limits<double>::max();      
            bit_to_check_msg_min_2 = std::numeric_limits<double>::max();  

            for (size_t k = 0; k < bit_to_check_msg[j].size(); ++k)         // Searching for two minimum values
            {
                if (bit_to_check_msg[j][k] < 0) ++negative_count;           // Count negative numbers to determine the sign of the product
                
                current_abs = std::abs(bit_to_check_msg[j][k]);
                if (current_abs < bit_to_check_msg_min_1) 
                {
                    // The current number becomes the new minimum number (modulo)
                    bit_to_check_msg_min_2 = bit_to_check_msg_min_1;
                    bit_to_check_msg_min_1 = current_abs;
                } 
                else if (current_abs < bit_to_check_msg_min_2) 
                {
                    // Update second minimum modulo number
                    bit_to_check_msg_min_2 = current_abs;
                }
            }
            sign_prod *= (negative_count % 2 == 0) ? 1. : -1.;        // Product of all signs in a row

            for (size_t k = 0; k < check_nodes[j].size(); ++k)   
            {
                decision_syndrome[j] ^= bit_array_out[check_nodes[j][k]];       // Calculation of j-check value
            }
            if (decision_syndrome[j] != syndrome[j])        // Messages from this node are considered unreliable and are suppressed
            {
                scaling_factor = sigma;
                syndromes_equal = false;
            }
            else
            {
                scaling_factor = beta;
            }
            
            for (size_t i = 0; i < check_nodes[j].size(); ++i)
            {
                prod = sign_prod * ((bit_to_check_msg[j][i] > 0) ? 1. : -1.);       // Exclusion of the i-th sign
                curr_bit_pos = matrix.check_nodes[j][i];
                // If the i-th modulo value is equal to the found minimum value, then we use the second minimum value
                diff = ((std::abs(bit_to_check_msg[j][i]) == bit_to_check_msg_min_1) ? bit_to_check_msg_min_2 : bit_to_check_msg_min_1) - scaling_factor;
                check_to_bit_msg[curr_bit_pos][check_pos_idx[curr_bit_pos]] = prod * ((diff < 0.) ? 0. : diff);        // If the absolute value of the selected constant 'beta' exceeds the minimum absolute value 
                ++check_pos_idx[curr_bit_pos];
            }
        }
        
        if (syndromes_equal)
        {
            if (CFG.TRACE_DECODING_ALG_LLR)
                fmt::print(fg(fmt::color::blue), "\nMAX_LLR = {}\n", max_llr);

            return {.iterations_num = (curr_iteration + 1), .syndromes_match = true};
        }

        if (CFG.ENABLE_DECODING_ALG_MSG_LLR_THRESHOLD)
            threshold_matrix(check_to_bit_msg, msg_threshold);

        for (size_t i = 0; i < num_bit_nodes; ++i)
        {
            total_bit_llr[i] = std::accumulate(check_to_bit_msg[i].begin(), check_to_bit_msg[i].end(), bit_array_llr[i]);

            if (total_bit_llr[i] <= 0)
                bit_array_out[i] = 1;
            else
                bit_array_out[i] = 0;
        }

        if (CFG.TRACE_DECODING_ALG)
        {
            fmt::print(fg(fmt::color::blue), "\n\nIteration: {}\n", curr_iteration + 1);
            fmt::print(fg(fmt::color::blue), "\nE:\n");
            print_matrix(check_to_bit_msg);
            fmt::print(fg(fmt::color::blue), "\nL:\n");
            print_array(total_bit_llr);
            fmt::print(fg(fmt::color::blue), "\nz:\n");
            print_array(bit_array_out);
            fmt::print(fg(fmt::color::blue), "\ns:\n");
            print_array(decision_syndrome);
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
                ++bit_pos_idx[curr_check_pos];
            }
        }

        if (CFG.ENABLE_DECODING_ALG_MSG_LLR_THRESHOLD)
            threshold_matrix(bit_to_check_msg, msg_threshold);

        if (CFG.TRACE_DECODING_ALG)
        {
            fmt::print(fg(fmt::color::blue), "\nM:\n");
            print_matrix(bit_to_check_msg);
        }
        if (CFG.TRACE_DECODING_ALG_LLR)
        {
            max_llr_c2b = get_max_llr(check_to_bit_msg);
            max_llr_b2c = get_max_llr(bit_to_check_msg);
            max_llr = std::max({max_llr, max_llr_c2b, max_llr_b2c});
        }

        ++curr_iteration;
    }

    if (CFG.TRACE_DECODING_ALG_LLR)
        fmt::print(fg(fmt::color::blue), "\nMAX_LLR = {}\n", max_llr);

    return {.iterations_num = max_num_iterations, .syndromes_match = false};
}

LDPC_result QKD_LDPC(
    const H_matrix &matrix,
    const std::vector<int> &alice_bit_array, 
    const std::vector<int> &bob_bit_array, 
    const double &QBER,
    const decoding_scaling_factors &scaling_factors,
    const H_matrix_params &matrix_params
)
{
    size_t num_bit_nodes = matrix.bit_nodes.size();
    size_t num_check_nodes = matrix.check_nodes.size();

    double log_p = log((1. - QBER) / QBER);         
    std::vector<double> apriori_llr(num_bit_nodes);

    for (size_t i = 0; i < num_bit_nodes; ++i)
    {
        apriori_llr[i] = (bob_bit_array[i] ? -log_p : log_p);
    }

    std::vector<int> alice_syndrome(num_check_nodes);
    calculate_syndrome(alice_bit_array, matrix, alice_syndrome);

    std::vector<int> bob_solution(num_bit_nodes);
    LDPC_result ldpc_res;
    if (CFG.DECODING_ALGORITHM == DEC_SPA)
    {
        ldpc_res.decoding_res = sum_product_decoding(apriori_llr, matrix, alice_syndrome, CFG.DECODING_ALG_MAX_ITERATIONS,
            CFG.DECODING_ALG_MSG_LLR_THRESHOLD, bob_solution);
    }
    else if (CFG.DECODING_ALGORITHM == DEC_SPA_APPROX)
    {
        ldpc_res.decoding_res = sum_product_linear_approx_decoding(apriori_llr, matrix, alice_syndrome, CFG.DECODING_ALG_MAX_ITERATIONS,
            CFG.DECODING_ALG_MSG_LLR_THRESHOLD, bob_solution);
    }
    else if (CFG.DECODING_ALGORITHM == DEC_NMSA)
    {
        ldpc_res.decoding_res = min_sum_normalized_decoding(apriori_llr, matrix, alice_syndrome, CFG.DECODING_ALG_MAX_ITERATIONS, 
            scaling_factors.primary, CFG.DECODING_ALG_MSG_LLR_THRESHOLD, bob_solution);
    }
    else if (CFG.DECODING_ALGORITHM == DEC_OMSA)
    {
        ldpc_res.decoding_res = min_sum_offset_decoding(apriori_llr, matrix, alice_syndrome, CFG.DECODING_ALG_MAX_ITERATIONS, 
            scaling_factors.primary, CFG.DECODING_ALG_MSG_LLR_THRESHOLD, bob_solution);
    }
    else if (CFG.DECODING_ALGORITHM == DEC_ANMSA)
    {
        ldpc_res.decoding_res = adaptive_min_sum_normalized_decoding(apriori_llr, matrix, alice_syndrome, CFG.DECODING_ALG_MAX_ITERATIONS, 
            scaling_factors.primary, scaling_factors.secondary, CFG.DECODING_ALG_MSG_LLR_THRESHOLD, bob_solution);
    }
    else if (CFG.DECODING_ALGORITHM == DEC_AOMSA)
    {
        ldpc_res.decoding_res = adaptive_min_sum_offset_decoding(apriori_llr, matrix, alice_syndrome, CFG.DECODING_ALG_MAX_ITERATIONS, 
            scaling_factors.primary, scaling_factors.secondary, CFG.DECODING_ALG_MSG_LLR_THRESHOLD, bob_solution);
    }

    ldpc_res.keys_match = arrays_equal(alice_bit_array, bob_solution);

    std::vector<int> alice_bit_array_pm;
    std::vector<int> bob_bit_array_pm;
    if (CFG.ENABLE_PRIVACY_MAINTENANCE)
        remove_bits(matrix_params.bits_to_remove, alice_bit_array, bob_solution, alice_bit_array_pm, bob_bit_array_pm);

    if (CFG.TRACE_QKD_LDPC)
    {
        fmt::print(fg(fmt::color::blue), "\nAlice bit array:\n");
        print_array(alice_bit_array);
        fmt::print(fg(fmt::color::blue), "\nBob bit array with errors:\n");
        print_array(bob_bit_array);
        fmt::print(fg(fmt::color::blue), "\nr:\n");
        print_array(apriori_llr);
        fmt::print(fg(fmt::color::blue), "\nAlice syndrome:\n");
        print_array(alice_syndrome);
        fmt::print(fg(fmt::color::blue), "\nBob corrected bit array:\n");
        print_array(bob_solution);
        if (CFG.ENABLE_PRIVACY_MAINTENANCE)
        {
            fmt::print(fg(fmt::color::blue), "\nAlice bit array after privacy maintenance:\n");
            print_array(alice_bit_array_pm);
            fmt::print(fg(fmt::color::blue), "\nBob bit array after privacy maintenance:\n");
            print_array(bob_bit_array_pm);
        }
        fmt::print(fg(fmt::color::blue), "\n\nIterations performed: {}\n", ldpc_res.decoding_res.iterations_num);
        fmt::print(fg(fmt::color::blue), "Syndromes matched: {}\n", ((ldpc_res.decoding_res.syndromes_match) ? "YES" : "NO"));
        fmt::print(fg(fmt::color::blue), "Keys matched: {}\n", ((ldpc_res.keys_match) ? "YES" : "NO"));
    }

    return ldpc_res;
}

LDPC_result QKD_LDPC_RATE_ADAPT(
    const H_matrix &matrix,
    const std::vector<int> &alice_bit_array, 
    const std::vector<int> &bob_bit_array, 
    const double &QBER,
    const decoding_scaling_factors &scaling_factors,
    const H_matrix_params &matrix_params,
    XoshiroCpp::Xoshiro256PlusPlus &prng
)
{
    size_t num_bit_nodes = matrix.bit_nodes.size();
    size_t num_check_nodes = matrix.check_nodes.size();

    double log_p = log((1. - QBER) / QBER);         
    std::vector<double> apriori_llr(num_bit_nodes);

    std::vector<int> alice_bit_array_extended(num_bit_nodes);
    std::vector<int> bob_bit_array_extended(num_bit_nodes);

    std::uniform_int_distribution<int> distribution(0, 1);
    
    size_t p = 0;   // p + s + n = num_bit_nodes
    size_t s = 0;
    size_t n = 0;
    size_t num_punct_bits = matrix_params.punctured_bits.size();
    size_t num_short_bits = matrix_params.shortened_bits.size();

    for (int i = 0; i < num_bit_nodes; ++i)
    {
        if (p < num_punct_bits && matrix_params.punctured_bits[p] == i)
        {
            // Punctured bits come from RNG, independently of the both sides.
            alice_bit_array_extended[i] = distribution(prng);
            bob_bit_array_extended[i] = distribution(prng);
            apriori_llr[i] = ALMOST_ZERO;   // To avoid division by zero.
            ++p;
        }
        else if (s < num_short_bits && matrix_params.shortened_bits[s] == i)
        {
            // Shortened symbols are the ones which have values exactly
            // known by Alice and Bob.
            alice_bit_array_extended[i] = 0;
            bob_bit_array_extended[i] = 0;
            apriori_llr[i] = std::numeric_limits<double>::max();    // +inf
            ++s;
        }
        else
        {
            alice_bit_array_extended[i] = alice_bit_array[n];
            bob_bit_array_extended[i] = bob_bit_array[n];
            apriori_llr[i] = (bob_bit_array[n] ? -log_p : log_p);
            ++n;
        }
    }
    // std::vector<int> alice_remaining_bits_buffer(p + s);
    // std::vector<int> bob_remaining_bits_buffer(p + s);
    // alice_remaining_bits_buffer.assign(alice_bit_array.begin() + n, alice_bit_array.end());
    // bob_remaining_bits_buffer.assign(bob_bit_array.begin() + n, bob_bit_array.end());

    std::vector<int> alice_syndrome(num_check_nodes);
    calculate_syndrome(alice_bit_array_extended, matrix, alice_syndrome);

    std::vector<int> bob_solution(num_bit_nodes);
    LDPC_result ldpc_res;
    if (CFG.DECODING_ALGORITHM == DEC_SPA)
    {
        ldpc_res.decoding_res = sum_product_decoding(apriori_llr, matrix, alice_syndrome, CFG.DECODING_ALG_MAX_ITERATIONS,
            CFG.DECODING_ALG_MSG_LLR_THRESHOLD, bob_solution);
    }
    else if (CFG.DECODING_ALGORITHM == DEC_SPA_APPROX)
    {
        ldpc_res.decoding_res = sum_product_linear_approx_decoding(apriori_llr, matrix, alice_syndrome, CFG.DECODING_ALG_MAX_ITERATIONS,
            CFG.DECODING_ALG_MSG_LLR_THRESHOLD, bob_solution);
    }
    else if (CFG.DECODING_ALGORITHM == DEC_NMSA)
    {
        ldpc_res.decoding_res = min_sum_normalized_decoding(apriori_llr, matrix, alice_syndrome, CFG.DECODING_ALG_MAX_ITERATIONS, 
            scaling_factors.primary, CFG.DECODING_ALG_MSG_LLR_THRESHOLD, bob_solution);
    }
    else if (CFG.DECODING_ALGORITHM == DEC_OMSA)
    {
        ldpc_res.decoding_res = min_sum_offset_decoding(apriori_llr, matrix, alice_syndrome, CFG.DECODING_ALG_MAX_ITERATIONS, 
            scaling_factors.primary, CFG.DECODING_ALG_MSG_LLR_THRESHOLD, bob_solution);
    }
    else if (CFG.DECODING_ALGORITHM == DEC_ANMSA)
    {
        ldpc_res.decoding_res = adaptive_min_sum_normalized_decoding(apriori_llr, matrix, alice_syndrome, CFG.DECODING_ALG_MAX_ITERATIONS, 
            scaling_factors.primary, scaling_factors.secondary, CFG.DECODING_ALG_MSG_LLR_THRESHOLD, bob_solution);
    }
    else if (CFG.DECODING_ALGORITHM == DEC_AOMSA)
    {
        ldpc_res.decoding_res = adaptive_min_sum_offset_decoding(apriori_llr, matrix, alice_syndrome, CFG.DECODING_ALG_MAX_ITERATIONS, 
            scaling_factors.primary, scaling_factors.secondary, CFG.DECODING_ALG_MSG_LLR_THRESHOLD, bob_solution);
    }

    ldpc_res.keys_match = arrays_equal(alice_bit_array_extended, bob_solution);

    std::vector<int> alice_bit_array_rb;
    std::vector<int> bob_bit_array_rb;
    remove_bits(matrix_params.bits_to_remove, alice_bit_array_extended, bob_solution, alice_bit_array_rb, bob_bit_array_rb);

    if (CFG.TRACE_QKD_LDPC)
    {
        fmt::print(fg(fmt::color::blue), "\nAlice bit array:\n");
        print_array(alice_bit_array);
        fmt::print(fg(fmt::color::blue), "\nBob bit array with errors:\n");
        print_array(bob_bit_array);
        fmt::print(fg(fmt::color::blue), "\nAlice extended bit array:\n");
        print_array(alice_bit_array_extended);
        fmt::print(fg(fmt::color::blue), "\nBob extended bit array:\n");
        print_array(bob_bit_array_extended);
        fmt::print(fg(fmt::color::blue), "\nr:\n");
        print_array(apriori_llr);
        fmt::print(fg(fmt::color::blue), "\nAlice syndrome:\n");
        print_array(alice_syndrome);
        fmt::print(fg(fmt::color::blue), "\nBob corrected bit array:\n");
        print_array(bob_solution);
        if (CFG.ENABLE_PRIVACY_MAINTENANCE)
        {
            fmt::print(fg(fmt::color::blue), "\nAlice bit array after privacy maintenance:\n");
            print_array(alice_bit_array_rb);
            fmt::print(fg(fmt::color::blue), "\nBob bit array after privacy maintenance:\n");
            print_array(bob_bit_array_rb);
        }
        else
        {
            fmt::print(fg(fmt::color::blue), "\nAlice bit array after removing punctured and shortened bits:\n");
            print_array(alice_bit_array_rb);
            fmt::print(fg(fmt::color::blue), "\nBob bit array after removing punctured and shortened bits:\n");
            print_array(bob_bit_array_rb);
        }
        fmt::print(fg(fmt::color::blue), "\n\nIterations performed: {}\n", ldpc_res.decoding_res.iterations_num);
        fmt::print(fg(fmt::color::blue), "Syndromes matched: {}\n", ((ldpc_res.decoding_res.syndromes_match) ? "YES" : "NO"));
        fmt::print(fg(fmt::color::blue), "Keys matched: {}\n", ((ldpc_res.keys_match) ? "YES" : "NO"));
    }

    return ldpc_res;
}
