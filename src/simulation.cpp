#include "simulation.hpp"

// Records the results of the simulation in a ".csv" format file.
fs::path write_file(
    const std::vector<sim_result> &data,
    std::string sim_duration,
    fs::path directory
)
{
    // Custom locale settings
    // class custom_numpunct : public std::numpunct<char> 
    // {
    //     protected:
    //         char do_decimal_point() const override 
    //         {
    //             return ','; 
    //         }

    //         std::string do_grouping() const override 
    //         {
    //             return ""; 
    //         }
    // };

    try
    {
        if (!fs::exists(directory))
            fs::create_directories(directory);

        bool is_nmsa_omsa = false, is_anmsa_aomsa = false;
        std::string dec_alg_name, scaling_factor_name = "";
        if (CFG.DECODING_ALGORITHM == DEC_SPA)
        { 
            dec_alg_name = "SPA";
        }
        else if (CFG.DECODING_ALGORITHM == DEC_SPA_APPROX)
        {
            dec_alg_name = "SPA-LIN-APPROX";
        }
        else if (CFG.DECODING_ALGORITHM == DEC_NMSA) 
        {
            dec_alg_name = "NMSA";
            scaling_factor_name = ";ALPHA";
            is_nmsa_omsa = true;
        }
        else if (CFG.DECODING_ALGORITHM == DEC_OMSA)
        {
            dec_alg_name = "OMSA";
            scaling_factor_name = ";BETA";
            is_nmsa_omsa = true;
        }
        else if (CFG.DECODING_ALGORITHM == DEC_ANMSA) 
        {
            dec_alg_name = "ANMSA";
            scaling_factor_name = ";ALPHA;NU";
            is_anmsa_aomsa = true;
        }
        else if (CFG.DECODING_ALGORITHM == DEC_AOMSA)
        {
            dec_alg_name = "AOMSA";
            scaling_factor_name = ";BETA;SIGMA";
            is_anmsa_aomsa = true;
        }

        std::string rate_adapt;
        if (CFG.ENABLE_CODE_RATE_ADAPTATION)
        {
            rate_adapt = "ON[punct=";
            if (CFG.ENABLE_UNTAINTED_PUNCTURING)
                rate_adapt += "untainted]";
            else
                rate_adapt += "random]";
        }
        else
            rate_adapt = "OFF";
        
        std::string rtt_part = "";
        if (CFG.ENABLE_THROUGHPUT_MEASUREMENT && CFG.CONSIDER_RTT) 
            rtt_part = ",RTT=" + fmt::format("{:.3f}", CFG.RTT) + "ms";

        std::string base_filename =
            "ldpc("
            "trial_num=" + std::to_string(CFG.TRIALS_NUMBER) + "," 
            "dec_alg=" + dec_alg_name + "," 
            "max_dec_alg_iters=" + std::to_string(CFG.DECODING_ALG_MAX_ITERATIONS) + "," 
            "priv_maint=" + std::string(CFG.ENABLE_PRIVACY_MAINTENANCE ? "ON" : "OFF") + "," 
            "rate_adapt=" + rate_adapt + 
            rtt_part + ","
            "seed=" + std::to_string(CFG.SIMULATION_SEED) + "," 
            "sim_duration=" + sim_duration + 
            ")";

        std::string extension = ".csv";
        fs::path result_file_path = directory / (base_filename + extension);

        size_t file_count = 1;
        while (fs::exists(result_file_path))
        {
            result_file_path = directory / (base_filename + "_" + std::to_string(file_count) + extension);
            file_count++;
        }

        std::fstream fout;
        // std::locale custom_locale(std::locale(""), new custom_numpunct());
        // fout.imbue(custom_locale);
        
        fout.open(result_file_path, std::ios::out | std::ios::trunc);
        fout << "#;MATRIX_FILENAME;TYPE;R;M;N;QBER;ITER_SUCCESS_MEAN;ITER_SUCCESS_STD;ITER_SUCCESS_MIN;ITER_SUCCESS_MAX;"
        "RATIO_SUCCESS_DEC;RATIO_SUCCESS_LDPC;FER";
        if (CFG.ENABLE_CODE_RATE_ADAPTATION)
            fout << ";DELTA;EFFICIENCY;PUNCT_FRACTION;SHORT_FRACTION;R_ADAPTED";
        if (CFG.ENABLE_THROUGHPUT_MEASUREMENT)
            fout << ";THROUGHPUT_MEAN;THROUGHPUT_STD;THROUGHPUT_MIN;THROUGHPUT_MAX";
        fout << scaling_factor_name << "\n";

        for (size_t i = 0; i < data.size(); i++)
        {
            double FER = 1. - data[i].ratio_trials_success_ldpc;
            FER = (trunc(FER * static_cast<double>(CFG.TRIALS_NUMBER)) / static_cast<double>(CFG.TRIALS_NUMBER));
            
            std::string line = fmt::format(
            "{};{};{};{:.3f};{};{};{:.4f};{:.2f};{:.2f};{};{};{};{};{}",
            data[i].sim_number,
            data[i].matrix_filename,
            (data[i].is_regular ? "regular" : "irregular"),
            1. - (static_cast<double>(data[i].num_check_nodes) / static_cast<double>(data[i].num_bit_nodes)),
            data[i].num_check_nodes,
            data[i].num_bit_nodes,
            data[i].accurate_QBER,
            data[i].iter_success_dec_alg_mean,
            data[i].iter_success_dec_alg_std_dev,
            data[i].iter_success_dec_alg_min,
            data[i].iter_success_dec_alg_max,
            data[i].ratio_trials_success_dec_alg,
            data[i].ratio_trials_success_ldpc,
            FER
            );

            if (CFG.ENABLE_CODE_RATE_ADAPTATION) 
            {
                line += fmt::format(
                    ";{:.3f};{:.3f};{:.3f};{:.3f};{:.3f}",
                    data[i].delta,
                    data[i].efficiency,
                    data[i].punctured_fraction,
                    data[i].shortened_fraction,
                    data[i].adapted_code_rate
                );
            }

            if (CFG.ENABLE_THROUGHPUT_MEASUREMENT) 
            {
                line += fmt::format(
                    ";{};{};{};{}",
                    data[i].throughput_mean,
                    data[i].throughput_std_dev,
                    data[i].throughput_min,
                    data[i].throughput_max
                );
            }

            if (is_nmsa_omsa || is_anmsa_aomsa) 
                line += fmt::format(";{:.3f}", data[i].scaling_factors.primary);
            if (is_anmsa_aomsa)
                line += fmt::format(";{:.3f}", data[i].scaling_factors.secondary);
            fout << line << "\n";
        }
        fout.close();
        return result_file_path;
    }
    catch (const std::exception &ex)
    {
        fmt::print(stderr, fg(fmt::color::red), "An error occurred while writing to the file.\n");
        throw;
    }
}

// Get QBER range based on code rate of matrix. 
// R_QBER_maps must be sorted. Looks for the first set of parameters
// where the code rate is less than or equal to the specified rate, 
// and uses these parameters to generate a range of QBER values.
std::vector<double> get_rate_based_QBER_range(
    const double code_rate,
    const std::vector<R_QBER_map> &R_QBER_maps
)
{
    std::vector<double> QBER;
    for (size_t i = 0; i < R_QBER_maps.size(); i++)
    {
        if (code_rate <= R_QBER_maps[i].code_rate)
        {
            if (R_QBER_maps[i].QBER_begin == R_QBER_maps[i].QBER_end)   // Use only one specified value.
            {
                QBER.push_back(R_QBER_maps[i].QBER_begin);
                break;
            }
            
            size_t steps = static_cast<size_t>(round((R_QBER_maps[i].QBER_end - R_QBER_maps[i].QBER_begin) / R_QBER_maps[i].QBER_step)) + 1;    // including 'end' value 
            double value {};
            for (size_t j = 0; j < steps; j++) 
            {
                value = R_QBER_maps[i].QBER_begin + static_cast<double>(j) * R_QBER_maps[i].QBER_step;
                QBER.push_back(value);
            }
            break;
        }
    }
    if (QBER.empty())
    {
        throw std::runtime_error("An error occurred while generating a QBER range based on code rate(R).");
    }
    return QBER;
}

// Get delta(δ) and efficiency(f_EC) ranges based on code rate of matrix. 
// R_adapt_param_maps must be sorted. Looks for the first set of parameters
// where the code rate is less than or equal to the specified rate, 
// and uses these parameters to generate a range of delta and efficiency values.
R_adaptation_parameters_range get_rate_based_adapt_parameter_range(
    const double code_rate,
    const std::vector<R_adaptation_parameters_map> &R_adapt_param_maps
)
{
    std::vector<double> delta{};
    for (size_t i = 0; i < R_adapt_param_maps.size(); i++)
    {
        if (code_rate <= R_adapt_param_maps[i].code_rate)
        {
            if (R_adapt_param_maps[i].delta_begin == R_adapt_param_maps[i].delta_end)   // Use only one specified value.
            {
                delta.push_back(R_adapt_param_maps[i].delta_begin);
                break;
            }
            
            size_t steps = static_cast<size_t>(round((R_adapt_param_maps[i].delta_end - R_adapt_param_maps[i].delta_begin) 
            / R_adapt_param_maps[i].delta_step)) + 1;    // including 'end' value 
            double value {};
            for (size_t j = 0; j < steps; j++) 
            {
                value = R_adapt_param_maps[i].delta_begin + static_cast<double>(j) * R_adapt_param_maps[i].delta_step;
                delta.push_back(value);
            }
            break;
        }
    }
    if (delta.empty())
    {
        throw std::runtime_error("An error occurred while generating a delta range based on code rate(R).");
    }

    std::vector<double> efficiency;
    for (size_t i = 0; i < R_adapt_param_maps.size(); i++)
    {
        if (code_rate <= R_adapt_param_maps[i].code_rate)
        {
            if (R_adapt_param_maps[i].efficiency_begin == R_adapt_param_maps[i].efficiency_end)   // Use only one specified value.
            {
                efficiency.push_back(R_adapt_param_maps[i].efficiency_begin);
                break;
            }
            
            size_t steps = static_cast<size_t>(round((R_adapt_param_maps[i].efficiency_end - R_adapt_param_maps[i].efficiency_begin) 
            / R_adapt_param_maps[i].efficiency_step)) + 1;    // including 'end' value 
            double value {};
            for (size_t j = 0; j < steps; j++) 
            {
                value = R_adapt_param_maps[i].efficiency_begin + static_cast<double>(j) * R_adapt_param_maps[i].efficiency_step;
                efficiency.push_back(value);
            }
            break;
        }
    }
    if (efficiency.empty())
    {
        throw std::runtime_error("An error occurred while generating a efficiency(f_EC) range based on code rate(R).");
    }

    return {.delta = delta, .efficiency = efficiency};
}

// Get all scaling factor range values used for all matrices regardless 
// of their code rate(R).
std::vector<double> get_scaling_factor_range_values(const scaling_factor_range &scaling_factor_range)
{
    std::vector<double> scaling_factors;
    if (scaling_factor_range.begin == scaling_factor_range.end)   // Use only one specified value.
        scaling_factors.push_back(scaling_factor_range.begin);
    else
    {
        size_t steps = static_cast<size_t>(round((scaling_factor_range.end - scaling_factor_range.begin) / scaling_factor_range.step)) + 1;   // including 'end' value 
        for (size_t i = 0; i < steps; i++) 
        {
            scaling_factors.push_back(scaling_factor_range.begin + static_cast<double>(i) * scaling_factor_range.step);
        }
    }

    if (scaling_factors.empty())
        throw std::runtime_error("An error occurred while generating vector of scaling factor values.");

    return scaling_factors;
}

// Get scaling factor value based on code rate of matrix. 
// R_scaling_factor_maps must be sorted. Looks for the first of parameters
// where the code rate is less than or equal to the specified rate.
double get_rate_based_scaling_factor_value(
    const double code_rate,
    const std::vector<R_scaling_factor_map> &R_scaling_factor_maps
)
{
    int param_idx = -1;
    for (int i = 0; i < R_scaling_factor_maps.size(); i++)
    {
        if (code_rate <= R_scaling_factor_maps[i].code_rate)
        {
            param_idx = i;
            break;
        }
    }
    if (param_idx == -1)
        throw std::runtime_error("An error occurred while searching scaling factor value on the basis of code rate(R).");

    return R_scaling_factor_maps[param_idx].scaling_factor;
}

// Prepares input data for batch simulation.
std::vector<sim_input> prepare_sim_inputs(const std::vector<fs::path> &matrix_paths)
{
    XoshiroCpp::Xoshiro256PlusPlus prng(CFG.SIMULATION_SEED);
    std::vector<sim_input> sim_inputs(matrix_paths.size());
    for (size_t i = 0; i < matrix_paths.size(); i++)
    {
        // Reading the matrix depending on the format
        if (CFG.MATRIX_FORMAT == MAT_UNCOMPRESSED)
            sim_inputs[i].matrix = read_uncompressed_matrix(matrix_paths[i]);
        else if (CFG.MATRIX_FORMAT == MAT_SPARSE_ALIST)
            sim_inputs[i].matrix = read_sparse_matrix_alist(matrix_paths[i]);
        else if (CFG.MATRIX_FORMAT == MAT_SPARSE_1)
            sim_inputs[i].matrix = read_sparse_matrix_1(matrix_paths[i]);
        else if (CFG.MATRIX_FORMAT == MAT_SPARSE_2)
            sim_inputs[i].matrix = read_sparse_matrix_2(matrix_paths[i]);

        sim_inputs[i].matrix_path = matrix_paths[i];

        double code_rate = 1. - static_cast<double>(sim_inputs[i].matrix.check_nodes.size()) / static_cast<double>(sim_inputs[i].matrix.bit_nodes.size());
        std::vector<double> QBER_values = get_rate_based_QBER_range(code_rate, CFG.R_QBER_MAPS);

        // List of pairs (QBER, H_matrix_params)
        std::vector<std::pair<double, H_matrix_params>> qber_mat_params;

        if(CFG.ENABLE_CODE_RATE_ADAPTATION)
        {
            R_adaptation_parameters_range rate_adapt_range = get_rate_based_adapt_parameter_range(code_rate, CFG.R_ADAPT_PARAMS_MAPS);
            if (CFG.ENABLE_UNTAINTED_PUNCTURING)
                sim_inputs[i].matrix.punctured_bits_untainted = get_punctured_bits_untainted(matrix_paths[i], prng, sim_inputs[i].matrix);
                
            for (double QBER : QBER_values)
            {
                for (size_t j = 0; j < rate_adapt_range.delta.size(); j++)
                {
                    double delta = rate_adapt_range.delta[j];
                    for (size_t k = 0; k < rate_adapt_range.efficiency.size(); k++)
                    {
                        double efficiency = rate_adapt_range.efficiency[k];
                        H_matrix_params mat_params = adapt_code_rate(prng, sim_inputs[i].matrix, QBER, delta, efficiency);
                        // Skip parameters combination. It will not be used in simulations.
                        if (mat_params.punctured_bits.empty() && mat_params.shortened_bits.empty())
                            continue;

                        if (CFG.ENABLE_PRIVACY_MAINTENANCE) 
                            mat_params.bits_to_remove = get_bits_positions_to_remove_rate_adapt(sim_inputs[i].matrix, mat_params);
                        else 
                        {
                            mat_params.bits_to_remove.reserve(mat_params.punctured_bits.size() + mat_params.shortened_bits.size());
                            std::merge(mat_params.punctured_bits.begin(), mat_params.punctured_bits.end(),
                                       mat_params.shortened_bits.begin(), mat_params.shortened_bits.end(),
                                       std::back_inserter(mat_params.bits_to_remove));
                        }
                        qber_mat_params.emplace_back(QBER, mat_params);
                    }
                }
            }
        }
        else
        {
            H_matrix_params mat_params{};
            if (CFG.ENABLE_PRIVACY_MAINTENANCE)
                mat_params.bits_to_remove = get_bits_positions_to_remove(sim_inputs[i].matrix);
            for (double QBER : QBER_values) 
            {
                qber_mat_params.emplace_back(QBER, mat_params);
            }
        }
        
        std::vector<decoding_scaling_factors> scaling_combinations{};
        std::vector<double> primary_values{};
        if (CFG.DECODING_ALGORITHM == DEC_NMSA || CFG.DECODING_ALGORITHM == DEC_OMSA)
        {
            if (CFG.DECODING_ALG_PARAMS.primary.use_range)
                primary_values = get_scaling_factor_range_values(CFG.DECODING_ALG_PARAMS.primary.range);
            else
            {
                double scaling_factor = get_rate_based_scaling_factor_value(code_rate, CFG.DECODING_ALG_PARAMS.primary.maps);
                primary_values.push_back(scaling_factor);      // Contains only one value.
            }
            for (double primary : primary_values) 
            {
                decoding_scaling_factors sf{};
                sf.primary = primary;
                scaling_combinations.push_back(sf);
            }
        }
        else if (CFG.DECODING_ALGORITHM == DEC_ANMSA || CFG.DECODING_ALGORITHM == DEC_AOMSA)
        {
            if (CFG.DECODING_ALG_PARAMS.primary.use_range)
                primary_values = get_scaling_factor_range_values(CFG.DECODING_ALG_PARAMS.primary.range);
            else 
            {
                double scaling_factor = get_rate_based_scaling_factor_value(code_rate, CFG.DECODING_ALG_PARAMS.primary.maps);
                primary_values.push_back(scaling_factor);
            }

            std::vector<double> secondary_values;
            if (CFG.DECODING_ALG_PARAMS.secondary.use_range)
                secondary_values = get_scaling_factor_range_values(CFG.DECODING_ALG_PARAMS.secondary.range);
            else 
            {
                double scaling_factor = get_rate_based_scaling_factor_value(code_rate, CFG.DECODING_ALG_PARAMS.secondary.maps);
                secondary_values.push_back(scaling_factor);
            }
            for (double primary : primary_values) 
            {
                for (double secondary : secondary_values) 
                {
                    decoding_scaling_factors sf;
                    sf.primary = primary;
                    sf.secondary = secondary;
                    scaling_combinations.push_back(sf);
                }
            }
        }
        else    // For algorithms without scaling_factors
        {
            decoding_scaling_factors sf{};
            scaling_combinations.push_back(sf);
        }

        // Generation of all combinations
        sim_inputs[i].combinations.reserve(qber_mat_params.size() * scaling_combinations.size());
        for (const auto& qmp : qber_mat_params) 
        {
            for (const auto& sf : scaling_combinations) 
            {
                sim_combination comb;
                comb.QBER = qmp.first;
                comb.matrix_params = qmp.second;
                comb.scaling_factors = sf;
                sim_inputs[i].combinations.push_back(comb);
            }
        }
    }
    return sim_inputs;
}

// Runs a single QKD LDPC trial.
trial_result run_trial(
    const H_matrix &matrix, 
    double QBER, 
    size_t seed,
    const H_matrix_params &matrix_params,
    const decoding_scaling_factors &scaling_factors
)
{
    trial_result result;
    XoshiroCpp::Xoshiro256PlusPlus prng(seed);
    size_t num_bit_nodes = matrix.bit_nodes.size();
    std::vector<int> alice_bit_array(num_bit_nodes);
    std::vector<int> bob_bit_array(num_bit_nodes);

    fill_random_bits(prng, alice_bit_array);
    result.accurate_QBER = inject_errors(prng, alice_bit_array, QBER, bob_bit_array);
    if (result.accurate_QBER == 0.)
        throw std::runtime_error("Key size '" + std::to_string(num_bit_nodes) + "' is too small for QBER.");

    if (CFG.ENABLE_THROUGHPUT_MEASUREMENT)
    {
        auto start = std::chrono::high_resolution_clock::now();
        if (CFG.ENABLE_CODE_RATE_ADAPTATION)
            result.ldpc_res = QKD_LDPC_RATE_ADAPT(matrix, alice_bit_array, bob_bit_array, result.accurate_QBER, scaling_factors, matrix_params, prng);
        else
            result.ldpc_res = QKD_LDPC(matrix, alice_bit_array, bob_bit_array, result.accurate_QBER, scaling_factors, matrix_params);
        auto end = std::chrono::high_resolution_clock::now();
        result.runtime = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    }
    else
    {
        if (CFG.ENABLE_CODE_RATE_ADAPTATION)
            result.ldpc_res = QKD_LDPC_RATE_ADAPT(matrix, alice_bit_array, bob_bit_array, result.accurate_QBER, scaling_factors, matrix_params, prng);
        else
            result.ldpc_res = QKD_LDPC(matrix, alice_bit_array, bob_bit_array, result.accurate_QBER, scaling_factors, matrix_params);
    }
    return result;
}

// Calculation of statistical characteristics and recording of the simulation result.
void process_trials_results(
    const std::vector<trial_result> &trial_results, 
    const H_matrix &matrix,
    const H_matrix_params &matrix_params,
    sim_result &result
)
{
    size_t trials_successful_decoding = 0;
    size_t trials_successful_ldpc = 0;
    size_t iter_success_dec_alg_max = 0;
    size_t iter_success_dec_alg_min = std::numeric_limits<size_t>::max();
    double iter_success_dec_alg_mean = 0;
    double iter_success_dec_alg_std_dev = 0;   //standard deviation
    size_t curr_decoding_iterations_num = 0;
    for (size_t i = 0; i < trial_results.size(); ++i)
    {
        if (trial_results[i].ldpc_res.decoding_res.syndromes_match)
        {
            trials_successful_decoding++;
            curr_decoding_iterations_num = trial_results[i].ldpc_res.decoding_res.iterations_num;
            if (iter_success_dec_alg_max < curr_decoding_iterations_num)
                iter_success_dec_alg_max = curr_decoding_iterations_num;
            if (iter_success_dec_alg_min > curr_decoding_iterations_num)
                iter_success_dec_alg_min = curr_decoding_iterations_num;
            if (trial_results[i].ldpc_res.keys_match)
                trials_successful_ldpc++;

            iter_success_dec_alg_mean += static_cast<double>(curr_decoding_iterations_num);
        }
    }

    if (trials_successful_decoding > 0)
    {
        iter_success_dec_alg_mean /= static_cast<double>(trials_successful_decoding);
        for (size_t i = 0; i < trial_results.size(); ++i)
        {
            if (trial_results[i].ldpc_res.decoding_res.syndromes_match)
            {
                curr_decoding_iterations_num = trial_results[i].ldpc_res.decoding_res.iterations_num;
                iter_success_dec_alg_std_dev += pow((static_cast<double>(curr_decoding_iterations_num) - iter_success_dec_alg_mean), 2);
            }
        }
        iter_success_dec_alg_std_dev /= static_cast<double>(trials_successful_decoding);
        iter_success_dec_alg_std_dev = sqrt(iter_success_dec_alg_std_dev);
    }

    if (CFG.ENABLE_THROUGHPUT_MEASUREMENT)
    {
        size_t num_bit_nodes = matrix.bit_nodes.size();
        size_t num_bits_to_remove = matrix_params.bits_to_remove.size();
        double out_key_length{};
        if (CFG.ENABLE_CODE_RATE_ADAPTATION || CFG.ENABLE_PRIVACY_MAINTENANCE)
            out_key_length = static_cast<double>(num_bit_nodes - num_bits_to_remove);
        else
            out_key_length = static_cast<double>(num_bit_nodes);
        
        const double MICROSECONDS_IN_SECOND = 1000000.;
        const double MICROSECONDS_IN_MILLISECOND = 1000.;
        double curr_throughput{};
        double throughput_max = 0;
        double throughput_min = std::numeric_limits<double>::max();
        double throughput_mean = 0;
        double throughput_std_dev = 0;

        for (size_t i = 0; i < trial_results.size(); ++i)
        {
            if (CFG.CONSIDER_RTT)
            {
                curr_throughput = out_key_length * MICROSECONDS_IN_SECOND / 
                (static_cast<double>(trial_results[i].runtime.count()) + CFG.RTT * MICROSECONDS_IN_MILLISECOND);  // bits/s
            }
            else
                curr_throughput = out_key_length * MICROSECONDS_IN_SECOND / static_cast<double>(trial_results[i].runtime.count());  // bits/s
            
            throughput_mean += curr_throughput;
            if (curr_throughput > throughput_max)
                throughput_max = curr_throughput;
            if (curr_throughput < throughput_min)
                throughput_min = curr_throughput;
        }
        throughput_mean /= static_cast<double>(CFG.TRIALS_NUMBER);
        
        for (size_t i = 0; i < trial_results.size(); ++i)
        {
            if (CFG.CONSIDER_RTT)
            {
                curr_throughput = out_key_length * MICROSECONDS_IN_SECOND / 
                (static_cast<double>(trial_results[i].runtime.count()) + CFG.RTT * MICROSECONDS_IN_MILLISECOND);  // bits/s
            }
            else
                curr_throughput = out_key_length * MICROSECONDS_IN_SECOND / static_cast<double>(trial_results[i].runtime.count());  // bits/s

            throughput_std_dev += pow((curr_throughput - throughput_mean), 2);
        }
        throughput_std_dev /= static_cast<double>(CFG.TRIALS_NUMBER);
        throughput_std_dev = sqrt(throughput_std_dev);

        result.throughput_max = static_cast<size_t>(throughput_max);
        result.throughput_min = static_cast<size_t>(throughput_min);
        result.throughput_mean = static_cast<size_t>(throughput_mean);
        result.throughput_std_dev = static_cast<size_t>(throughput_std_dev);        
    }
    
    result.iter_success_dec_alg_max = iter_success_dec_alg_max;
    result.iter_success_dec_alg_min = (iter_success_dec_alg_min == std::numeric_limits<size_t>::max()) ? 0 : static_cast<size_t>(iter_success_dec_alg_min);
    result.iter_success_dec_alg_mean = iter_success_dec_alg_mean;
    result.iter_success_dec_alg_std_dev = iter_success_dec_alg_std_dev;

    result.ratio_trials_success_ldpc = static_cast<double>(trials_successful_ldpc) / static_cast<double>(CFG.TRIALS_NUMBER);
    result.ratio_trials_success_dec_alg = static_cast<double>(trials_successful_decoding) / static_cast<double>(CFG.TRIALS_NUMBER);            
}

// Distributes all combinations of the experiment evenly across the CPU threads and runs it.
std::vector<sim_result> QKD_LDPC_batch_simulation(const std::vector<sim_input> &sim_in)
{
    using namespace indicators;
    size_t sim_total = 0;
    for (size_t i = 0; i < sim_in.size(); ++i) 
    {
        sim_total += sim_in[i].combinations.size();
    }
    size_t trials_total = sim_total * CFG.TRIALS_NUMBER;

    indicators::show_console_cursor(false);
    indicators::ProgressBar bar{
        option::BarWidth{50}, option::Start{" ["}, option::Fill{"="}, option::Lead{">"},
        option::Remainder{"-"}, option::End{"]"}, option::PrefixText{"PROGRESS"},
        option::ForegroundColor{Color::cyan}, option::ShowElapsedTime{true},
        option::ShowRemainingTime{true}, option::FontStyles{std::vector<FontStyle>{FontStyle::bold}},
        option::MaxProgress{trials_total}};

    std::vector<sim_result> sim_results(sim_total);
    std::vector<trial_result> trial_results(CFG.TRIALS_NUMBER);
    XoshiroCpp::Xoshiro256PlusPlus prng(CFG.SIMULATION_SEED);
    std::uniform_int_distribution<size_t> distribution(0, std::numeric_limits<size_t>::max());
    std::vector<size_t> seeds(CFG.TRIALS_NUMBER);
    for (size_t i = 0; i < seeds.size(); ++i) 
    {
        seeds[i] = distribution(prng);
    }

    BS::thread_pool pool(CFG.THREADS_NUMBER);
    size_t curr_sim = 0;
    size_t iteration = 0;

    for (size_t i = 0; i < sim_in.size(); ++i) 
    {
        const H_matrix &matrix = sim_in[i].matrix;
        std::string matrix_filename = sim_in[i].matrix_path.filename().string();

        for (auto& comb : sim_in[i].combinations) 
        {
            double QBER = comb.QBER;
            H_matrix_params matrix_params = comb.matrix_params;
            decoding_scaling_factors scaling_factors = comb.scaling_factors;

            iteration += CFG.TRIALS_NUMBER;
            bar.set_option(option::PostfixText{
                std::to_string(iteration) + "/" + std::to_string(trials_total)});

            pool.detach_loop<size_t>(0, CFG.TRIALS_NUMBER,
                [&matrix, &QBER, &matrix_params, &scaling_factors, &trial_results, &seeds, &curr_sim, &bar](size_t n) 
                {
                    trial_results[n] = run_trial(matrix, QBER, (seeds[n] + curr_sim), matrix_params, scaling_factors);
                    bar.tick();
                });
            pool.wait();

            sim_results[curr_sim].sim_number = curr_sim;
            sim_results[curr_sim].matrix_filename = matrix_filename;
            sim_results[curr_sim].is_regular = matrix.is_regular;
            sim_results[curr_sim].num_bit_nodes = matrix.bit_nodes.size();
            sim_results[curr_sim].num_check_nodes = matrix.check_nodes.size();
            sim_results[curr_sim].delta = matrix_params.delta;
            sim_results[curr_sim].efficiency = matrix_params.efficiency;
            sim_results[curr_sim].punctured_fraction = matrix_params.punctured_fraction;
            sim_results[curr_sim].shortened_fraction = matrix_params.shortened_fraction;
            sim_results[curr_sim].adapted_code_rate = matrix_params.adapted_code_rate;
            sim_results[curr_sim].accurate_QBER = trial_results[0].accurate_QBER;
            sim_results[curr_sim].scaling_factors = scaling_factors;

            process_trials_results(trial_results, matrix, matrix_params, sim_results[curr_sim]);
            ++curr_sim;
        }
    }
    indicators::show_console_cursor(true);
    return sim_results;
}
