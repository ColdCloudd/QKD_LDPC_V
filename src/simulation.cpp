#include "simulation.hpp"

// Records the results of the simulation in a ".csv" format file
fs::path write_file(const std::vector<sim_result> &data,
                    std::string sim_duration,
                    fs::path directory)
{
    // Custom locale settings
    class custom_numpunct : public std::numpunct<char> 
    {
        protected:
            char do_decimal_point() const override 
            {
                return ','; 
            }

            std::string do_grouping() const override 
            {
                return ""; 
            }
    };

    try
    {
        if (!fs::exists(directory))
            fs::create_directories(directory);

        bool is_nmsa_omsa = false, is_anmsa_aomsa = false;
        std::string dec_alg_name, scaling_factor_name = "";
        if (CFG.DECODING_ALGORITHM == SPA) 
            dec_alg_name = "SPA";
        else if (CFG.DECODING_ALGORITHM == SPA_APPROX) 
            dec_alg_name = "SPA-LIN-APPROX";
        else if (CFG.DECODING_ALGORITHM == NMSA) 
        {
            dec_alg_name = "NMSA";
            scaling_factor_name = ";ALPHA";
            is_nmsa_omsa = true;
        }
        else if (CFG.DECODING_ALGORITHM == OMSA)
        {
            dec_alg_name = "OMSA";
            scaling_factor_name = ";BETA";
            is_nmsa_omsa = true;
        }
        else if (CFG.DECODING_ALGORITHM == ANMSA) 
        {
            dec_alg_name = "ANMSA";
            scaling_factor_name = ";ALPHA;NU";
            is_anmsa_aomsa = true;
        }
        else if (CFG.DECODING_ALGORITHM == AOMSA)
        {
            dec_alg_name = "AOMSA";
            scaling_factor_name = ";BETA;SIGMA";
            is_anmsa_aomsa = true;
        }
        
        std::string base_filename  = "ldpc(trial_num=" + std::to_string(CFG.TRIALS_NUMBER) + ",dec_alg=" + 
                               dec_alg_name + ",max_dec_alg_iters=" +
                               std::to_string(CFG.DECODING_ALG_MAX_ITERATIONS) + ",priv_maint=" + 
                               ((CFG.ENABLE_PRIVACY_MAINTENANCE)?"on":"off") + 
                               ((CFG.ENABLE_THROUGHPUT_MEASUREMENT && CFG.CONSIDER_RTT)?(",RTT=" + std::to_string(CFG.RTT)):"") + 
                               ",seed=" + std::to_string(CFG.SIMULATION_SEED) + ",sim_duration=" + sim_duration + ")";

        std::string extension = ".csv";
        fs::path result_file_path = directory / (base_filename + extension);

        size_t file_count = 1;
        while (fs::exists(result_file_path))
        {
            result_file_path = directory / (base_filename + "_" + std::to_string(file_count) + extension);
            file_count++;
        }

        std::locale custom_locale(std::locale("ru_RU.utf8"), new custom_numpunct());
        std::fstream fout;
        fout.imbue(custom_locale);
        
        fout.open(result_file_path, std::ios::out | std::ios::trunc);
        fout << "#;MATRIX_FILENAME;TYPE;CODE_RATE;M;N;QBER;ITERATIONS_SUCCESSFUL_DEC_ALG_MEAN;ITERATIONS_SUCCESSFUL_DEC_ALG_STD_DEV;ITERATIONS_SUCCESSFUL_DEC_ALG_MIN;ITERATIONS_SUCCESSFUL_DEC_ALG_MAX;" << 
        "RATIO_TRIALS_SUCCESSFUL_DEC_ALG;RATIO_TRIALS_SUCCESSFUL_LDPC;FER" << (CFG.ENABLE_THROUGHPUT_MEASUREMENT ? ";THROUGHPUT_MEAN;THROUGHPUT_STD_DEV;THROUGHPUT_MIN;THROUGHPUT_MAX" : "") << 
        scaling_factor_name << "\n";

        std::string prim_scal_fact_val = "", sec_scal_fact_val = "";
        for (size_t i = 0; i < data.size(); i++)
        {
            if (is_nmsa_omsa)
                prim_scal_fact_val = (";" + std::to_string(data[i].scaling_factors.primary));
            else if (is_anmsa_aomsa)
            {
                prim_scal_fact_val = (";" + std::to_string(data[i].scaling_factors.primary));
                sec_scal_fact_val = (";" + std::to_string(data[i].scaling_factors.secondary));
            }
            
            fout << data[i].sim_number << ";" << data[i].matrix_filename << ";" << (data[i].is_regular ? "regular" : "irregular") << ";" 
                 << 1. - (static_cast<double>(data[i].num_check_nodes) / data[i].num_bit_nodes) << ";" << data[i].num_check_nodes << ";" 
                 << data[i].num_bit_nodes << ";" << data[i].initial_QBER << ";" << data[i].iter_success_dec_alg_mean << ";" 
                 << data[i].iter_success_dec_alg_std_dev << ";" << data[i].iter_success_dec_alg_min << ";" 
                 << data[i].iter_success_dec_alg_max << ";" << data[i].ratio_trials_success_dec_alg << ";" 
                 << data[i].ratio_trials_success_ldpc << ";" << 1. - data[i].ratio_trials_success_ldpc 
                 << (CFG.ENABLE_THROUGHPUT_MEASUREMENT ? (";" + std::to_string(data[i].throughput_mean) + ";"
                 + std::to_string(data[i].throughput_std_dev) + ";" + std::to_string(data[i].throughput_min) + ";"
                 + std::to_string(data[i].throughput_max)):"") << prim_scal_fact_val << sec_scal_fact_val << "\n";
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

// Get QBER range based on code rate of matrix. R_QBER_maps must be sorted. Looks for the first set of parameters
// where the code rate is less than or equal to the specified rate, and uses these parameters to generate a range of QBER values.
std::vector<double> get_rate_based_QBER_range(const double code_rate,
                                              const std::vector<R_QBER_map> &R_QBER_maps)
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

// Get all scaling factor range values used for all matrices regardless of their code rate(R).
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

// Get scaling factor value based on code rate of matrix. R_scaling_factor_maps must be sorted. Looks for the first of parameters
// where the code rate is less than or equal to the specified rate.
double get_rate_based_scaling_factor_value(const double code_rate,
                                           const std::vector<R_scaling_factor_map> &R_scaling_factor_maps)
{
    double scaling_factor = -1.;
    for (size_t i = 0; i < R_scaling_factor_maps.size(); i++)
    {
        if (code_rate <= R_scaling_factor_maps[i].code_rate)
        {
            scaling_factor = R_scaling_factor_maps[i].scaling_factor;
            break;
        }
    }
    if (scaling_factor == -1.)
        throw std::runtime_error("An error occurred while searching scaling factor value on the basis of code rate(R).");

    return scaling_factor;
}

// Interactive simulation of quantum key distribution (QKD) using LDPC codes.
void QKD_LDPC_interactive_simulation(fs::path matrix_dir_path)
{
    H_matrix matrix;
    std::vector<fs::path> matrix_paths = get_file_paths_in_directory(matrix_dir_path);
    fs::path matrix_path = select_matrix_file(matrix_paths);

    if (CFG.MATRIX_FORMAT == DENSE_MAT)
        matrix = read_dense_matrix(matrix_path);
    else if (CFG.MATRIX_FORMAT == SPARSE_MAT_ALIST)
        matrix = read_sparse_matrix_alist(matrix_path);
    else if (CFG.MATRIX_FORMAT == SPARSE_MAT_1)
        matrix = read_sparse_matrix_1(matrix_path);
    else if (CFG.MATRIX_FORMAT == SPARSE_MAT_2)
        matrix = read_sparse_matrix_2(matrix_path);

    fmt::print(fg(fmt::color::green), "{}\n", ((matrix.is_regular) ? "Matrix H is regular." : "Matrix H is irregular."));

    size_t num_check_nodes = matrix.check_nodes.size();
    size_t num_bit_nodes = matrix.bit_nodes.size();
    std::vector<int> alice_bit_array(num_bit_nodes);
    std::vector<int> bob_bit_array(num_bit_nodes);

    XoshiroCpp::Xoshiro256PlusPlus prng(CFG.SIMULATION_SEED); 

    double code_rate = 1. - static_cast<double>(num_check_nodes) / static_cast<double>(num_bit_nodes);
    std::vector<double> QBER = get_rate_based_QBER_range(code_rate, CFG.R_QBER_MAPS);
    
    decoding_scaling_factors scaling_factors{};
    if (CFG.DECODING_ALGORITHM == SPA)
    {
        fmt::print(fg(fmt::color::green), "Sum-product decoding algorithm selected.\n\n");
    }
    else if (CFG.DECODING_ALGORITHM == SPA_APPROX)
    {
        fmt::print(fg(fmt::color::green), "Sum-product decoding algorithm with linear approximation of tanh and atanh selected.\n\n");
    }
    else
    {
        if (CFG.DECODING_ALG_PARAMS.primary.use_range || CFG.DECODING_ALG_PARAMS.secondary.use_range)
            throw std::runtime_error("Scaling factor ranges are not available in the interactive simulation. Use code rate and scaling factor correspondences!");
        
        if (CFG.DECODING_ALGORITHM == NMSA)
        {
            scaling_factors.primary = get_rate_based_scaling_factor_value(code_rate, CFG.DECODING_ALG_PARAMS.primary.maps);
            fmt::print(fg(fmt::color::green), "Min-sum normalized decoding algorithm selected. Alpha = {}\n\n", std::to_string(scaling_factors.primary));
        }
        else if (CFG.DECODING_ALGORITHM == OMSA)
        {
            scaling_factors.primary = get_rate_based_scaling_factor_value(code_rate, CFG.DECODING_ALG_PARAMS.primary.maps);
            fmt::print(fg(fmt::color::green), "Min-sum offset decoding algorithm selected. Beta = {}\n\n", std::to_string(scaling_factors.primary));
        }
        else if (CFG.DECODING_ALGORITHM == ANMSA)
        {
            scaling_factors.primary = get_rate_based_scaling_factor_value(code_rate, CFG.DECODING_ALG_PARAMS.primary.maps);
            scaling_factors.secondary = get_rate_based_scaling_factor_value(code_rate, CFG.DECODING_ALG_PARAMS.secondary.maps);
            fmt::print(fg(fmt::color::green), "Adaptive min-sum normalized decoding algorithm selected. Alpha = {}, Nu = {}\n\n", 
                std::to_string(scaling_factors.primary), std::to_string(scaling_factors.secondary));
        }
        else if (CFG.DECODING_ALGORITHM == AOMSA)
        {
            scaling_factors.primary = get_rate_based_scaling_factor_value(code_rate, CFG.DECODING_ALG_PARAMS.primary.maps);
            scaling_factors.secondary = get_rate_based_scaling_factor_value(code_rate, CFG.DECODING_ALG_PARAMS.secondary.maps);
            fmt::print(fg(fmt::color::green), "Adaptive min-sum offset decoding algorithm selected. Beta = {}, Sigma = {}\n\n", 
                std::to_string(scaling_factors.primary), std::to_string(scaling_factors.secondary));
        }
    }

    for (size_t i = 0; i < QBER.size(); ++i)
    {
        fmt::print(fg(fmt::color::green), "№:{}\n", i + 1);

        fill_random_bits(prng, alice_bit_array);
        double initial_QBER = introduce_errors(prng, alice_bit_array, QBER[i], bob_bit_array);
        fmt::print(fg(fmt::color::green), "Actual QBER: {}\n", initial_QBER);

        if (initial_QBER == 0.)
            throw std::runtime_error("Key size '" + std::to_string(num_bit_nodes) + "' is too small for QBER.");

        int error_num = 0;
        for (size_t j = 0; j < num_bit_nodes; ++j)
        {
            error_num += alice_bit_array[j] ^ bob_bit_array[j];
        }
        fmt::print(fg(fmt::color::green), "Number of errors in a key: {}\n", error_num);

        LDPC_result try_result;
        if (CFG.DECODING_ALGORITHM > SPA_APPROX)
            try_result = QKD_LDPC(matrix, alice_bit_array, bob_bit_array, initial_QBER, scaling_factors);
        else 
            try_result = QKD_LDPC(matrix, alice_bit_array, bob_bit_array, initial_QBER);

        fmt::print(fg(fmt::color::green), "Iterations performed: {}\n", try_result.decoding_res.iterations_num);
        fmt::print(fg(fmt::color::green), "{}\n\n", ((try_result.keys_match && try_result.decoding_res.syndromes_match) ? "Error reconciliation SUCCESSFUL" : "Error reconciliation FAILED"));
    }
}

// Prepares input data for batch simulation.
std::vector<sim_input> prepare_sim_inputs(const std::vector<fs::path> &matrix_paths)
{
    std::vector<sim_input> sim_inputs(matrix_paths.size());
    for (size_t i = 0; i < matrix_paths.size(); i++)
    {
        if (CFG.MATRIX_FORMAT == DENSE_MAT)
            sim_inputs[i].matrix = read_dense_matrix(matrix_paths[i]);
        else if (CFG.MATRIX_FORMAT == SPARSE_MAT_ALIST)
            sim_inputs[i].matrix = read_sparse_matrix_alist(matrix_paths[i]);
        else if (CFG.MATRIX_FORMAT == SPARSE_MAT_1)
            sim_inputs[i].matrix = read_sparse_matrix_1(matrix_paths[i]);
        else if (CFG.MATRIX_FORMAT == SPARSE_MAT_2)
            sim_inputs[i].matrix = read_sparse_matrix_2(matrix_paths[i]);

        sim_inputs[i].matrix_path = matrix_paths[i];

        double code_rate = 1. - static_cast<double>(sim_inputs[i].matrix.check_nodes.size()) / static_cast<double>(sim_inputs[i].matrix.bit_nodes.size());
        sim_inputs[i].QBER = get_rate_based_QBER_range(code_rate, CFG.R_QBER_MAPS);
        
        if (CFG.DECODING_ALGORITHM == NMSA || CFG.DECODING_ALGORITHM == OMSA)
        {
            if (CFG.DECODING_ALG_PARAMS.primary.use_range)
                sim_inputs[i].primary_scaling_factor = get_scaling_factor_range_values(CFG.DECODING_ALG_PARAMS.primary.range);
            else
            {
                double scaling_factor = get_rate_based_scaling_factor_value(code_rate, CFG.DECODING_ALG_PARAMS.primary.maps);
                sim_inputs[i].primary_scaling_factor.push_back(scaling_factor);   // Contains only one value.
            }
        }
        else if (CFG.DECODING_ALGORITHM == ANMSA || CFG.DECODING_ALGORITHM == AOMSA)
        {
            double scaling_factor;
            if (CFG.DECODING_ALG_PARAMS.primary.use_range && CFG.DECODING_ALG_PARAMS.secondary.use_range)           // Use only ranges.
            {
                sim_inputs[i].primary_scaling_factor = get_scaling_factor_range_values(CFG.DECODING_ALG_PARAMS.primary.range);
                sim_inputs[i].secondary_scaling_factor = get_scaling_factor_range_values(CFG.DECODING_ALG_PARAMS.secondary.range);
            }
            else if (CFG.DECODING_ALG_PARAMS.primary.use_range && !CFG.DECODING_ALG_PARAMS.secondary.use_range)     // Primary scaling factor - range, secondary - map.
            {
                sim_inputs[i].primary_scaling_factor = get_scaling_factor_range_values(CFG.DECODING_ALG_PARAMS.primary.range);
                scaling_factor = get_rate_based_scaling_factor_value(code_rate, CFG.DECODING_ALG_PARAMS.secondary.maps);
                sim_inputs[i].secondary_scaling_factor.push_back(scaling_factor); 
            }
            else if (!CFG.DECODING_ALG_PARAMS.primary.use_range && CFG.DECODING_ALG_PARAMS.secondary.use_range)     // Primary scaling factor - map, secondary - range.
            {
                scaling_factor = get_rate_based_scaling_factor_value(code_rate, CFG.DECODING_ALG_PARAMS.primary.maps);
                sim_inputs[i].primary_scaling_factor.push_back(scaling_factor);
                sim_inputs[i].secondary_scaling_factor = get_scaling_factor_range_values(CFG.DECODING_ALG_PARAMS.secondary.range);
            }
            else    // Use only maps.
            {
                scaling_factor = get_rate_based_scaling_factor_value(code_rate, CFG.DECODING_ALG_PARAMS.primary.maps);
                sim_inputs[i].primary_scaling_factor.push_back(scaling_factor);
                scaling_factor = get_rate_based_scaling_factor_value(code_rate, CFG.DECODING_ALG_PARAMS.secondary.maps);
                sim_inputs[i].secondary_scaling_factor.push_back(scaling_factor);   
            }
        }
    }
    return sim_inputs;
}

// Runs a single QKD LDPC trial.
trial_result run_trial(const H_matrix &matrix, 
                       double QBER, 
                       size_t seed,
                       const decoding_scaling_factors &scaling_factors)
{
    trial_result result;
    XoshiroCpp::Xoshiro256PlusPlus prng(seed);
    size_t num_bit_nodes = matrix.bit_nodes.size();
    std::vector<int> alice_bit_array(num_bit_nodes);
    std::vector<int> bob_bit_array(num_bit_nodes);

    fill_random_bits(prng, alice_bit_array);
    result.initial_QBER = introduce_errors(prng, alice_bit_array, QBER, bob_bit_array);
    if (result.initial_QBER == 0.)
        throw std::runtime_error("Key size '" + std::to_string(num_bit_nodes) + "' is too small for QBER.");

    if (CFG.ENABLE_THROUGHPUT_MEASUREMENT)
    {
        auto start = std::chrono::high_resolution_clock::now();
        result.ldpc_res = QKD_LDPC(matrix, alice_bit_array, bob_bit_array, result.initial_QBER, scaling_factors);
        auto end = std::chrono::high_resolution_clock::now();
        result.runtime = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    }
    else
        result.ldpc_res = QKD_LDPC(matrix, alice_bit_array, bob_bit_array, result.initial_QBER, scaling_factors);

    return result;
}

// Calculation of statistical characteristics and recording of the simulation result
void process_trials_results(const std::vector<trial_result> &trial_results, 
                            const size_t &num_bit_nodes, 
                            const size_t &num_check_nodes,
                            sim_result &result)
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
        double out_key_length = static_cast<double>(num_bit_nodes - num_check_nodes);
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
                (static_cast<double>(trial_results[i].runtime.count()) + static_cast<double>(CFG.RTT) * MICROSECONDS_IN_MILLISECOND);  // bits/s
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
                (static_cast<double>(trial_results[i].runtime.count()) + static_cast<double>(CFG.RTT) * MICROSECONDS_IN_MILLISECOND);  // bits/s
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

    const bool is_nmsa_omsa = (CFG.DECODING_ALGORITHM == NMSA || CFG.DECODING_ALGORITHM == OMSA);
    const bool is_anmsa_aomsa = (CFG.DECODING_ALGORITHM == ANMSA || CFG.DECODING_ALGORITHM == AOMSA);

    for (size_t i = 0; i < sim_in.size(); ++i)
    {     
        if (is_nmsa_omsa)
            sim_total += sim_in[i].QBER.size() * sim_in[i].primary_scaling_factor.size();    // For each QBER value, tests are performed with all alpha (or beta) values from the vector
        else if (is_anmsa_aomsa)
            sim_total += sim_in[i].QBER.size() * sim_in[i].primary_scaling_factor.size() * sim_in[i].secondary_scaling_factor.size();     // Each scaling factor (alpha or beta) from primary with each scaling factor (nu or sigma) from secondary 
        else
            sim_total += sim_in[i].QBER.size();     // For each matrix, keys are generated with error rates given in the QBER vector
    }

    size_t trials_total = sim_total * CFG.TRIALS_NUMBER;    
    indicators::show_console_cursor(false);
    indicators::ProgressBar bar{
        option::BarWidth{50}, option::Start{" ["}, option::Fill{"="}, option::Lead{">"},
        option::Remainder{"-"}, option::End{"]"}, option::PrefixText{"PROGRESS"},
        option::ForegroundColor{Color::cyan}, option::ShowElapsedTime{true},
        option::ShowRemainingTime{true}, option::FontStyles{std::vector<FontStyle>{FontStyle::bold}},
        option::MaxProgress{trials_total}};

    size_t curr_sim = 0;
    size_t iteration = 0;
    decoding_scaling_factors scaling_factors{};
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
    for (size_t i = 0; i < sim_in.size(); ++i)
    {
        const H_matrix &matrix = sim_in[i].matrix;
        const size_t num_bit_nodes = matrix.bit_nodes.size();
        const size_t num_check_nodes = matrix.check_nodes.size();
        const std::string matrix_filename = sim_in[i].matrix_path.filename().string();

        for (size_t j = 0; j < sim_in[i].QBER.size(); ++j)
        {
            double QBER = sim_in[i].QBER[j];
            if (is_nmsa_omsa)
            {
                for (size_t k = 0; k < sim_in[i].primary_scaling_factor.size(); ++k)
                {
                    scaling_factors.primary = sim_in[i].primary_scaling_factor[k];

                    iteration += CFG.TRIALS_NUMBER;
                    bar.set_option(option::PostfixText{
                        std::to_string(iteration) + "/" + std::to_string(trials_total)});

                    // TRIALS_NUMBER of trials are performed with each combination to calculate the mean values
                    pool.detach_loop<size_t>(0, CFG.TRIALS_NUMBER,
                                            [&matrix, &QBER, &scaling_factors, &trial_results, &seeds, &curr_sim, &bar](size_t n)
                                            {
                                                trial_results[n] = run_trial(matrix, QBER, (seeds[n] + curr_sim), scaling_factors);
                                                bar.tick(); // For correct time estimation
                                            });
                    pool.wait();

                    sim_results[curr_sim].sim_number = curr_sim;
                    sim_results[curr_sim].matrix_filename = matrix_filename;
                    sim_results[curr_sim].is_regular = matrix.is_regular;
                    sim_results[curr_sim].num_bit_nodes = num_bit_nodes;
                    sim_results[curr_sim].num_check_nodes = num_check_nodes;
                    sim_results[curr_sim].initial_QBER = trial_results[0].initial_QBER;
                    sim_results[curr_sim].scaling_factors = scaling_factors;

                    process_trials_results(trial_results, num_bit_nodes, num_check_nodes, sim_results[curr_sim]);
                    ++curr_sim;
                }
            }
            else if (is_anmsa_aomsa)
            {
                for (size_t k = 0; k < sim_in[i].primary_scaling_factor.size(); ++k)
                {
                    scaling_factors.primary = sim_in[i].primary_scaling_factor[k];
                    for (size_t l = 0; l < sim_in[i].secondary_scaling_factor.size(); ++l)
                    {
                        scaling_factors.secondary = sim_in[i].secondary_scaling_factor[l];

                        iteration += CFG.TRIALS_NUMBER;
                        bar.set_option(option::PostfixText{
                            std::to_string(iteration) + "/" + std::to_string(trials_total)});

                        // TRIALS_NUMBER of trials are performed with each combination to calculate the mean values
                        pool.detach_loop<size_t>(0, CFG.TRIALS_NUMBER,
                                                [&matrix, &QBER, &scaling_factors, &trial_results, &seeds, &curr_sim, &bar](size_t n)
                                                {
                                                    trial_results[n] = run_trial(matrix, QBER, (seeds[n] + curr_sim), scaling_factors);
                                                    bar.tick(); // For correct time estimation
                                                });
                        pool.wait();

                        sim_results[curr_sim].sim_number = curr_sim;
                        sim_results[curr_sim].matrix_filename = matrix_filename;
                        sim_results[curr_sim].is_regular = matrix.is_regular;
                        sim_results[curr_sim].num_bit_nodes = num_bit_nodes;
                        sim_results[curr_sim].num_check_nodes = num_check_nodes;
                        sim_results[curr_sim].initial_QBER = trial_results[0].initial_QBER;
                        sim_results[curr_sim].scaling_factors = scaling_factors;

                        process_trials_results(trial_results, num_bit_nodes, num_check_nodes, sim_results[curr_sim]);
                        ++curr_sim;
                    }
                }
            }
            else
            {
                iteration += CFG.TRIALS_NUMBER;
                bar.set_option(option::PostfixText{
                    std::to_string(iteration) + "/" + std::to_string(trials_total)});

                // TRIALS_NUMBER of trials are performed with each combination to calculate the mean values
                pool.detach_loop<size_t>(0, CFG.TRIALS_NUMBER,
                                        [&matrix, &QBER, &trial_results, &seeds, &curr_sim, &bar](size_t n)
                                        {
                                            trial_results[n] = run_trial(matrix, QBER, (seeds[n] + curr_sim));
                                            bar.tick(); // For correct time estimation (move outside this loop to speed up execution time)
                                        });
                pool.wait();

                sim_results[curr_sim].sim_number = curr_sim;
                sim_results[curr_sim].matrix_filename = matrix_filename;
                sim_results[curr_sim].is_regular = matrix.is_regular;
                sim_results[curr_sim].num_bit_nodes = num_bit_nodes;
                sim_results[curr_sim].num_check_nodes = num_check_nodes;
                sim_results[curr_sim].initial_QBER = trial_results[0].initial_QBER;

                process_trials_results(trial_results, num_bit_nodes, num_check_nodes, sim_results[curr_sim]);
                ++curr_sim;
            }
        }
    }
    return sim_results;
}
