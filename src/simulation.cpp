#include "simulation.hpp"

// Records the results of the simulation in a ".csv" format file
void write_file(const std::vector<sim_result> &data,
                fs::path directory)
{
    try
    {
        if (!fs::exists(directory))
        {
            fs::create_directories(directory);
        }
        std::string base_filename  = "ldpc(trial_num=" + std::to_string(CFG.TRIALS_NUMBER) + ",decoding_alg=" + 
                               ((CFG.USE_MIN_SUM_DECODING_ALG)?"MSA":"SPA") + ",max_decoding_alg_iters=" +
                               std::to_string(CFG.DECODING_ALG_MAX_ITERATIONS) + ",seed=" + std::to_string(CFG.SIMULATION_SEED) + ")";
        std::string extension = ".csv";
        fs::path result_file_path = directory / (base_filename + extension);

        size_t file_count = 1;
        while (fs::exists(result_file_path))
        {
            result_file_path = directory / (base_filename + "_" + std::to_string(file_count) + extension);
            file_count++;
        }

        std::fstream fout;
        fout.open(result_file_path, std::ios::out | std::ios::trunc);
        fout << "№;MATRIX_FILENAME;TYPE;CODE_RATE;M;N;QBER;ITERATIONS_SUCCESSFUL_DEC_ALG_MEAN;ITERATIONS_SUCCESSFUL_DEC_ALG_STD_DEV;ITERATIONS_SUCCESSFUL_DEC_ALG_MIN;ITERATIONS_SUCCESSFUL_DEC_ALG_MAX;" << 
        "RATIO_TRIALS_SUCCESSFUL_DEC_ALG;RATIO_TRIALS_SUCCESSFUL_LDPC;FER\n";
        for (size_t i = 0; i < data.size(); i++)
        {
            fout << data[i].sim_number << ";" << data[i].matrix_filename << ";" << (data[i].is_regular ? "regular" : "irregular") << ";" 
                 << 1. - (static_cast<double>(data[i].num_check_nodes) / data[i].num_bit_nodes) << ";" << data[i].num_check_nodes << ";" 
                 << data[i].num_bit_nodes << ";" << data[i].initial_QBER << ";" << data[i].iter_success_dec_alg_mean << ";" 
                 << data[i].iter_success_dec_alg_std_dev << ";" << data[i].iter_success_dec_alg_min << ";" 
                 << data[i].iter_success_dec_alg_max << ";" << data[i].ratio_trials_success_dec_alg << ";" 
                 << data[i].ratio_trials_success_ldpc << ";" << 1. - data[i].ratio_trials_success_ldpc << "\n";
        }
        fout.close();
    }
    catch (const std::exception &ex)
    {
        fmt::print(stderr, fg(fmt::color::red), "An error occurred while writing to the file.\n");
        throw;
    }
}

// Get QBER range based on code rate of matrix. R_QBER_parameters must be sorted. Looks for the first set of parameters
// where the code rate is less than or equal to the specified rate, and uses these parameters to generate a range of QBER values.
std::vector<double> get_rate_based_QBER_range(const double code_rate,
                                              const std::vector<R_QBER_params> &R_QBER_parameters)
{
    std::vector<double> QBER;
    for (size_t i = 0; i < R_QBER_parameters.size(); i++)
    {
        if (code_rate <= R_QBER_parameters[i].code_rate)
        {
            size_t steps = round((R_QBER_parameters[i].QBER_end - R_QBER_parameters[i].QBER_begin) / R_QBER_parameters[i].QBER_step);
            double value {};
            for (size_t j = 0; j < steps; j++) 
            {
                value = R_QBER_parameters[i].QBER_begin + j * R_QBER_parameters[i].QBER_step;
                QBER.push_back(value);
            }
            break;
        }
    }
    if (QBER.empty())
    {
        throw std::runtime_error("An error occurred when generating a QBER range based on code rate.");
    }
    return QBER;
}

// Interactive simulation of quantum key distribution (QKD) using LDPC codes.
void QKD_LDPC_interactive_simulation(fs::path matrix_dir_path)
{
    H_matrix matrix;
    std::vector<fs::path> matrix_paths = get_file_paths_in_directory(matrix_dir_path);
    fs::path matrix_path = select_matrix_file(matrix_paths);

    if (CFG.MATRIX_FORMAT == 0)
        matrix = read_dense_matrix(matrix_path);
    else if (CFG.MATRIX_FORMAT == 1)
        matrix = read_sparse_alist_matrix(matrix_path);
    else if (CFG.MATRIX_FORMAT == 2)
        matrix = read_sparse_matrix(matrix_path);

    fmt::print(fg(fmt::color::green), "{}\n", ((matrix.is_regular) ? "Matrix H is regular." : "Matrix H is irregular."));

    size_t num_check_nodes = matrix.check_nodes.size();
    size_t num_bit_nodes = matrix.bit_nodes.size();
    std::vector<int> alice_bit_array(num_bit_nodes);
    std::vector<int> bob_bit_array(num_bit_nodes);

    XoshiroCpp::Xoshiro256PlusPlus prng(CFG.SIMULATION_SEED); 

    double code_rate = 1. - (static_cast<double>(num_check_nodes) / num_bit_nodes);
    std::vector<double> QBER = get_rate_based_QBER_range(code_rate, CFG.R_QBER_PARAMETERS);

    for (size_t i = 0; i < QBER.size(); ++i)
    {
        fmt::print(fg(fmt::color::green), "№:{}\n", i + 1);

        fill_random_bits(prng, alice_bit_array);
        double initial_QBER = introduce_errors(prng, alice_bit_array, QBER[i], bob_bit_array);
        fmt::print(fg(fmt::color::green), "Actual QBER: {}\n", initial_QBER);

        if (initial_QBER == 0.)
        {
            throw std::runtime_error("Key size '" + std::to_string(num_bit_nodes) + "' is too small for QBER.");
        }

        int error_num = 0;
        for (size_t i = 0; i < num_bit_nodes; ++i)
        {
            error_num += alice_bit_array[i] ^ bob_bit_array[i];
        }
        fmt::print(fg(fmt::color::green), "Number of errors in a key: {}\n", error_num);

        LDPC_result try_result;
        try_result = QKD_LDPC(alice_bit_array, bob_bit_array, initial_QBER, matrix);
        
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
        if (CFG.MATRIX_FORMAT == 0)
            sim_inputs[i].matrix = read_dense_matrix(matrix_paths[i]);
        else if (CFG.MATRIX_FORMAT == 1)
            sim_inputs[i].matrix = read_sparse_alist_matrix(matrix_paths[i]);
        else if (CFG.MATRIX_FORMAT == 2)
            sim_inputs[i].matrix = read_sparse_matrix(matrix_paths[i]);

        sim_inputs[i].matrix_path = matrix_paths[i];

        double code_rate = 1. - (static_cast<double>(sim_inputs[i].matrix.check_nodes.size()) / sim_inputs[i].matrix.bit_nodes.size());
        sim_inputs[i].QBER = get_rate_based_QBER_range(code_rate, CFG.R_QBER_PARAMETERS);
    }
    return sim_inputs;
}

// Runs a single QKD LDPC trial.
trial_result run_trial(const H_matrix &matrix, 
                       double QBER, 
                       size_t seed)
{
    trial_result result;
    XoshiroCpp::Xoshiro256PlusPlus prng(seed);
    size_t num_bit_nodes = matrix.bit_nodes.size();
    std::vector<int> alice_bit_array(num_bit_nodes);
    std::vector<int> bob_bit_array(num_bit_nodes);

    fill_random_bits(prng, alice_bit_array);
    result.initial_QBER = introduce_errors(prng, alice_bit_array, QBER, bob_bit_array);
    if (result.initial_QBER == 0.)
    {
        throw std::runtime_error("Key size '" + std::to_string(num_bit_nodes) + "' is too small for QBER.");
    }

    result.ldpc_res = QKD_LDPC(alice_bit_array, bob_bit_array, result.initial_QBER, matrix);

    return result;
}

// Distributes all combinations of the experiment evenly across the CPU threads and runs it.
std::vector<sim_result> QKD_LDPC_batch_simulation(const std::vector<sim_input> &sim_in)
{
    using namespace indicators;
    size_t sim_total = 0;
    for (size_t i = 0; i < sim_in.size(); i++)
    {
        sim_total += sim_in[i].QBER.size();     // For each matrix, keys are generated with error rates given in the QBER vector
    }

    size_t trials_total = sim_total * CFG.TRIALS_NUMBER;    
    indicators::show_console_cursor(false);
    indicators::ProgressBar bar{
        option::BarWidth{50},
        option::Start{" ["},
        option::Fill{"="},
        option::Lead{">"},
        option::Remainder{"-"},
        option::End{"]"},
        option::PrefixText{"PROGRESS"},
        option::ForegroundColor{Color::green},
        option::ShowElapsedTime{true},
        option::ShowRemainingTime{true},
        option::FontStyles{std::vector<FontStyle>{FontStyle::bold}},
        option::MaxProgress{trials_total}};

    size_t curr_sim = 0;
    size_t iteration = 0;
    std::vector<sim_result> sim_results(sim_total);
    std::vector<trial_result> trial_results(CFG.TRIALS_NUMBER);

    XoshiroCpp::Xoshiro256PlusPlus prng(CFG.SIMULATION_SEED);
    std::uniform_int_distribution<size_t> distribution(0, std::numeric_limits<size_t>::max());
    std::vector<size_t> seeds(CFG.TRIALS_NUMBER);
    for (size_t i = 0; i < seeds.size(); i++)
    {
        seeds[i] = distribution(prng);
    }

    BS::thread_pool pool(CFG.THREADS_NUMBER);
    for (size_t i = 0; i < sim_in.size(); i++)
    {
        const H_matrix &matrix = sim_in[i].matrix;
        double code_rate = 1. - (static_cast<double>(matrix.check_nodes.size()) / matrix.bit_nodes.size());
        std::string matrix_filename = sim_in[i].matrix_path.filename().string();
        for (size_t j = 0; j < sim_in[i].QBER.size(); j++)
        {
            iteration += CFG.TRIALS_NUMBER;
            bar.set_option(option::PostfixText{
                std::to_string(iteration) + "/" + std::to_string(trials_total)});

            double QBER = sim_in[i].QBER[j];
            // TRIALS_NUMBER of trials are performed with each combination to calculate the mean values
            pool.detach_loop<size_t>(0, CFG.TRIALS_NUMBER,
                                     [&matrix, &QBER, &trial_results, &seeds, &curr_sim, &bar](size_t k)
                                     {
                                         trial_results[k] = run_trial(matrix, QBER, (seeds[k] + curr_sim));
                                         bar.tick(); // For correct time estimation
                                     });
            pool.wait();

            size_t trials_successful_decoding = 0;
            size_t trials_successful_ldpc = 0;
            size_t iter_success_dec_alg_max = 0;
            size_t iter_success_dec_alg_min = CFG.DECODING_ALG_MAX_ITERATIONS;
            double iter_success_dec_alg_mean = 0;
            double iter_success_dec_alg_std_dev = 0;   //standard deviation
            size_t curr_decoding_iterations_num = 0;
            for (size_t k = 0; k < trial_results.size(); k++)
            {
                if (trial_results[k].ldpc_res.decoding_res.syndromes_match)
                {
                    trials_successful_decoding++;
                    curr_decoding_iterations_num = trial_results[k].ldpc_res.decoding_res.iterations_num;
                    if (iter_success_dec_alg_max < curr_decoding_iterations_num)
                    {
                        iter_success_dec_alg_max = curr_decoding_iterations_num;
                    }
                    if (iter_success_dec_alg_min > curr_decoding_iterations_num)
                    {
                        iter_success_dec_alg_min = curr_decoding_iterations_num;
                    }
                    if (trial_results[k].ldpc_res.keys_match)
                    {
                        trials_successful_ldpc++;
                    }

                    iter_success_dec_alg_mean += static_cast<double>(curr_decoding_iterations_num);
                }
            }

            if (trials_successful_decoding > 0)
            {
                iter_success_dec_alg_mean /= static_cast<double>(trials_successful_decoding);
                for (size_t k = 0; k < trial_results.size(); k++)
                {
                    if (trial_results[k].ldpc_res.decoding_res.syndromes_match)
                    {
                        curr_decoding_iterations_num = trial_results[k].ldpc_res.decoding_res.iterations_num;
                        iter_success_dec_alg_std_dev += pow((static_cast<double>(curr_decoding_iterations_num) - iter_success_dec_alg_mean), 2);
                    }
                }
                iter_success_dec_alg_std_dev /= static_cast<double>(trials_successful_decoding);
                iter_success_dec_alg_std_dev = sqrt(iter_success_dec_alg_std_dev);
            }

            sim_results[curr_sim].sim_number = curr_sim;

            sim_results[curr_sim].matrix_filename = matrix_filename;
            sim_results[curr_sim].is_regular = matrix.is_regular;
            sim_results[curr_sim].num_bit_nodes = matrix.bit_nodes.size();
            sim_results[curr_sim].num_check_nodes = matrix.check_nodes.size();

            sim_results[curr_sim].initial_QBER = trial_results[0].initial_QBER;
            sim_results[curr_sim].iter_success_dec_alg_max = iter_success_dec_alg_max;
            sim_results[curr_sim].iter_success_dec_alg_min = (iter_success_dec_alg_min == CFG.DECODING_ALG_MAX_ITERATIONS) ? 0. : iter_success_dec_alg_min;
            sim_results[curr_sim].iter_success_dec_alg_mean = iter_success_dec_alg_mean;
            sim_results[curr_sim].iter_success_dec_alg_std_dev = iter_success_dec_alg_std_dev;

            sim_results[curr_sim].ratio_trials_success_ldpc = static_cast<double>(trials_successful_ldpc) / CFG.TRIALS_NUMBER;
            sim_results[curr_sim].ratio_trials_success_dec_alg = static_cast<double>(trials_successful_decoding) / CFG.TRIALS_NUMBER;
            curr_sim++;
        }
    }
    return sim_results;
}
