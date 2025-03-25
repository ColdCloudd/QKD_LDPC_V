#include "config.hpp"

scaling_factor_range parse_scaling_factor_range(const json& scaling_factor_range) 
{
    double begin = scaling_factor_range["begin"].template get<double>();
    double end = scaling_factor_range["end"].template get<double>();
    double step = scaling_factor_range["step"].template get<double>();
    if (begin <= 0. || end <= 0. || step <= 0.)
        throw std::runtime_error("Scaling factor range begin, end, step must be > 0!");
    if (begin > end)
        throw std::runtime_error("Scaling factor range begin cannot be larger than end!");
    if (begin != end)
    {
        if (step - EPSILON > end - begin)
            throw std::runtime_error("Scaling factor range step is too large!");
    }

    return {begin, end, step};
}

std::vector<R_scaling_factor_map> parse_scaling_factor_maps(const json& scaling_factor_maps,
                                                            const std::string& key) 
{
    std::vector<R_scaling_factor_map> maps;
    double code_rate;
    double scaling_factor;
    for (const auto& m : scaling_factor_maps) 
    {
        code_rate = m["code_rate"].template get<double>();
        scaling_factor = m[key].template get<double>();
        if (code_rate <= 0. || code_rate >= 1.)
            throw std::runtime_error("Code rate(R) must be: 0 < R < 1!");
        if (scaling_factor <= 0.)
            throw std::runtime_error("Scaling factor must be > 0!");

        maps.push_back({code_rate, scaling_factor});
    }
    if (maps.empty())
        throw std::runtime_error("Array with code rate(R) and scaling factor maps is empty!");

    std::sort(maps.begin(), maps.end(),
        [](R_scaling_factor_map &a, R_scaling_factor_map &b)
        {
              return (a.code_rate < b.code_rate);
        });

    return maps;
}

void print_config_info(config_data cfg, 
                       std::string cfg_name,
                       size_t cfg_number)
{
    fmt::print(fg(fmt::color::yellow), "------------------------- CONFIG #{} INFO --------------------------\n", cfg_number);
    fmt::print(fg(fmt::color::yellow), "Config name: {}\n", fmt::styled(cfg_name, fg(fmt::color::crimson)));
    fmt::print(fg(fmt::color::yellow), "Simulation MODE: {}\n", fmt::styled(cfg.INTERACTIVE_MODE ? "INTERACTIVE" : "BATCH", fg(fmt::color::crimson)));
    fmt::print(fg(fmt::color::yellow), "Threads number: {}\n", fmt::styled(cfg.THREADS_NUMBER, fg(fmt::color::crimson)));
    fmt::print(fg(fmt::color::yellow), "Trials number: {}\n", fmt::styled(cfg.TRIALS_NUMBER, fg(fmt::color::crimson)));
    fmt::print(fg(fmt::color::yellow), "Simulation seed: {}\n", fmt::styled(cfg.SIMULATION_SEED, fg(fmt::color::crimson)));
    fmt::print(fg(fmt::color::yellow), "Privacy maintenance: {}\n", fmt::styled((cfg.ENABLE_PRIVACY_MAINTENANCE ? "Enabled" : "Disabled"), fg(fmt::color::crimson)));
    fmt::print(fg(fmt::color::yellow), "Throughput measurement: {}\n", fmt::styled((cfg.ENABLE_THROUGHPUT_MEASUREMENT ? ("Enabled, RTT = " + std::to_string(cfg.RTT) + " ms") : "Disabled"), fg(fmt::color::crimson)));

    std::string alg_name =
        (cfg.DECODING_ALGORITHM == 0) ? "SPA" :
        (cfg.DECODING_ALGORITHM == 1) ? "SPA (lin approx)" :
        (cfg.DECODING_ALGORITHM == 2) ? "NMSA" :
        (cfg.DECODING_ALGORITHM == 3) ? "OMSA" :
        (cfg.DECODING_ALGORITHM == 4) ? "ANMSA" :
        (cfg.DECODING_ALGORITHM == 5) ? "AOMSA" : "Unknown";
    fmt::print(fg(fmt::color::yellow), "Decoding algorithm: {}\n", fmt::styled(alg_name, fg(fmt::color::crimson)));
    fmt::print(fg(fmt::color::yellow), "Decoding algorithm maximum iterations: {}\n", fmt::styled(cfg.DECODING_ALG_MAX_ITERATIONS, fg(fmt::color::crimson)));

    std::string mat_format =
        (cfg.MATRIX_FORMAT == 0) ? "Dense" :
        (cfg.MATRIX_FORMAT == 1) ? "Sparse (alist)" :
        (cfg.MATRIX_FORMAT == 2) ? "Sparse (1)" :
        (cfg.MATRIX_FORMAT == 3) ? "Sparse (2)" : "Unknown";
    fmt::print(fg(fmt::color::yellow), "Parity-check matrix format: {}\n", fmt::styled(mat_format, fg(fmt::color::crimson)));
    fmt::print(fg(fmt::color::yellow), "--------------------------------------------------------------------\n");

}

// Reads user-defined configuration parameters from a .json file.
config_data parse_config_data(fs::path config_path)
{
    if (!fs::exists(config_path))
        throw std::runtime_error("Configuration file not found: " + config_path.string());

    if (config_path.extension() != ".json")
        throw std::runtime_error("Configuration file must have a .json extension: " + config_path.string());

    std::ifstream config_file(config_path);
    if (!config_file.is_open())
        throw std::runtime_error("Failed to open configuration file: " + config_path.string());

    json config = json::parse(config_file);
    config_file.close();
    if (config.empty())
        throw std::runtime_error("Configuration file is empty: " + config_path.string());

    try
    {
        config_data cfg{};
        cfg.THREADS_NUMBER = config["threads_number"].template get<size_t>();
        if (cfg.THREADS_NUMBER < 1)
            throw std::runtime_error("Number of threads must be >= 1!");

        cfg.TRIALS_NUMBER = config["trials_number"].template get<size_t>();
        if (cfg.TRIALS_NUMBER < 1)
            throw std::runtime_error("Number of trials must be >= 1!");

        if (config["use_config_simulation_seed"].template get<bool>())
            cfg.SIMULATION_SEED = config["simulation_seed"].template get<size_t>();
        else
            cfg.SIMULATION_SEED = time(nullptr);

        cfg.INTERACTIVE_MODE = config["interactive_mode"].template get<bool>();
        cfg.ENABLE_PRIVACY_MAINTENANCE = config["enable_privacy_maintenance"].template get<bool>();
        cfg.ENABLE_THROUGHPUT_MEASUREMENT = config["enable_throughput_measurement"].template get<bool>();
        if (cfg.ENABLE_THROUGHPUT_MEASUREMENT)
        {
            fmt::print(fg(fmt::color::purple), "WARNING: Throughput measurement is enabled. It is recommended to perform experiments in single-threaded mode.\n");
            
            const auto &tm_params = config["throughput_measurement_parameters"];
            cfg.CONSIDER_RTT = tm_params["consider_RTT"].template get<bool>();
            if (cfg.CONSIDER_RTT)
                cfg.RTT = tm_params["RTT"].template get<size_t>();
        }

        cfg.DECODING_ALGORITHM = config["decoding_algorithm"].template get<size_t>();
        if (cfg.DECODING_ALGORITHM > AOMSA)
            throw std::runtime_error("Only six options are available: \n0 - SPA;\n1 - SPA (with linear approximation of tanh and atanh);\n2 - NMSA;\n3 - OMSA;\n4 - ANMSA;\n5 - AOMSA.");

        if (cfg.DECODING_ALGORITHM > SPA_APPROX)
        {   
            nlohmann::json alg_params;
            if (cfg.DECODING_ALGORITHM == NMSA)
            {
                alg_params = config["min_sum_normalized_parameters"];
                cfg.DECODING_ALG_PARAMS.primary.use_range = alg_params["use_alpha_range"];
                if (cfg.DECODING_ALG_PARAMS.primary.use_range)
                    cfg.DECODING_ALG_PARAMS.primary.range = parse_scaling_factor_range(alg_params["alpha_range"]);
                else
                    cfg.DECODING_ALG_PARAMS.primary.maps = parse_scaling_factor_maps(alg_params["code_rate_alpha_maps"], "alpha");
            }
            else if (cfg.DECODING_ALGORITHM == OMSA)
            {
                alg_params = config["min_sum_offset_parameters"];
                cfg.DECODING_ALG_PARAMS.primary.use_range = alg_params["use_beta_range"];
                if (cfg.DECODING_ALG_PARAMS.primary.use_range)
                    cfg.DECODING_ALG_PARAMS.primary.range = parse_scaling_factor_range(alg_params["beta_range"]);
                else
                    cfg.DECODING_ALG_PARAMS.primary.maps = parse_scaling_factor_maps(alg_params["code_rate_beta_maps"], "beta");
            }
            else 
            {
                if (cfg.DECODING_ALGORITHM == ANMSA)
                {
                    alg_params = config["adaptive_min_sum_normalized_parameters"];
                    cfg.DECODING_ALG_PARAMS.primary.use_range = alg_params["use_alpha_range"];
                    if (cfg.DECODING_ALG_PARAMS.primary.use_range)
                        cfg.DECODING_ALG_PARAMS.primary.range = parse_scaling_factor_range(alg_params["alpha_range"]);
                    else
                        cfg.DECODING_ALG_PARAMS.primary.maps = parse_scaling_factor_maps(alg_params["code_rate_alpha_maps"], "alpha");
                    
                    cfg.DECODING_ALG_PARAMS.secondary.use_range = alg_params["use_nu_range"];
                    if (cfg.DECODING_ALG_PARAMS.secondary.use_range)
                        cfg.DECODING_ALG_PARAMS.secondary.range = parse_scaling_factor_range(alg_params["nu_range"]);
                    else
                        cfg.DECODING_ALG_PARAMS.secondary.maps = parse_scaling_factor_maps(alg_params["code_rate_nu_maps"], "nu");
                }
                else if (cfg.DECODING_ALGORITHM == AOMSA)
                {
                    alg_params = config["adaptive_min_sum_offset_parameters"];
                    cfg.DECODING_ALG_PARAMS.primary.use_range = alg_params["use_beta_range"];
                    if (cfg.DECODING_ALG_PARAMS.primary.use_range)
                        cfg.DECODING_ALG_PARAMS.primary.range = parse_scaling_factor_range(alg_params["beta_range"]);
                    else
                        cfg.DECODING_ALG_PARAMS.primary.maps = parse_scaling_factor_maps(alg_params["code_rate_beta_maps"], "beta");
                    
                    cfg.DECODING_ALG_PARAMS.secondary.use_range = alg_params["use_sigma_range"];
                    if (cfg.DECODING_ALG_PARAMS.secondary.use_range)
                        cfg.DECODING_ALG_PARAMS.secondary.range = parse_scaling_factor_range(alg_params["sigma_range"]);
                    else
                        cfg.DECODING_ALG_PARAMS.secondary.maps = parse_scaling_factor_maps(alg_params["code_rate_sigma_maps"], "sigma");
                }

                if (!(cfg.DECODING_ALG_PARAMS.primary.use_range || cfg.DECODING_ALG_PARAMS.secondary.use_range))    // Both scaling factors from maps.
                {
                    std::string alg_name, prim_scal_factor, second_scal_factor;
                    if (cfg.DECODING_ALGORITHM == ANMSA)
                    {
                        alg_name = "ANMSA";
                        prim_scal_factor = "alpha";
                        second_scal_factor = "nu";
                    }
                    else if (cfg.DECODING_ALGORITHM == AOMSA)
                    {
                        alg_name = "AOMSA";
                        prim_scal_factor = "beta";
                        second_scal_factor = "sigma";
                    }  
                    
                    if (cfg.DECODING_ALG_PARAMS.primary.maps.size() != cfg.DECODING_ALG_PARAMS.secondary.maps.size())
                    {
                        throw std::runtime_error(fmt::format(
                            "{}: The sizes of code_rate_{}_maps and code_rate_{}_maps vectors must match! ({} vs {})",
                            alg_name, prim_scal_factor, second_scal_factor, cfg.DECODING_ALG_PARAMS.primary.maps.size(),
                            cfg.DECODING_ALG_PARAMS.secondary.maps.size()
                        ));
                    }
                    
                    for (size_t i = 0; i < cfg.DECODING_ALG_PARAMS.primary.maps.size(); ++i) 
                    {
                        const auto& primary_map = cfg.DECODING_ALG_PARAMS.primary.maps[i];      // Both sorted by code_rate.
                        const auto& secondary_map = cfg.DECODING_ALG_PARAMS.secondary.maps[i];
                        
                        if (std::abs(primary_map.code_rate - secondary_map.code_rate) > EPSILON) 
                        {
                            throw std::runtime_error(fmt::format( 
                                "{}: Mismatch of code_rate in {} and {} maps: {:.3f} vs {:.3f}\nAll code_rate values, from code_rate_{}_maps must also be in code_rate_{}_maps!",
                                alg_name, prim_scal_factor, second_scal_factor, primary_map.code_rate,
                                secondary_map.code_rate, prim_scal_factor, second_scal_factor
                            ));
                        }
                    }    
                }
            }
        }
        
        cfg.DECODING_ALG_MAX_ITERATIONS = config["decoding_algorithm_max_iterations"].template get<size_t>();
        if (cfg.DECODING_ALG_MAX_ITERATIONS < 1)
            throw std::runtime_error("Minimum number of decoding algorithm iterations must be >= 1!");

        cfg.MATRIX_FORMAT = config["matrix_format"].template get<size_t>();
        if (cfg.MATRIX_FORMAT > SPARSE_MAT_2)
            throw std::runtime_error("Only four options are available: \n0 - dense;\n1 - sparse alist;\n2 - sparse_1;\n3 - sparse_2.");

        cfg.TRACE_QKD_LDPC = config["trace_qkd_ldpc"].template get<bool>();
        cfg.TRACE_DECODING_ALG = config["trace_decoding_algorithm"].template get<bool>();
        cfg.TRACE_DECODING_ALG_LLR = config["trace_decoding_algorithm_llr"].template get<bool>();
        cfg.ENABLE_DECODING_ALG_MSG_LLR_THRESHOLD = config["enable_decoding_algorithm_msg_llr_threshold"].template get<bool>();

        if (cfg.ENABLE_DECODING_ALG_MSG_LLR_THRESHOLD)
        {
            cfg.DECODING_ALG_MSG_LLR_THRESHOLD = config["decoding_algorithm_msg_llr_threshold"].template get<double>();
            if (cfg.DECODING_ALG_MSG_LLR_THRESHOLD <= 0.)
                throw std::runtime_error("Sum-product message LLR threshold must be > 0!");
        }

        const auto &r_qber_maps = config["code_rate_QBER_maps"];
        for (const auto &m : r_qber_maps)
        {
            cfg.R_QBER_MAPS.push_back({m["code_rate"].template get<double>(), m["QBER_begin"].template get<double>(),
                                       m["QBER_end"].template get<double>(), m["QBER_step"].template get<double>()});
        }

        if (cfg.R_QBER_MAPS.empty())
            throw std::runtime_error("Array with code rate(R) and QBER maps is empty!");
        for (size_t i = 0; i < cfg.R_QBER_MAPS.size(); i++)
        {
            if (cfg.R_QBER_MAPS[i].code_rate <= 0. || cfg.R_QBER_MAPS[i].code_rate >= 1.)
                throw std::runtime_error("Code rate(R) must be: 0 < R < 1!");
            if (cfg.R_QBER_MAPS[i].QBER_begin <= 0. || cfg.R_QBER_MAPS[i].QBER_begin >= 1. || cfg.R_QBER_MAPS[i].QBER_end <= 0. || cfg.R_QBER_MAPS[i].QBER_end >= 1. || cfg.R_QBER_MAPS[i].QBER_begin > cfg.R_QBER_MAPS[i].QBER_end)
                throw std::runtime_error("Invalid QBER begin or end parameters. QBER must be: 0 < QBER < 1, and begin cannot be larger than end!");
            if (cfg.R_QBER_MAPS[i].QBER_step <= 0.)
                throw std::runtime_error("QBER step must be > 0!");
            if (cfg.R_QBER_MAPS[i].QBER_begin != cfg.R_QBER_MAPS[i].QBER_end)
            {
                if (cfg.R_QBER_MAPS[i].QBER_step - EPSILON > cfg.R_QBER_MAPS[i].QBER_end - cfg.R_QBER_MAPS[i].QBER_begin)
                    throw std::runtime_error("QBER step is too large.");
            }
        }
        std::sort(cfg.R_QBER_MAPS.begin(), cfg.R_QBER_MAPS.end(),
            [](R_QBER_map &a, R_QBER_map &b)
            {
                return (a.code_rate < b.code_rate);
            });

        return cfg;
    }
    catch (const std::exception &e)
    {
        fmt::print(stderr, fg(fmt::color::red), "An error occurred while reading a configuration parameter.\n");
        throw;
    }
}
