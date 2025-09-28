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

    return {.begin = begin, .end = end, .step = step};
}

std::vector<R_scaling_factor_map> parse_scaling_factor_maps(
    const json& scaling_factor_maps,
    const std::string& key
) 
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

void print_config_info(
    config_data cfg, 
    std::string cfg_name,
    size_t cfg_number
)
{
    fmt::print(fg(fmt::color::yellow), "------------------------- CONFIG #{} INFO --------------------------\n", cfg_number);
    fmt::print(fg(fmt::color::yellow), "Config name: {}\n", fmt::styled(cfg_name, fg(fmt::color::crimson)));
    fmt::print(fg(fmt::color::yellow), "Threads number: {}\n", fmt::styled(cfg.THREADS_NUMBER, fg(fmt::color::crimson)));
    fmt::print(fg(fmt::color::yellow), "Trials number: {}\n", fmt::styled(cfg.TRIALS_NUMBER, fg(fmt::color::crimson)));
    fmt::print(fg(fmt::color::yellow), "Simulation seed: {}\n", fmt::styled(cfg.SIMULATION_SEED, fg(fmt::color::crimson)));
    fmt::print(fg(fmt::color::yellow), "Privacy maintenance: {}\n", fmt::styled((cfg.ENABLE_PRIVACY_MAINTENANCE ? "Enabled" : "Disabled"), fg(fmt::color::crimson)));
    fmt::print(fg(fmt::color::yellow), "Throughput measurement: {}\n", fmt::styled((cfg.ENABLE_THROUGHPUT_MEASUREMENT ? ("Enabled, RTT = " + fmt::format("{:.3f}", CFG.RTT) + " ms") : "Disabled"), fg(fmt::color::crimson)));

    std::string alg_name =
        (cfg.DECODING_ALGORITHM == DEC_SPA) ? "SPA" :
        (cfg.DECODING_ALGORITHM == DEC_SPA_APPROX) ? "SPA(lin approx)" :
        (cfg.DECODING_ALGORITHM == DEC_NMSA) ? "NMSA" :
        (cfg.DECODING_ALGORITHM == DEC_OMSA) ? "OMSA" :
        (cfg.DECODING_ALGORITHM == DEC_ANMSA) ? "ANMSA" :
        (cfg.DECODING_ALGORITHM == DEC_AOMSA) ? "AOMSA" : "Unknown";
    fmt::print(fg(fmt::color::yellow), "Decoding algorithm: {}\n", fmt::styled(alg_name, fg(fmt::color::crimson)));
    fmt::print(fg(fmt::color::yellow), "Decoding algorithm maximum iterations: {}\n", fmt::styled(cfg.DECODING_ALG_MAX_ITERATIONS, fg(fmt::color::crimson)));

    std::string mat_format =
        (cfg.MATRIX_FORMAT == MAT_SPARSE_UNCOMPRESSED) ? "Sparse (uncompressed)" :
        (cfg.MATRIX_FORMAT == MAT_SPARSE_ALIST) ? "Sparse (alist)" :
        (cfg.MATRIX_FORMAT == MAT_SPARSE_1) ? "Sparse (1)" :
        (cfg.MATRIX_FORMAT == MAT_SPARSE_2) ? "Sparse (2)" : "Unknown";
    fmt::print(fg(fmt::color::yellow), "Parity-check matrix format: {}\n", fmt::styled(mat_format, fg(fmt::color::crimson)));
    std::string rate_adapt_params_option = (cfg.USE_ADAPTATION_PARAMETERS_RANGES) ? " (ranges)" : " (maps)";
    fmt::print(fg(fmt::color::yellow), "Code rate adaptation: {}\n", fmt::styled((cfg.ENABLE_CODE_RATE_ADAPTATION ? ("Enabled" + rate_adapt_params_option): "Disabled"), fg(fmt::color::crimson)));
    fmt::print(fg(fmt::color::yellow), "Untainted puncturing: {}\n", fmt::styled((cfg.ENABLE_UNTAINTED_PUNCTURING ? "Enabled" : "Disabled"), fg(fmt::color::crimson)));
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

        cfg.ENABLE_PRIVACY_MAINTENANCE = config["enable_privacy_maintenance"].template get<bool>();
        cfg.ENABLE_THROUGHPUT_MEASUREMENT = config["enable_throughput_measurement"].template get<bool>();
        if (cfg.ENABLE_THROUGHPUT_MEASUREMENT)
        {
            fmt::print(fg(fmt::color::purple), "WARNING: Throughput measurement is enabled. It is recommended to perform experiments in single-threaded mode.\n");
            
            const auto &tm_params = config["throughput_measurement_parameters"];
            cfg.CONSIDER_RTT = tm_params["consider_RTT"].template get<bool>();
            if (cfg.CONSIDER_RTT)
            {
                cfg.RTT = tm_params["RTT"].template get<double>();
                if (cfg.RTT < 0.)
                    throw std::runtime_error("Round-Trip Time (RTT) must be >= 0!");
            }
        }

        cfg.DECODING_ALGORITHM = config["decoding_algorithm"].template get<size_t>();
        if (cfg.DECODING_ALGORITHM > DEC_AOMSA)
            throw std::runtime_error("Only six options are available: \n0 - SPA;\n1 - SPA (with linear approximation of tanh and atanh);\n2 - NMSA;\n3 - OMSA;\n4 - ANMSA;\n5 - AOMSA.");

        if (cfg.DECODING_ALGORITHM > DEC_SPA_APPROX)
        {   
            nlohmann::json alg_params;
            if (cfg.DECODING_ALGORITHM == DEC_NMSA)
            {
                alg_params = config["min_sum_normalized_parameters"];
                cfg.DECODING_ALG_PARAMS.primary.use_range = alg_params["use_alpha_range"];
                if (cfg.DECODING_ALG_PARAMS.primary.use_range)
                    cfg.DECODING_ALG_PARAMS.primary.range = parse_scaling_factor_range(alg_params["alpha_range"]);
                else
                    cfg.DECODING_ALG_PARAMS.primary.maps = parse_scaling_factor_maps(alg_params["code_rate_alpha_maps"], "alpha");
            }
            else if (cfg.DECODING_ALGORITHM == DEC_OMSA)
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
                if (cfg.DECODING_ALGORITHM == DEC_ANMSA)
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
                else if (cfg.DECODING_ALGORITHM == DEC_AOMSA)
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
                    if (cfg.DECODING_ALGORITHM == DEC_ANMSA)
                    {
                        alg_name = "ANMSA";
                        prim_scal_factor = "alpha";
                        second_scal_factor = "nu";
                    }
                    else if (cfg.DECODING_ALGORITHM == DEC_AOMSA)
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
        if (cfg.MATRIX_FORMAT > MAT_SPARSE_2)
            throw std::runtime_error("Only four options are available: \n0 - uncompressed;\n1 - sparse alist;\n2 - sparse_1;\n3 - sparse_2.");

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

        const auto &r_qber_ranges = config["code_rate_QBER_ranges"];
        for (const auto &r : r_qber_ranges)
        {
            const auto &qber_vals = r["QBER"];
            R_QBER_range r_qber_rng{};
            r_qber_rng.code_rate = r["code_rate"].template get<double>();
            r_qber_rng.QBER_begin = qber_vals["begin"].template get<double>();
            r_qber_rng.QBER_end = qber_vals["end"].template get<double>();
            r_qber_rng.QBER_step = qber_vals["step"].template get<double>();

            cfg.R_QBER_RANGES.push_back(r_qber_rng);
        }

        if (cfg.R_QBER_RANGES.empty())
            throw std::runtime_error("Array with code rate(R) and QBER ranges is empty!");
        for (size_t i = 0; i < cfg.R_QBER_RANGES.size(); i++)
        {
            if (cfg.R_QBER_RANGES[i].code_rate <= 0. || cfg.R_QBER_RANGES[i].code_rate >= 1.)
                throw std::runtime_error("Code rate(R) must be: 0 < R < 1!");
            if (cfg.R_QBER_RANGES[i].QBER_begin <= 0. || cfg.R_QBER_RANGES[i].QBER_begin >= 1. || cfg.R_QBER_RANGES[i].QBER_end <= 0. 
                || cfg.R_QBER_RANGES[i].QBER_end >= 1. || cfg.R_QBER_RANGES[i].QBER_begin > cfg.R_QBER_RANGES[i].QBER_end)
                throw std::runtime_error("Invalid QBER begin or end parameters. QBER must be: 0 < QBER < 1, and begin cannot be larger than end!");
            if (cfg.R_QBER_RANGES[i].QBER_step <= 0.)
                throw std::runtime_error("QBER step must be > 0!");
            if (cfg.R_QBER_RANGES[i].QBER_begin != cfg.R_QBER_RANGES[i].QBER_end)
            {
                if (cfg.R_QBER_RANGES[i].QBER_step - EPSILON > cfg.R_QBER_RANGES[i].QBER_end - cfg.R_QBER_RANGES[i].QBER_begin)
                    throw std::runtime_error("QBER step is too large.");
            }
        }
        std::sort(cfg.R_QBER_RANGES.begin(), cfg.R_QBER_RANGES.end(),
            [](R_QBER_range &a, R_QBER_range &b)
            {
                return (a.code_rate < b.code_rate);
            });

        cfg.ENABLE_CODE_RATE_ADAPTATION = config["enable_code_rate_adaptation"].template get<bool>();
        if(cfg.ENABLE_CODE_RATE_ADAPTATION)
        {
            const auto &r_adapt_params = config["code_rate_adaptation_parameters"];

            cfg.ENABLE_UNTAINTED_PUNCTURING = r_adapt_params["enable_untainted_puncturing"].template get<bool>();

            cfg.USE_ADAPTATION_PARAMETERS_RANGES = r_adapt_params["use_adaptation_parameters_ranges"].template get<bool>();
            if (cfg.USE_ADAPTATION_PARAMETERS_RANGES)
            {
                const auto &r_adapt_params_ranges = r_adapt_params["code_rate_adaptation_parameters_ranges"];
                for (const auto &r : r_adapt_params_ranges)
                {
                    const auto &delta_vals = r["delta"];
                    const auto &efficiency_vals = r["efficiency"];

                    R_adaptation_parameters_range r_adapt_par_rng{};
                    r_adapt_par_rng.code_rate = r["code_rate"].template get<double>();

                    r_adapt_par_rng.delta_begin = delta_vals["begin"].template get<double>();
                    r_adapt_par_rng.delta_end = delta_vals["end"].template get<double>();
                    r_adapt_par_rng.delta_step = delta_vals["step"].template get<double>();

                    r_adapt_par_rng.efficiency_begin = efficiency_vals["begin"].template get<double>();
                    r_adapt_par_rng.efficiency_end = efficiency_vals["end"].template get<double>();
                    r_adapt_par_rng.efficiency_step = efficiency_vals["step"].template get<double>();

                    cfg.R_ADAPT_PARAMS_RANGES.push_back(r_adapt_par_rng);
                }

                if (cfg.R_ADAPT_PARAMS_RANGES.empty())
                    throw std::runtime_error("Array with code rate(R) and adaptation parameters ranges is empty!");
                for (size_t i = 0; i < cfg.R_ADAPT_PARAMS_RANGES.size(); i++)
                {
                    if (cfg.R_ADAPT_PARAMS_RANGES[i].code_rate <= 0. || cfg.R_ADAPT_PARAMS_RANGES[i].code_rate >= 1.)
                        throw std::runtime_error("Code rate(R) must be: 0 < R < 1!");

                    if (cfg.R_ADAPT_PARAMS_RANGES[i].delta_begin <= 0. || cfg.R_ADAPT_PARAMS_RANGES[i].delta_begin >= 1. 
                        || cfg.R_ADAPT_PARAMS_RANGES[i].delta_end <= 0. || cfg.R_ADAPT_PARAMS_RANGES[i].delta_end >= 1. 
                        || cfg.R_ADAPT_PARAMS_RANGES[i].delta_begin > cfg.R_ADAPT_PARAMS_RANGES[i].delta_end)
                        throw std::runtime_error("Invalid delta begin or end parameters. Delta must be: 0 < delta < 1, and begin cannot be larger than end!");
                    if (cfg.R_ADAPT_PARAMS_RANGES[i].delta_step <= 0.)
                        throw std::runtime_error("Delta step must be > 0!");
                    if (cfg.R_ADAPT_PARAMS_RANGES[i].delta_begin != cfg.R_ADAPT_PARAMS_RANGES[i].delta_end)
                    {
                        if (cfg.R_ADAPT_PARAMS_RANGES[i].delta_step - EPSILON > cfg.R_ADAPT_PARAMS_RANGES[i].delta_end - cfg.R_ADAPT_PARAMS_RANGES[i].delta_begin)
                            throw std::runtime_error("Delta step is too large.");
                    }

                    if (cfg.R_ADAPT_PARAMS_RANGES[i].efficiency_begin < 1. || cfg.R_ADAPT_PARAMS_RANGES[i].efficiency_end < 1. 
                        || cfg.R_ADAPT_PARAMS_RANGES[i].efficiency_begin > cfg.R_ADAPT_PARAMS_RANGES[i].efficiency_end)
                        throw std::runtime_error("Invalid efficiency begin or end parameters. Efficiency(f_EC) must be: f_EC >= 1, and begin cannot be larger than end!");
                    if (cfg.R_ADAPT_PARAMS_RANGES[i].efficiency_step <= 0.)
                        throw std::runtime_error("Efficiency step must be > 0!");
                    if (cfg.R_ADAPT_PARAMS_RANGES[i].efficiency_begin != cfg.R_ADAPT_PARAMS_RANGES[i].efficiency_end)
                    {
                        if (cfg.R_ADAPT_PARAMS_RANGES[i].efficiency_step - EPSILON > cfg.R_ADAPT_PARAMS_RANGES[i].efficiency_end - cfg.R_ADAPT_PARAMS_RANGES[i].efficiency_begin)
                            throw std::runtime_error("Efficiency step is too large.");
                    }
                }
                std::sort(cfg.R_ADAPT_PARAMS_RANGES.begin(), cfg.R_ADAPT_PARAMS_RANGES.end(),
                    [](R_adaptation_parameters_range &a, R_adaptation_parameters_range &b)
                    {
                        return (a.code_rate < b.code_rate);
                    });
            }
            else
            {
                const auto &r_qber_adapt_params_maps = r_adapt_params["code_rate_QBER_adaptation_parameters_maps"];
                for (const auto &m : r_qber_adapt_params_maps)
                {
                    R_QBER_adaptation_parameters_map r_qber_adapt_par_map{};
                    
                    r_qber_adapt_par_map.code_rate = m["code_rate"].template get<double>();
                    r_qber_adapt_par_map.QBER_adapt_params.QBER = m["QBER"].template get<double>();
                    r_qber_adapt_par_map.QBER_adapt_params.delta = m["delta"].template get<double>();
                    r_qber_adapt_par_map.QBER_adapt_params.efficiency= m["efficiency"].template get<double>();

                    cfg.R_QBER_ADAPT_PARAMS_MAPS.push_back(r_qber_adapt_par_map);
                }

                if (cfg.R_QBER_ADAPT_PARAMS_MAPS.empty())
                    throw std::runtime_error("Array with code rate(R), QBER and adaptation parameters maps is empty!");
                for (size_t i = 0; i < cfg.R_QBER_ADAPT_PARAMS_MAPS.size(); i++)
                {
                    if (cfg.R_QBER_ADAPT_PARAMS_MAPS[i].code_rate <= 0. || cfg.R_QBER_ADAPT_PARAMS_MAPS[i].code_rate >= 1.)
                        throw std::runtime_error("Code rate(R) must be: 0 < R < 1!");
                    if (cfg.R_QBER_ADAPT_PARAMS_MAPS[i].QBER_adapt_params.QBER <= 0. || cfg.R_QBER_ADAPT_PARAMS_MAPS[i].QBER_adapt_params.QBER >= 1.)
                        throw std::runtime_error("Invalid QBER parameter. QBER must be: 0 < QBER < 1!");
                    if (cfg.R_QBER_ADAPT_PARAMS_MAPS[i].QBER_adapt_params.delta <= 0. || cfg.R_QBER_ADAPT_PARAMS_MAPS[i].QBER_adapt_params.delta >= 1.)
                        throw std::runtime_error("Invalid delta parameter. Delta must be: 0 < delta < 1!");
                    if (cfg.R_QBER_ADAPT_PARAMS_MAPS[i].QBER_adapt_params.efficiency < 1.)
                        throw std::runtime_error("Invalid efficiency parameter. Efficiency(f_EC) must be: f_EC >= 1!");
                }
                std::sort(cfg.R_QBER_ADAPT_PARAMS_MAPS.begin(), cfg.R_QBER_ADAPT_PARAMS_MAPS.end(),
                    [](R_QBER_adaptation_parameters_map &a, R_QBER_adaptation_parameters_map &b)
                    {
                        return (a.code_rate < b.code_rate);
                    });
            }
        }
        return cfg;
    }
    catch (const std::exception &e)
    {
        fmt::print(stderr, fg(fmt::color::red), "An error occurred while reading a configuration parameter.\n");
        throw;
    }
}
