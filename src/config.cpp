#include "config.hpp"

// Reads user-defined configuration parameters from a .json file.
config_data get_config_data(fs::path config_path)
{
    if (!fs::exists(config_path))
    {
        throw std::runtime_error("Configuration file not found: " + config_path.string());
    }

    std::ifstream config_file(config_path);
    if (!config_file.is_open())
    {
        throw std::runtime_error("Failed to open configuration file: " + config_path.string());
    }

    json config = json::parse(config_file);
    config_file.close();
    if (config.empty())
    {
        throw std::runtime_error("Configuration file is empty: " + config_path.string());
    }

    try
    {
        config_data cfg{};
        cfg.THREADS_NUMBER = config["threads_number"].template get<size_t>();
        if (cfg.THREADS_NUMBER < 1)
        {
            throw std::runtime_error("Number of threads must be >= 1!");
        }

        cfg.TRIALS_NUMBER = config["trials_number"].template get<size_t>();
        if (cfg.TRIALS_NUMBER < 1)
        {
            throw std::runtime_error("Number of trials must be >= 1!");
        }

        if (config["use_config_simulation_seed"].template get<bool>())
        {
            cfg.SIMULATION_SEED = config["simulation_seed"].template get<size_t>();
        }
        else
        {
            cfg.SIMULATION_SEED = time(nullptr);
        }

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
        if (cfg.DECODING_ALGORITHM > 3)
        {
            throw std::runtime_error("Only four options are available: \n0 - SPA;\n1 - SPA (with linear approximation of tanh and atanh);\n2 - NMSA;\n3 - OMSA.");
        }
        if (cfg.DECODING_ALGORITHM == 2 || cfg.DECODING_ALGORITHM == 3)
        {   
            nlohmann::json alg_params;
            if (cfg.DECODING_ALGORITHM == 2)
            {
                alg_params = config["min_sum_normalized_parameters"];
                cfg.USE_SCALING_FACTOR_RANGE = alg_params["use_alpha_range"].template get<bool>();
            }
            else
            {
                alg_params = config["min_sum_offset_parameters"];
                cfg.USE_SCALING_FACTOR_RANGE = alg_params["use_beta_range"].template get<bool>();
            }
            
            if (cfg.USE_SCALING_FACTOR_RANGE)
            {
                nlohmann::json scaling_factor_range;
                if (cfg.DECODING_ALGORITHM == 2)
                    scaling_factor_range = alg_params["alpha_range"];
                else
                    scaling_factor_range = alg_params["beta_range"];
                
                cfg.SCALING_FACTOR_RANGE = {scaling_factor_range["begin"].template get<double>(), scaling_factor_range["end"].template get<double>(), scaling_factor_range["step"].template get<double>()};
                if (cfg.SCALING_FACTOR_RANGE.begin <= 0. || cfg.SCALING_FACTOR_RANGE.end <= 0. || cfg.SCALING_FACTOR_RANGE.step <= 0.)
                {
                    throw std::runtime_error("Alpha (or beta) range begin, end, step must be > 0!");
                }
                if (cfg.SCALING_FACTOR_RANGE.begin > cfg.SCALING_FACTOR_RANGE.end)
                {
                    throw std::runtime_error("Alpha (or beta) range begin cannot be larger than end!");
                }
                if (cfg.SCALING_FACTOR_RANGE.begin != cfg.SCALING_FACTOR_RANGE.end)
                {
                    if (cfg.SCALING_FACTOR_RANGE.step - EPSILON > cfg.SCALING_FACTOR_RANGE.end - cfg.SCALING_FACTOR_RANGE.begin)
                    {
                        throw std::runtime_error("Alpha (or beta) range step is too large!");
                    }
                }
            }
            else
            {
                nlohmann::json r_scaling_factor_maps;
                if (cfg.DECODING_ALGORITHM == 2)
                {
                    r_scaling_factor_maps = alg_params["code_rate_alpha_maps"];
                    for (const auto &m : r_scaling_factor_maps)
                    {
                        cfg.R_SCALING_FACTOR_MAPS.push_back({m["code_rate"].template get<double>(), m["alpha"].template get<double>()});
                    }
                }
                else
                {
                    r_scaling_factor_maps = alg_params["code_rate_beta_maps"];
                    for (const auto &m : r_scaling_factor_maps)
                    {
                        cfg.R_SCALING_FACTOR_MAPS.push_back({m["code_rate"].template get<double>(), m["beta"].template get<double>()});
                    }
                }
                
                if (cfg.R_SCALING_FACTOR_MAPS.empty())
                {
                    throw std::runtime_error("Array with code rate(R) and alpha (or beta) maps is empty!");
                }
                for (size_t i = 0; i < cfg.R_SCALING_FACTOR_MAPS.size(); i++)
                {
                    if (cfg.R_SCALING_FACTOR_MAPS[i].code_rate <= 0. || cfg.R_SCALING_FACTOR_MAPS[i].code_rate >= 1.)
                    {
                        throw std::runtime_error("Code rate(R) must be: 0 < R < 1!");
                    }
                    if (cfg.R_SCALING_FACTOR_MAPS[i].scaling_factor <= 0.)
                    {
                        throw std::runtime_error("Alpha (or beta) must be > 0!");
                    }
                }
                std::sort(cfg.R_SCALING_FACTOR_MAPS.begin(), cfg.R_SCALING_FACTOR_MAPS.end(),
                  [](R_scaling_factor_map &a, R_scaling_factor_map &b)
                  {
                        return (a.code_rate < b.code_rate);
                  });
            }
        }
        
        cfg.DECODING_ALG_MAX_ITERATIONS = config["decoding_algorithm_max_iterations"].template get<size_t>();
        if (cfg.DECODING_ALG_MAX_ITERATIONS < 1)
        {
            throw std::runtime_error("Minimum number of decoding algorithm iterations must be >= 1!");
        }

        cfg.MATRIX_FORMAT = config["matrix_format"].template get<size_t>();
        if (cfg.MATRIX_FORMAT > 3)
        {
            throw std::runtime_error("Only four options are available: \n0 - dense;\n1 - sparse alist;\n2 - sparse_1;\n3 - sparse_2.");
        }

        cfg.TRACE_QKD_LDPC = config["trace_qkd_ldpc"].template get<bool>();
        cfg.TRACE_DECODING_ALG = config["trace_decoding_algorithm"].template get<bool>();
        cfg.TRACE_DECODING_ALG_LLR = config["trace_decoding_algorithm_llr"].template get<bool>();
        cfg.ENABLE_DECODING_ALG_MSG_LLR_THRESHOLD = config["enable_decoding_algorithm_msg_llr_threshold"].template get<bool>();

        if (cfg.ENABLE_DECODING_ALG_MSG_LLR_THRESHOLD)
        {
            cfg.DECODING_ALG_MSG_LLR_THRESHOLD = config["decoding_algorithm_msg_llr_threshold"].template get<double>();
            if (cfg.DECODING_ALG_MSG_LLR_THRESHOLD <= 0.)
            {
                throw std::runtime_error("Sum-product message LLR threshold must be > 0!");
            }
        }

        const auto &r_qber_maps = config["code_rate_QBER_maps"];
        for (const auto &m : r_qber_maps)
        {
            cfg.R_QBER_MAPS.push_back({m["code_rate"].template get<double>(), m["QBER_begin"].template get<double>(),
                                       m["QBER_end"].template get<double>(), m["QBER_step"].template get<double>()});
        }

        if (cfg.R_QBER_MAPS.empty())
        {
            throw std::runtime_error("Array with code rate(R) and QBER maps is empty!");
        }
        for (size_t i = 0; i < cfg.R_QBER_MAPS.size(); i++)
        {
            if (cfg.R_QBER_MAPS[i].code_rate <= 0. || cfg.R_QBER_MAPS[i].code_rate >= 1.)
            {
                throw std::runtime_error("Code rate(R) must be: 0 < R < 1!");
            }
            if (cfg.R_QBER_MAPS[i].QBER_begin <= 0. || cfg.R_QBER_MAPS[i].QBER_begin >= 1. || cfg.R_QBER_MAPS[i].QBER_end <= 0. || cfg.R_QBER_MAPS[i].QBER_end >= 1. || cfg.R_QBER_MAPS[i].QBER_begin > cfg.R_QBER_MAPS[i].QBER_end)
            {
                throw std::runtime_error("Invalid QBER begin or end parameters. QBER must be: 0 < QBER < 1, and begin cannot be larger than end!");
            }
            if (cfg.R_QBER_MAPS[i].QBER_step <= 0.)
            {
                throw std::runtime_error("QBER step must be > 0!");
            }
            if (cfg.R_QBER_MAPS[i].QBER_begin != cfg.R_QBER_MAPS[i].QBER_end)
            {
                if (cfg.R_QBER_MAPS[i].QBER_step - EPSILON > cfg.R_QBER_MAPS[i].QBER_end - cfg.R_QBER_MAPS[i].QBER_begin)
                {
                    throw std::runtime_error("QBER step is too large.");
                }
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
