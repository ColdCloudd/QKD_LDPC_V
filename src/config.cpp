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
        cfg.USE_MIN_SUM_DECODING_ALG = config["use_min_sum_decoding_algorithm"].template get<bool>();

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

        const auto &params = config["code_rate_QBER_parameters"];
        for (const auto &p : params)
        {
            cfg.R_QBER_PARAMETERS.push_back({p["code_rate"].template get<double>(), p["QBER_begin"].template get<double>(),
                                             p["QBER_end"].template get<double>(), p["QBER_step"].template get<double>()});
        }

        if (cfg.R_QBER_PARAMETERS.empty())
        {
            throw std::runtime_error("Array with code rate and QBER parameters is empty!");
        }
        for (size_t i = 0; i < cfg.R_QBER_PARAMETERS.size(); i++)
        {
            if (cfg.R_QBER_PARAMETERS[i].code_rate <= 0. || cfg.R_QBER_PARAMETERS[i].code_rate >= 1.)
            {
                throw std::runtime_error("Code rate(R) must be: 0 < R < 1!");
            }
            if (cfg.R_QBER_PARAMETERS[i].QBER_begin <= 0. || cfg.R_QBER_PARAMETERS[i].QBER_begin >= 1. || cfg.R_QBER_PARAMETERS[i].QBER_end <= 0. || cfg.R_QBER_PARAMETERS[i].QBER_end >= 1. || cfg.R_QBER_PARAMETERS[i].QBER_begin >= cfg.R_QBER_PARAMETERS[i].QBER_end)
            {
                throw std::runtime_error("Invalid QBER begin or end parameters. QBER must be: 0 < QBER < 1, and begin must be less than end.");
            }
            if (cfg.R_QBER_PARAMETERS[i].QBER_step <= 0.)
            {
                throw std::runtime_error("QBER step must be > 0!");
            }
            const double epsilon = 1e-6;
            if (cfg.R_QBER_PARAMETERS[i].QBER_step - epsilon > cfg.R_QBER_PARAMETERS[i].QBER_end - cfg.R_QBER_PARAMETERS[i].QBER_begin)
            {
                throw std::runtime_error("QBER step is too large.");
            }
        }
        std::sort(cfg.R_QBER_PARAMETERS.begin(), cfg.R_QBER_PARAMETERS.end(),
                  [](R_QBER_params &a, R_QBER_params &b)
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
