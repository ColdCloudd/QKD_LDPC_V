#include "utils.hpp"

void print_array(const std::vector<int> &array)
{
    for (size_t i = 0; i < array.size(); i++)
    {
        fmt::print(fg(fmt::color::blue), "{} ", array[i]);
    }
}

void print_array(const std::vector<double> &array)
{
    for (size_t i = 0; i < array.size(); i++)
    {
        fmt::print(fg(fmt::color::blue), "{:.4} ", array[i]);
    }
}

// Gets paths to files with the specified extension in the specified directory.
std::vector<fs::path> get_file_paths_in_directory(
    const fs::path &directory_path,
     const std::string &extension
)
{
    std::vector<fs::path> file_paths;
    try
    {
        if (!fs::exists(directory_path))
            throw std::runtime_error(fmt::format("Directory '{}' doesn't exist.", directory_path.string()));
        if (!fs::is_directory(directory_path))
            throw std::runtime_error(fmt::format("Path '{}' is not a directory.", directory_path.string()));

        std::string ext = extension;
        if (!ext.empty() && ext[0] != '.')
            ext = "." + ext;

        for (const auto &entry : fs::directory_iterator(directory_path))
        {
            if (fs::is_regular_file(entry.path()) && entry.path().extension().string() == ext)
                file_paths.push_back(entry.path());
        }

        if (file_paths.empty())
            throw std::runtime_error(fmt::format("No files with extension '{}' in the directory '{}'.", extension, directory_path.string()));
    }
    catch (const std::exception &e)
    {
        fmt::print(stderr, fg(fmt::color::red), "An error occurred while getting file paths in directory: {}\n", directory_path.string());
        throw;
    }

    return file_paths;
}
