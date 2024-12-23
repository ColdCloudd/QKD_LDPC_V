#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <filesystem>

#include <fmt/core.h>
#include <fmt/color.h>
#include <fmt/ranges.h>
#include <fmt/format.h>

namespace fs = std::filesystem;

void print_array(const std::vector<int> &array);

void print_array(const std::vector<double> &array);

std::vector<fs::path> get_file_paths_in_directory(const fs::path &directory_path);

fs::path select_matrix_file(const std::vector<fs::path> &matrix_paths);

//Outputs to the console matrix.
template <typename T>
void print_matrix(const std::vector<std::vector<T>> &matrix)
{
    for (size_t i = 0; i < matrix.size(); i++)
    {
        for (size_t j = 0; j < matrix[i].size(); j++)
        {
            fmt::print(fg(fmt::color::blue), "{:.4} ", matrix[i][j]);
        }
        fmt::print("\n");
    }
}
