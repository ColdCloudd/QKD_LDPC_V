#include "array_and_matrix_operations.hpp"

// Convert a dense parity check matrix into an array containing information about bit nodes and associated check nodes (sparse matrix).
std::vector<std::vector<int>> get_bit_nodes(const std::vector<std::vector<int8_t>> &matrix, 
                                            const std::vector<int> &bit_nodes_weight)
{
    size_t num_bit_nodes = matrix[0].size();
    size_t num_check_nodes = matrix.size();

    std::vector<std::vector<int>> bit_nodes(num_bit_nodes);

    for (int i = 0; i < num_bit_nodes; i++)
    {
        bit_nodes[i].reserve(bit_nodes_weight[i]);
        for (int j = 0; j < num_check_nodes; j++)
        {
            if (matrix[j][i] == 1)
            {
                bit_nodes[i].push_back(j);
            }
        }
    }
    return bit_nodes;
}

// Convert a dense parity check matrix into an array containing information about check nodes and associated bit nodes (sparse matrix).
std::vector<std::vector<int>> get_check_nodes(const std::vector<std::vector<int8_t>> &matrix, 
                                              const std::vector<int> &check_nodes_weight)
{
    size_t num_bit_nodes = matrix[0].size();
    size_t num_check_nodes = matrix.size();

    std::vector<std::vector<int>> check_nodes(num_check_nodes);

    for (int i = 0; i < num_check_nodes; i++)
    {
        check_nodes[i].reserve(check_nodes_weight[i]);
        for (int j = 0; j < num_bit_nodes; j++)
        {
            if (matrix[i][j] == 1)
            {
                check_nodes[i].push_back(j);
            }
        }
    }
    return check_nodes;
}

// Finding the maximum modulo LLR value in a given matrix.
double get_max_llr(const std::vector<std::vector<double>> &matrix)
{
    double max_abs_llr = 0;
    double curr_abs_llr = 0;
    for (size_t i = 0; i < matrix.size(); ++i)
    {
        for (size_t j = 0; j < matrix[i].size(); ++j)
        {
            curr_abs_llr = std::abs(matrix[i][j]);
            if (curr_abs_llr > max_abs_llr)
            {
                max_abs_llr = curr_abs_llr;
            }
        }
    }
    return max_abs_llr;
}

bool arrays_equal(const std::vector<int> &array1,
                  const std::vector<int> &array2)
{
    if (array1.size() != array2.size())
    {
        return false;
    }
    for (size_t i = 0; i < array1.size(); i++)
    {
        if (array1[i] != array2[i])
        {
            return false;
        }
    }
    return true;
}

// Find the first element in array1 that is not present in array2
int find_available_index(const std::vector<int>& array1, 
                         const std::vector<int>& array2) 
{
    for (int item : array1) 
    {
        // Check if the current item is not found in array2
        if (std::find(array2.begin(), array2.end(), item) == array2.end()) 
        {
            return item;
        }
    }
    // If all elements in array1 are present in array2, return -1
    return -1;
}

// Determine positions of bits to be removed from key based on parity-check matrix for privacy maintenance.
void get_bits_positions_to_remove(H_matrix &matrix)
{
    auto& bit_nodes = matrix.bit_nodes;

    // Vector of pairs: (index of the bit node, corresponding check node indices)
    std::vector<std::pair<int, std::vector<int>>> indexed_bit_nodes;
    indexed_bit_nodes.reserve(bit_nodes.size());
    for (size_t i = 0; i < bit_nodes.size(); ++i) 
    {
        indexed_bit_nodes.emplace_back(i, bit_nodes[i]);
    }

    // Sort bit nodes by the number of their connections (ascending order). Sort by column weight in H.
    sort(indexed_bit_nodes.begin(), indexed_bit_nodes.end(),
        [](const std::pair<int, std::vector<int>>& a, const std::pair<int, std::vector<int>>& b) 
        {
            return a.second.size() < b.second.size();
        });
    
    size_t num_check_nodes = matrix.check_nodes.size();

    // Indices of bits to be removed from key
    matrix.bits_to_remove.reserve(num_check_nodes);

    // For tracking check nodes that already processed
    std::vector<int> check_node_indexes;
    check_node_indexes.reserve(num_check_nodes);

    for (const auto& bit_node : indexed_bit_nodes) 
    {
        int bit_node_idx = bit_node.first;
        const std::vector<int>& check_idx = bit_node.second;
        int idx = find_available_index(check_idx, check_node_indexes);
        if (idx != -1) 
        {
            matrix.bits_to_remove.push_back(bit_node_idx);
            check_node_indexes.push_back(idx);      // Mark the check node index as used
        }
    }
    // Sort in ascending order
    std::sort(matrix.bits_to_remove.begin(), matrix.bits_to_remove.end());
}

// Maintaining privacy by removing m bits from the Alice and Bob keys.
void privacy_maintenance(const H_matrix &matrix, 
                         const std::vector<int>& array1, 
                         const std::vector<int>& array2,
                         std::vector<int>& array1_out,
                         std::vector<int>& array2_out)
{
    auto& bits_to_remove = matrix.bits_to_remove;

    size_t btr_len = bits_to_remove.size();
    size_t new_arr_len = array1.size() - btr_len;
    array1_out.resize(new_arr_len);
    array2_out.resize(new_arr_len);

    size_t n = 0;
    size_t m = 0;
    for (size_t i = 0; i < array1.size(); ++i)
    {
        if (n < btr_len && bits_to_remove[n] == i)
        {
            ++n;
        }
        else
        {
            array1_out[m] = array1[i];
            array2_out[m] = array2[i];
            ++m;
        }
    }
}

// Function for reading a sparse matrix from a file in alist format (https://rptu.de/channel-codes/matrix-file-formats).
H_matrix read_sparse_alist_matrix(const fs::path &matrix_path)
{
    std::vector<std::string> line_vec;
    std::ifstream file(matrix_path);

    if (!file.is_open())
    {
        throw std::runtime_error("Failed to open file: " + matrix_path.string());
    }
    std::string line;

    while (getline(file, line))
    {
        line_vec.push_back(line);
    }
    file.close();

    if (line_vec.empty())
    {
        throw std::runtime_error("File is empty or cannot be read properly: " + matrix_path.string());
    }

    std::vector<std::vector<int>> vec_int;
    try
    {
        for (const auto &line : line_vec)
        {
            std::istringstream iss(line);
            std::vector<int> numbers;
            int number;
            while (iss >> number)
            {
                numbers.push_back(number);
            }
            vec_int.push_back(numbers);
        }
    }
    catch (const std::exception &e)
    {
        fmt::print(stderr, fg(fmt::color::red), "An error occurred while parsing file: {}\n", matrix_path.string());
        throw;
    }

    if (vec_int.size() < 4)
    {
        throw std::runtime_error("Insufficient data in the file: " + matrix_path.string());
    }

    if (vec_int[0].size() != 2 || vec_int[1].size() != 2)
    {
        throw std::runtime_error("File format does not match the alist format: " + matrix_path.string());
    }

    size_t col_num = vec_int[0][0];     // n
    size_t row_num = vec_int[0][1];     // m

    size_t max_col_weight = vec_int[1][0];      // d^(v)_max
    size_t max_row_weight = vec_int[1][1];      // d^(c)_max

    size_t num_bit_nodes = vec_int[2].size();       // n
    size_t num_check_nodes = vec_int[3].size();     // m

    size_t curr_line = 4;

    if (vec_int.size() < curr_line + num_bit_nodes + num_check_nodes)
    {
        throw std::runtime_error("Insufficient data in the file: " + matrix_path.string());
    }

    if (col_num != num_bit_nodes)
    {
        throw std::runtime_error("Number of columns '" + std::to_string(col_num) + "' is not the same as the length of the third line '" + std::to_string(num_bit_nodes) + "'. File: " + matrix_path.string());
    }
    else if (row_num != num_check_nodes)
    {
        throw std::runtime_error("Number of rows '" + std::to_string(row_num) + "' is not the same as the length of the fourth line '" + std::to_string(num_check_nodes) + "'. File: " + matrix_path.string());
    }

    // Initialization of node weights and check for matrix regularity
    bool is_regular = true;
    std::vector<int> bit_nodes_weight(num_bit_nodes);
    std::vector<int> check_nodes_weight(num_check_nodes);
    for (size_t i = 0; i < num_bit_nodes; i++)
    {
        bit_nodes_weight[i] = vec_int[2][i];
        if (vec_int[2][i] != vec_int[2][0])
        {
            is_regular = false;
        }
    }
    for (size_t i = 0; i < num_check_nodes; i++)
    {
        check_nodes_weight[i] = vec_int[3][i];
        if (vec_int[3][i] != vec_int[3][0])
        {
            is_regular = false;
        }
    }

    // Check if the number of non-zero elements corresponds to the node weights
    size_t non_zero_num;
    for (size_t i = 0; i < num_bit_nodes; i++)
    {
        non_zero_num = 0;
        for (size_t j = 0; j < vec_int[curr_line + i].size(); j++)
        {
            if (vec_int[curr_line + i][j] != 0)
            {
                non_zero_num++;
            }
        }
        if (non_zero_num != vec_int[2][i])
        {
            throw std::runtime_error("Number of non-zero elements '" + std::to_string(non_zero_num) + "' in the line '" + std::to_string(curr_line + i + 1) + "' does not match the weight in the third line '" + std::to_string(vec_int[2][i]) + "'. File: " + matrix_path.string());
        }
    }

    curr_line += num_bit_nodes;
    for (size_t i = 0; i < num_check_nodes; i++)
    {
        non_zero_num = 0;
        for (size_t j = 0; j < vec_int[curr_line + i].size(); j++)
        {
            if (vec_int[curr_line + i][j] != 0)
            {
                non_zero_num++;
            }
        }
        if (non_zero_num != vec_int[3][i])
        {
            throw std::runtime_error("Number of non-zero elements '" + std::to_string(non_zero_num) + "' in the line '" + std::to_string(curr_line + i + 1) + "' does not match the weight in the fourth line '" + std::to_string(vec_int[3][i]) + "'. File: " + matrix_path.string());
        }
    }

    // Filling the matrix of bit nodes
    H_matrix matrix_out{};
    try
    {
        curr_line = 4;
        matrix_out.bit_nodes.resize(num_bit_nodes);
        for (size_t i = 0; i < num_bit_nodes; ++i)
        {
            matrix_out.bit_nodes[i].reserve(bit_nodes_weight[i]);
            for (size_t j = 0; j < bit_nodes_weight[i]; ++j)
            {
                matrix_out.bit_nodes[i].push_back(vec_int[curr_line + i][j] - 1);
            }
        }
    }
    catch (const std::exception &e)
    {
        fmt::print(stderr, fg(fmt::color::red), "An error occurred while creating 'bit_nodes' matrix from file: {}\n", matrix_path.string());
        throw;
    }

    // Filling the matrix of check nodes
    try
    {
        curr_line += num_bit_nodes;
        matrix_out.check_nodes.resize(num_check_nodes);
        for (size_t i = 0; i < num_check_nodes; ++i)
        {
            matrix_out.check_nodes[i].reserve(check_nodes_weight[i]);
            for (size_t j = 0; j < check_nodes_weight[i]; ++j)
            {
                matrix_out.check_nodes[i].push_back(vec_int[curr_line + i][j] - 1);
            }
        }
    }
    catch (const std::exception &e)
    {
        fmt::print(stderr, fg(fmt::color::red), "An error occurred while creating 'check_nodes' matrix from file: {}\n", matrix_path.string());
        throw;
    }

    matrix_out.is_regular = is_regular;
    if (CFG.ENABLE_PRIVACY_MAINTENANCE)
    {
        get_bits_positions_to_remove(matrix_out);
    }
    return matrix_out;
}

// Read dense matrix from file.
H_matrix read_dense_matrix(const fs::path &matrix_path)
{
    std::vector<std::string> line_vec;
    std::ifstream file(matrix_path);

    if (!file.is_open())
    {
        throw std::runtime_error("Failed to open file: " + matrix_path.string());
    }
    std::string line;

    while (getline(file, line))
    {
        line_vec.push_back(line);
    }
    file.close();

    if (line_vec.empty())
    {
        throw std::runtime_error("File is empty or cannot be read properly: " + matrix_path.string());
    }

    std::vector<std::vector<int8_t>> dense_matrix;
    try
    {
        for (const auto &line : line_vec)
        {
            std::istringstream iss(line);
            std::vector<int8_t> numbers;
            int number;
            while (iss >> number)
            {

                if (number != 0 && number != 1)
                {
                    throw std::runtime_error("Parity check matrix can only take values ​​0 or 1.");
                }
                numbers.push_back(static_cast<int8_t>(number));
            }
            dense_matrix.push_back(numbers);
        }
    }
    catch (const std::exception &e)
    {
        fmt::print(stderr, fg(fmt::color::red), "An error occurred while parsing file: {}\n", matrix_path.string());
        throw;
    }

    for (size_t i = 0; i < dense_matrix.size(); i++)
    {
        if (dense_matrix[0].size() != dense_matrix[i].size())
        {
            throw std::runtime_error("Different lengths of rows in a matrix. File: " + matrix_path.string());
        }
    }

    size_t col_num = dense_matrix[0].size();    //n - bit_nodes number
    size_t row_num = dense_matrix.size();       //m - check nodes number

    std::vector<int> bit_nodes_weight(col_num);
    std::vector<int> check_nodes_weight(row_num);

    // Counting column weights and determining the maximum weight
    size_t curr_weight = 0;
    size_t max_col_weight = 0;
    for (size_t i = 0; i < col_num; i++)
    {
        curr_weight = 0;
        for (size_t j = 0; j < row_num; j++)
        {
            curr_weight += dense_matrix[j][i];
        }
        if (curr_weight <= 0)
        {
            throw std::runtime_error("Column '" + std::to_string(i + 1) + "' weight cannot be equal to or less than zero. File: " + matrix_path.string());
        }
        bit_nodes_weight[i] = static_cast<int>(curr_weight);
        if (curr_weight > max_col_weight)
        {
            max_col_weight = curr_weight;
        }
    }

    // Counting row weights and determining the maximum weight
    size_t max_row_weight = 0;
    for (size_t i = 0; i < row_num; i++)
    {
        curr_weight = accumulate(dense_matrix[i].begin(), dense_matrix[i].end(), 0);
        if (curr_weight <= 0)
        {
            throw std::runtime_error("Row '" + std::to_string(i + 1) + "' weight cannot be equal to or less than zero. File: " + matrix_path.string());
        }
        check_nodes_weight[i] = static_cast<int>(curr_weight);
        if (curr_weight > max_row_weight)
        {
            max_row_weight = curr_weight;
        }
    }

    bool is_regular = true;
    for (size_t i = 0; i < col_num; i++)
    {
        if (bit_nodes_weight[0] != bit_nodes_weight[i])
        {
            is_regular = false;
        }
    }

    for (size_t i = 0; i < row_num; i++)
    {
        if (check_nodes_weight[0] != check_nodes_weight[i])
        {
            is_regular = false;
        }
    }

    // Filling bit and check nodes 
    H_matrix matrix_out{};
    matrix_out.bit_nodes = get_bit_nodes(dense_matrix, bit_nodes_weight);
    matrix_out.check_nodes = get_check_nodes(dense_matrix, check_nodes_weight);
    if (CFG.ENABLE_PRIVACY_MAINTENANCE)
    {
        get_bits_positions_to_remove(matrix_out);
    }
    matrix_out.is_regular = is_regular;
    return matrix_out;
}

// Generates Alice's key.
void fill_random_bits(XoshiroCpp::Xoshiro256PlusPlus &prng,
                      std::vector<int> &bit_array)
{
    std::uniform_int_distribution<int> distribution(0, 1);

    // Generate random bits and fill the vector
    for (size_t i = 0; i < bit_array.size(); ++i)
    {
        bit_array[i] = distribution(prng);
    }
}

// Generates Bob's key by making errors in Alice's key. Generates the exact number of errors in the key and returns the exact QBER.
double introduce_errors(XoshiroCpp::Xoshiro256PlusPlus &prng, 
                        const std::vector<int> &bit_array, 
                        double QBER, 
                        std::vector<int> &bit_array_with_errors_out)
{
    size_t array_length = bit_array.size();
    size_t num_errors = static_cast<size_t>(array_length * QBER);

    bit_array_with_errors_out = bit_array;

    if (num_errors > 0)
    {
        std::vector<size_t> error_positions(array_length);
        for (size_t i = 0; i < array_length; ++i)
        {
            error_positions[i] = i;
        }

        std::shuffle(error_positions.begin(), error_positions.end(), prng);

        for (size_t i = 0; i < num_errors; ++i)
        {
            bit_array_with_errors_out[error_positions[i]] ^= 1;
        }
    }
    return static_cast<double>(num_errors) / array_length;
}

// Computes the key syndrome using parity check matrix.
void calculate_syndrome(const std::vector<int> &bit_array,
                        const H_matrix &matrix, 
                        std::vector<int> &syndrome_out)
{
    std::fill(syndrome_out.begin(), syndrome_out.end(), 0);
    for (size_t i = 0; i < matrix.check_nodes.size(); ++i)          // CN number
    {
        for (size_t j = 0; j < matrix.check_nodes[i].size(); ++j)   // Their weights
        {
            syndrome_out[i] ^= bit_array[matrix.check_nodes[i][j]];
        }
    }
}

// Limiting the LLR values of messages in irregular matrix to a given threshold.
void threshold_matrix(std::vector<std::vector<double>> &matrix, 
                      const double &msg_threshold)
{
    for (size_t i = 0; i < matrix.size(); ++i)
    {
        for (size_t j = 0; j < matrix[i].size(); ++j)
        {
            if (matrix[i][j] > msg_threshold)
            {
                matrix[i][j] = msg_threshold;
            }
            else if (matrix[i][j] < -msg_threshold)
            {
                matrix[i][j] = -msg_threshold;
            }
        }
    }
}
