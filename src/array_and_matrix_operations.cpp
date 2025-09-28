#include "array_and_matrix_operations.hpp"

// Convert a uncompressed parity check matrix into an array containing information 
// about bit nodes and associated check nodes (sparse matrix).
std::vector<std::vector<int>> get_bit_nodes(
    const std::vector<std::vector<int8_t>> &matrix, 
    const std::vector<int> &bit_nodes_weight
)
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

// Convert a uncompressed parity check matrix into an array containing information 
// about check nodes and associated bit nodes (sparse matrix).
std::vector<std::vector<int>> get_check_nodes(
    const std::vector<std::vector<int8_t>> &matrix, 
    const std::vector<int> &check_nodes_weight
)
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

std::vector<std::vector<int>> get_bit_nodes_from_check_nodes(const std::vector<std::vector<int>> &check_nodes)
{
    size_t num_bit_nodes = 0;     // bit nodes number
    size_t row_max = 0;
    for (const auto& row : check_nodes) 
    {
        row_max = *std::max_element(row.begin(), row.end());
        if (row_max > num_bit_nodes) 
        {
            num_bit_nodes = row_max;
        }
    }
    ++num_bit_nodes;

    std::vector<std::vector<int>> bit_nodes(num_bit_nodes);
    for (size_t i = 0; i < num_bit_nodes; ++i)
    {
        for (size_t j = 0; j < check_nodes.size(); ++j)
        {
            for (size_t k = 0; k < check_nodes[j].size(); ++k)
            {
                if (check_nodes[j][k] == i)
                {
                    bit_nodes[i].push_back(static_cast<int>(j));
                }
            }
        }
    }
    return bit_nodes;
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

bool arrays_equal(
    const std::vector<int> &array1,
    const std::vector<int> &array2
)
{
    for (size_t i = 0; i < array1.size(); ++i)
    {
        if (array1[i] != array2[i])
        {
            return false;
        }
    }
    return true;
}

// Find the first element in array1 that is not present in array2
int find_available_index(
    const std::vector<int> &array1, 
    const std::vector<int> &array2
) 
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

// Determine positions of bits to be removed from key based on
// parity-check matrix for privacy maintenance.
std::vector<int> get_bits_positions_to_remove(const H_matrix &matrix)
{
    auto& bit_nodes = matrix.bit_nodes;

    // Vector of pairs: (index of the bit node, corresponding check node indices)
    std::vector<std::pair<int, std::vector<int>>> indexed_bit_nodes;
    indexed_bit_nodes.reserve(bit_nodes.size());
    for (size_t i = 0; i < bit_nodes.size(); ++i) 
    {
        indexed_bit_nodes.emplace_back(i, bit_nodes[i]);
    }

    // Sort bit nodes by the number of their connections (ascending order). 
    // Sort by column weight in H.
    sort(indexed_bit_nodes.begin(), indexed_bit_nodes.end(),
        [](const std::pair<int, std::vector<int>>& a, const std::pair<int, std::vector<int>>& b) 
        {
            return a.second.size() < b.second.size();
        });
    
    size_t num_check_nodes = matrix.check_nodes.size();

    // Indices of bits to be removed from key
    std::vector<int> bits_to_remove;
    bits_to_remove.reserve(num_check_nodes);

    // For tracking check nodes that already processed
    std::vector<int> marked_check_node_idxs;
    marked_check_node_idxs.reserve(num_check_nodes);

    for (const auto& bit_node : indexed_bit_nodes) 
    {
        int bit_node_idx = bit_node.first;
        const std::vector<int>& check_idxs = bit_node.second;
        int idx = find_available_index(check_idxs, marked_check_node_idxs);
        if (idx != -1) 
        {
            bits_to_remove.push_back(bit_node_idx);
            marked_check_node_idxs.push_back(idx);      // Mark the check node index as used
        }
    }
    // Sort in ascending order
    std::sort(bits_to_remove.begin(), bits_to_remove.end());

    return bits_to_remove;
}

// Determine positions of bits (including punctured and shortened bits)
// to be removed from key based on parity-check matrix for privacy maintenance.
std::vector<int> get_bits_positions_to_remove_rate_adapt(
    const H_matrix &matrix, 
    const H_matrix_params &mat_params
)
{
    auto& bit_nodes = matrix.bit_nodes;
    size_t num_bit_nodes = bit_nodes.size();
    size_t num_check_nodes = matrix.check_nodes.size();
    size_t num_punct_bits = mat_params.punctured_bits.size();
    size_t num_short_bits = mat_params.shortened_bits.size();

    // Indices of bits to be removed from key.
    std::vector<int> bits_to_remove;
    bits_to_remove.reserve(num_short_bits + num_check_nodes);   // num_check_nodes > num_punct_bits

    // Vector of pairs: (index of the bit node, corresponding check node indices).
    std::vector<std::pair<int, std::vector<int>>> indexed_bit_nodes;
    indexed_bit_nodes.reserve(num_bit_nodes - num_short_bits);  // Removing shortened bits from candidates.
    // For tracking check nodes that already processed.
    std::vector<int> marked_check_node_idxs;
    marked_check_node_idxs.reserve(num_check_nodes);
    size_t s = 0;
    size_t p = 0;
    int idx;
    for (int i = 0; i < num_bit_nodes; ++i) 
    {
        if (mat_params.shortened_bits[s] == i)
        {
            bits_to_remove.push_back(i);    // Delete all shortened bits.
            ++s;
        }
        else if (mat_params.punctured_bits[p] == i)
        {
            bits_to_remove.push_back(i);    // Delete all punctured bits.
            idx = find_available_index(bit_nodes[i], marked_check_node_idxs);
            if (idx != -1)  // Mark check node that attached to the current punctured bit node.
                marked_check_node_idxs.push_back(idx);  
            ++p;
        }
        else
            indexed_bit_nodes.emplace_back(i, bit_nodes[i]);
    }

    // Sort remaining bit nodes by the number of their connections (ascending order). 
    // Sort by column weight in H.
    sort(indexed_bit_nodes.begin(), indexed_bit_nodes.end(),
        [](const std::pair<int, std::vector<int>>& a, const std::pair<int, std::vector<int>>& b) 
        {
            return a.second.size() < b.second.size();
        });
    
    // Find the remaining m-p (num_check_nodes - num_punct_bits) bits.
    for (const auto& bit_node : indexed_bit_nodes) 
    {
        int bit_node_idx = bit_node.first;
        const std::vector<int>& check_idxs = bit_node.second;
        int idx = find_available_index(check_idxs, marked_check_node_idxs);
        if (idx != -1) 
        {
            bits_to_remove.push_back(bit_node_idx);
            marked_check_node_idxs.push_back(idx);      // Mark the check node index as used.
        }
    }
    // Sort in ascending order.
    std::sort(bits_to_remove.begin(), bits_to_remove.end());

    return bits_to_remove;
}

// Removing bits at positions specified in 'bits_to_remove' from the Alice and Bob keys.
void remove_bits(
    const std::vector<int> &bits_to_remove, 
    const std::vector<int> &array1, 
    const std::vector<int> &array2,
    std::vector<int> &array1_out,
    std::vector<int> &array2_out
)
{
    size_t btr_len = bits_to_remove.size();
    size_t new_arr_len = array1.size() - btr_len;
    array1_out.resize(new_arr_len);
    array2_out.resize(new_arr_len);

    size_t n = 0;
    size_t m = 0;
    for (int i = 0; i < array1.size(); ++i)
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

// Read sparse matrix from file in alist format 
// (https://rptu.de/channel-codes/matrix-file-formats).
H_matrix read_sparse_matrix_alist(const fs::path &matrix_path)
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
        for (const auto &l : line_vec)
        {
            std::istringstream iss(l);
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
        throw std::runtime_error("Wrong sparse alist matrix format: " + matrix_path.string());
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
    return matrix_out;
}

// Read sparse matrix from file in format:
// The first line contains the block length, N. The second line 
// defines the number of parity-checks, M. The third line defines
// the number of columns of the compressed parity-check matrix. 
// The following M lines are then the compressed parity-check matrix. 
// Each of the M rows contains the indices (1 ... N) of 1's in the 
// compressed row of parity-check matrix. If not all column entries 
// are used, the column is filled up with 0's.
H_matrix read_sparse_matrix_1(const fs::path &matrix_path)
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
        for (const auto &l : line_vec)
        {
            std::istringstream iss(l);
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

    if (vec_int.size() < 3)
    {
        throw std::runtime_error("Insufficient data in the file: " + matrix_path.string());
    }

    if (vec_int[0].size() != 1 || vec_int[1].size() != 1 || vec_int[2].size() != 1)
    {
        throw std::runtime_error("Wrong sparse matrix format: " + matrix_path.string());
    }

    size_t col_num = vec_int[0][0];             // n
    size_t row_num = vec_int[1][0];             // m
    size_t max_row_weight = vec_int[2][0];      // d^(c)_max

    size_t curr_line = 3;

    if (vec_int.size() < curr_line + row_num)
    {
        throw std::runtime_error("Insufficient data in the file: " + matrix_path.string());
    }

    bool max_weights_matched = false;
    size_t curr_row_weight = 0;

    H_matrix matrix_out{};
    // Filling the matrix of check nodes
    int bit_node_index{};
    try
    {
        matrix_out.check_nodes.resize(row_num);
        for (size_t i = 0; i < row_num; ++i)
        {
            curr_row_weight = vec_int[curr_line + i].size();
            for (size_t j = 0; j < curr_row_weight; ++j)
            {
                bit_node_index = vec_int[curr_line + i][j];
                if (bit_node_index < 0)
                {
                    throw std::runtime_error("Bit node index cannot be less than zero: " + std::to_string(bit_node_index) 
                    + ", row '"+ std::to_string(curr_line + i) + "'.");
                }
                if (curr_row_weight > max_row_weight)
                {
                    throw std::runtime_error("Actual weight '" + std::to_string(curr_row_weight) + "' of row '"
                    + std::to_string(curr_line + i) + "' exceeded the maximum specified weight '"
                    + std::to_string(max_row_weight) +"'."); 
                }
                if (bit_node_index != 0)
                {
                    matrix_out.check_nodes[i].push_back(bit_node_index - 1);
                }
                if (curr_row_weight == max_row_weight)
                {
                    max_weights_matched = true;
                }
            }
        }
    }
    catch (const std::exception &e)
    {
        fmt::print(stderr, fg(fmt::color::red), "An error occurred while creating 'check_nodes' matrix from file: {}\n", matrix_path.string());
        throw;
    }

    if (!max_weights_matched)
    {
        throw std::runtime_error("None of the row weights matched the specified maximum weight '" 
        + std::to_string(max_row_weight) + "'. File: " + matrix_path.string());
    }
    
    bool is_regular = true;
    for (size_t i = 0; i < matrix_out.check_nodes.size(); ++i)
    {
        if (matrix_out.check_nodes[0].size() != matrix_out.check_nodes[i].size())
        {
            is_regular = false;
            break;
        }   
    }

    matrix_out.is_regular = is_regular;
    try
    {
        matrix_out.bit_nodes = get_bit_nodes_from_check_nodes(matrix_out.check_nodes);
        if (matrix_out.bit_nodes.size() != col_num)
        {
            throw std::runtime_error("The actual number of bit nodes '"+ std::to_string(matrix_out.bit_nodes.size())
             + "' did not match the specified number '" + std::to_string(col_num) + "' of bit nodes.");
        }
    }
    catch(const std::exception& e)
    {
        fmt::print(stderr, fg(fmt::color::red), "An error occurred while creating 'bit_nodes' matrix from file: {}\n", matrix_path.string());
        throw;
    }
    return matrix_out;
}

// Read sparse matrix from file in format:
// The first line contains two numbers: the first is the block length (N)
// and the second is the number of parity-checks (M).The following M lines
// are then the compressed parity-check matrix. Each of the M rows contains 
// the indices (0 ... N-1) of 1's in the compressed row of parity-check matrix. 
// The next N lines contains the indices (0 ... M-1) of 1's in the compressed 
// column of parity-check matrix.
H_matrix read_sparse_matrix_2(const fs::path &matrix_path)
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
        for (const auto &l : line_vec)
        {
            std::istringstream iss(l);
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

    if (vec_int.size() < 2)
    {
        throw std::runtime_error("Insufficient data in the file: " + matrix_path.string());
    }

    if (vec_int[0].size() != 2)
    {
        throw std::runtime_error("Wrong sparse matrix format: " + matrix_path.string());
    }

    size_t col_num = vec_int[0][0];     // n
    size_t row_num = vec_int[0][1];     // m

    size_t curr_line = 1;

    if (vec_int.size() < curr_line + col_num + row_num)
    {
        throw std::runtime_error("Insufficient data in the file: " + matrix_path.string());
    }

    H_matrix matrix_out{};
    // Filling the matrix of check nodes
    int bit_node_index{};
    try
    {
        matrix_out.check_nodes.resize(row_num);
        for (size_t i = 0; i < row_num; ++i)
        {
            for (size_t j = 0; j < vec_int[curr_line + i].size(); ++j)
            {
                bit_node_index = vec_int[curr_line + i][j];
                if (bit_node_index < 0)
                {
                    throw std::runtime_error("Bit node index cannot be less than zero: " + std::to_string(bit_node_index) 
                    + ", row '"+ std::to_string(curr_line + i) + "'.");
                }
                matrix_out.check_nodes[i].push_back(bit_node_index);
            }
        }
    }
    catch (const std::exception &e)
    {
        fmt::print(stderr, fg(fmt::color::red), "An error occurred while creating 'check_nodes' matrix from file: {}\n", matrix_path.string());
        throw;
    }

    curr_line += row_num;

    int check_node_index{};
    try
    {
        matrix_out.bit_nodes.resize(col_num);
        for (size_t i = 0; i < col_num; ++i)
        {
            for (size_t j = 0; j < vec_int[curr_line + i].size(); ++j)
            {
                check_node_index = vec_int[curr_line + i][j];
                if (check_node_index < 0)
                {
                    throw std::runtime_error("Check node index cannot be less than zero: " + std::to_string(check_node_index) 
                    + ", row '"+ std::to_string(curr_line + i) + "'.");
                }
                matrix_out.bit_nodes[i].push_back(check_node_index);
            }
        }
    }
    catch (const std::exception &e)
    {
        fmt::print(stderr, fg(fmt::color::red), "An error occurred while creating 'bit_nodes' matrix from file: {}\n", matrix_path.string());
        throw;
    }

    bool is_regular = true;
    for (size_t i = 0; i < matrix_out.check_nodes.size(); ++i)
    {
        if (matrix_out.check_nodes[0].size() != matrix_out.check_nodes[i].size())
        {
            is_regular = false;
            break;
        }   
    }
    for (size_t i = 0; i < matrix_out.bit_nodes.size(); ++i)
    {
        if (matrix_out.bit_nodes[0].size() != matrix_out.bit_nodes[i].size())
        {
            is_regular = false;
            break;
        }   
    }

    matrix_out.is_regular = is_regular;
    return matrix_out;
}

// Read uncompressed sparse matrix from file.
H_matrix read_sparse_uncompressed_matrix(const fs::path &matrix_path)
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

    std::vector<std::vector<int8_t>> uncompr_mat;
    try
    {
        for (const auto &l : line_vec)
        {
            std::istringstream iss(l);
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
            uncompr_mat.push_back(numbers);
        }
    }
    catch (const std::exception &e)
    {
        fmt::print(stderr, fg(fmt::color::red), "An error occurred while parsing file: {}\n", matrix_path.string());
        throw;
    }

    for (size_t i = 0; i < uncompr_mat.size(); i++)
    {
        if (uncompr_mat[0].size() != uncompr_mat[i].size())
        {
            throw std::runtime_error("Different lengths of rows in a matrix. File: " + matrix_path.string());
        }
    }

    size_t col_num = uncompr_mat[0].size();    //n - bit_nodes number
    size_t row_num = uncompr_mat.size();       //m - check nodes number

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
            curr_weight += uncompr_mat[j][i];
        }
        if (curr_weight == 0)
        {
            throw std::runtime_error("Column '" + std::to_string(i + 1) + "' weight cannot be equal to zero. File: " + matrix_path.string());
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
        curr_weight = accumulate(uncompr_mat[i].begin(), uncompr_mat[i].end(), 0);
        if (curr_weight == 0)
        {
            throw std::runtime_error("Row '" + std::to_string(i + 1) + "' weight cannot be equal to zero. File: " + matrix_path.string());
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
    matrix_out.bit_nodes = get_bit_nodes(uncompr_mat, bit_nodes_weight);
    matrix_out.check_nodes = get_check_nodes(uncompr_mat, check_nodes_weight);
    matrix_out.is_regular = is_regular;
    return matrix_out;
}

// Generates Alice's key.
void fill_random_bits(
    XoshiroCpp::Xoshiro256PlusPlus &prng,
    std::vector<int> &bit_array
)
{
    std::uniform_int_distribution<int> distribution(0, 1);

    // Generate random bits and fill the vector
    for (size_t i = 0; i < bit_array.size(); ++i)
    {
        bit_array[i] = distribution(prng);
    }
}

// Generates Bob's key by making errors in Alice's key. Generates the 
// exact number of errors in the key and returns the exact QBER.
double inject_errors(
    XoshiroCpp::Xoshiro256PlusPlus &prng, 
    const std::vector<int> &bit_array, 
    double QBER, 
    std::vector<int> &bit_array_with_errors_out
)
{
    size_t array_length = bit_array.size();
    size_t num_errors = static_cast<size_t>(static_cast<double>(array_length) * QBER);

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
    return static_cast<double>(num_errors) / static_cast<double>(array_length);
}

// Computes the key syndrome using parity check matrix.
void calculate_syndrome(
    const std::vector<int> &bit_array,
    const H_matrix &matrix, 
    std::vector<int> &syndrome_out
)
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
void threshold_matrix(
    std::vector<std::vector<double>> &matrix, 
    const double &msg_threshold
)
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

// Returns the second-order neighbors for each bit node of parity check matrix.
std::vector<second_order_neighbors> get_second_order_neighbors(const H_matrix &matrix) 
{
    auto& check_nodes = matrix.check_nodes;
    auto& bit_nodes = matrix.bit_nodes;
    size_t num_bit_nodes = bit_nodes.size();
    std::vector<second_order_neighbors> N2_order_list(num_bit_nodes);

    for (int i = 0; i < num_bit_nodes; ++i) 
    {
        N2_order_list[i].bit_node_idx = i;

        // Set union: add all bit nodes connected to adjacent check nodes of the current bit node
        for (int check_node_idx : bit_nodes[i]) 
        {
            // Note: set can only store unique values, we can't add the same value multiple times
            N2_order_list[i].neighbors.insert(check_nodes[check_node_idx].begin(), check_nodes[check_node_idx].end());
        }

        // Delete index of the current bit node
        N2_order_list[i].neighbors.erase(i);
    }
    return N2_order_list;
}

// Returns puncturable bits (symbol nodes) indexes using untainted puncturing 
// (https://arxiv.org/abs/1103.6149). Returns positions of punctured bits for  
// the given matrix (maximum possible number of positions).
std::vector<int> select_punctured_bits_untainted(
    XoshiroCpp::Xoshiro256PlusPlus &prng,
    const H_matrix &matrix
) 
{
    size_t num_bit_nodes = matrix.bit_nodes.size();
    std::vector<second_order_neighbors> N2_order_list = get_second_order_neighbors(matrix);
    std::vector<int> punct_nodes;   // Set of punctured bit nodes
    std::set<int> X;                // Set of untainted bit nodes

    // Initialize X with all bit nodes from 0 to n-1
    for (int i = 0; i < num_bit_nodes; ++i) 
    {
        X.insert(i);
    }

    while (!X.empty()) 
    {
        // Step 1: Search for minimum number of second-order neighbors among the nodes in X
        size_t min_n = num_bit_nodes;
        for (int i : X) 
        {
            size_t cur_n = 0;
            // Calculation of the intersection size of sets N2_order_list[i].neighbors and X
            for (int neighbor : N2_order_list[i].neighbors) 
            {
                if (X.count(neighbor)) 
                    cur_n++;
            }
            if (cur_n < min_n) 
                min_n = cur_n;
        }

        // Step 2: Collect candidates with minimum number of second-order neighbors
        std::vector<int> candidate_nodes;
        for (int i : X) 
        {
            size_t cur_n = 0;
            for (int neighbor : N2_order_list[i].neighbors) 
            {
                if (X.count(neighbor))
                    cur_n++;
            }
            if (cur_n == min_n)
                candidate_nodes.push_back(i);
        }

        // Step 3: Select a random bit node from candidates
        std::uniform_int_distribution<size_t> distribution(0, candidate_nodes.size() - 1);
        size_t idx = distribution(prng);
        int chosen_idx = candidate_nodes[idx];
        const second_order_neighbors& punct_bit = N2_order_list[chosen_idx];

        // Add chosen node to set
        punct_nodes.push_back(punct_bit.bit_node_idx);

        // Remove the selected node from X
        X.erase(punct_bit.bit_node_idx);

        // Remove all its second-order neighbors from X
        for (int neighbor : punct_bit.neighbors) 
        {
            X.erase(neighbor);
        }
    }
    return punct_nodes;
}

// Gets a list of indexes of “untainted” bits intended for puncturing.
// This function attempts to read a previously saved list of indexes of 
// untainted bits from a file with the extension `.untp` located next to the
// specified matrix (extension `.mtrx`). If the file exists and contains valid
// indexes, they are used directly. Otherwise: a new list is generated using 
// select_punctured_bits_untainted() and saved to the .untp file.
std::vector<int> get_punctured_bits_untainted(
    const fs::path &matrix_path,
    XoshiroCpp::Xoshiro256PlusPlus &prng,
    const H_matrix &matrix
)
{
    fs::path punct_bits_unt_path = matrix_path;
    punct_bits_unt_path.replace_extension(".untp");

    std::vector<int> punct_bits_unt{};
    std::ifstream infile(punct_bits_unt_path);
    if (infile.is_open()) 
    {
        std::string line{};
        if (std::getline(infile, line)) 
        {
            std::istringstream iss(line);
            std::copy(std::istream_iterator<int>(iss), std::istream_iterator<int>(), std::back_inserter(punct_bits_unt));
        }
        infile.close();
    }

    for (size_t i = 0; i < punct_bits_unt.size(); ++i)
    {
        if (punct_bits_unt[i] < 0 || punct_bits_unt[i] >= matrix.bit_nodes.size())
            throw std::runtime_error(fmt::format("The punctured bit index '{}' is out of range [0,{}]. File: {}",
                punct_bits_unt[i],  (matrix.bit_nodes.size() - 1), punct_bits_unt_path.string()));
    }

    if (punct_bits_unt.empty()) 
    {
        fmt::print(fg(fmt::color::purple), "WARNING: No file with punctured untainted bits found: {}"
        " \nThis file will be automatically created. Wait...\n", punct_bits_unt_path.string());

        punct_bits_unt = select_punctured_bits_untainted(prng, matrix);
        std::ofstream outfile(punct_bits_unt_path);
        if (outfile.is_open()) 
        {
            std::copy(punct_bits_unt.begin(), punct_bits_unt.end(), std::ostream_iterator<int>(outfile, " "));
            outfile.close();
            fmt::print(fg(fmt::color::purple), "File created successfully.\n");
        } 
        else 
            throw std::runtime_error(fmt::format("Unable to open file for writing: {}", punct_bits_unt_path.string()));
    }

    return punct_bits_unt;
}

// Determines the number and positions of punctured and shortened bits for the original
// (pre-built) code (with rate R0) as a function of QBER, fraction δ, and efficiency (f_EC).  
// The indices of punctured and shortened bits are written in ascending order.
// https://arxiv.org/abs/1007.1616
H_matrix_params adapt_code_rate(
    XoshiroCpp::Xoshiro256PlusPlus &prng, 
    const H_matrix &matrix,
    double QBER,
    double delta,
    double efficiency
) 
{
    // Shannon binary entropy of the QBER
    double h_b = -QBER * std::log2(QBER) - (1. - QBER) * std::log2(1. - QBER);

    // Optimal code rate
    double optimal_R = 1. - efficiency * h_b;

    size_t num_bit_nodes = matrix.bit_nodes.size();
    size_t num_check_nodes = matrix.check_nodes.size();

    double original_R = 1. - static_cast<double>(num_check_nodes) / static_cast<double>(num_bit_nodes);

    // Optimal values for puncturing and shortening, p and s respectively
    int num_short_bits = static_cast<int>(std::ceil((original_R - optimal_R * (1. - delta)) * static_cast<double>(num_bit_nodes)));
    int num_punct_bits = static_cast<int>(delta * static_cast<double>(num_bit_nodes) - static_cast<double>(num_short_bits));

    H_matrix_params m_param{};
    double min_R = (original_R - delta)/(1. - delta);
    double max_R = (original_R)/(1. - delta);
    if (num_short_bits <= 0 || num_punct_bits <= 0)
    {
        fmt::print(fg(fmt::color::purple), "WARNING: R0 = {:.3f}, QBER = {:.4f}, delta = {:.3f}, f_EC = {:.3f}. "
        "Adapted code rate R = {:.3f} beyond the achievable rate range: Rmin = {:.3f}, Rmax = {:.3f}. "
        "This parameters will not be used in simulations.\n", 
        original_R, QBER, delta, efficiency, optimal_R, min_R, max_R);
        return m_param;
    }
    
    std::vector<int> bit_positions(num_bit_nodes);
    if (CFG.ENABLE_UNTAINTED_PUNCTURING)
    {
        const std::vector<int>& punct_bits = matrix.punctured_bits_untainted;
        if (num_punct_bits > punct_bits.size())
        {
            fmt::print(fg(fmt::color::purple), "WARNING: R0 = {:.3f}, QBER = {:.4f}, delta = {:.3f}, f_EC = {:.3f}, R = {:.3f}"
                ", Rmin = {:.3f}, Rmax = {:.3f}. The calculated number of punctured bits ({}) exceeds the number of bits "
                "produced by untainted algorithm ({}). These parameters will not be used in simulations.\n", 
                original_R, QBER, delta, efficiency, optimal_R, min_R, max_R, num_punct_bits, punct_bits.size()); 
            return m_param;
        }
        m_param.punctured_bits.assign(punct_bits.begin(), punct_bits.begin() + num_punct_bits);
    }
    else
    {
        for (int i = 0; i < num_bit_nodes; ++i)
        {
            bit_positions[i] = i;
        }
        std::shuffle(bit_positions.begin(), bit_positions.end(), prng);
        m_param.punctured_bits.assign(bit_positions.begin(), bit_positions.begin() + num_punct_bits);
    }
    // Sort by index (ascending order).
    sort(m_param.punctured_bits.begin(), m_param.punctured_bits.end(),
        [](const int& a, const int& b) 
        {
            return a < b;
        });

    for (int i = 0; i < num_bit_nodes; ++i)
    {
        bit_positions[i] = i;
    }

    std::vector<int> remaining_bits(num_bit_nodes - num_punct_bits);
    std::set_difference(
        bit_positions.begin(), bit_positions.end(), 
        m_param.punctured_bits.begin(), m_param.punctured_bits.end(), 
        remaining_bits.begin()
    );
    
    std::shuffle(remaining_bits.begin(), remaining_bits.end(), prng);
    m_param.shortened_bits.reserve(num_short_bits);
    m_param.shortened_bits.assign(remaining_bits.begin(), remaining_bits.begin() + num_short_bits);
    sort(m_param.shortened_bits.begin(), m_param.shortened_bits.end(),
        [](const int& a, const int& b) 
        {
            return a < b;
        });
    
    m_param.delta = delta;
    m_param.efficiency = efficiency;
    m_param.shortened_fraction = static_cast<double>(num_short_bits)/static_cast<double>(num_bit_nodes);
    m_param.punctured_fraction = static_cast<double>(num_punct_bits)/static_cast<double>(num_bit_nodes);
    m_param.adapted_code_rate = (static_cast<double>(num_bit_nodes - num_check_nodes - num_short_bits))/
        (static_cast<double>(num_bit_nodes - num_punct_bits - num_short_bits));

    return m_param;
}