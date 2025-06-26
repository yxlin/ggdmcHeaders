#pragma once
#include <RcppArmadillo.h>

using bool3D = std::vector<std::vector<std::vector<bool>>>;
using uint2D = std::vector<std::vector<unsigned int>>;

/* ------------- Type casting ------------- */
template <typename T>
std::map<std::string, std::vector<T>> list_to_map(const Rcpp::List &input_list)
{
    // Must stay in the hpp file, bcuz it's a template.
    std::map<std::string, std::vector<T>> result;
    Rcpp::CharacterVector names = input_list.names();

    for (int i = 0; i < input_list.size(); ++i)
    {
        std::string key = Rcpp::as<std::string>(names[i]);
        std::vector<T> values = Rcpp::as<std::vector<T>>(input_list[i]);
        result[key] = values;
    }

    return result;
}

inline std::map<std::string, std::map<std::string, std::string>>
nested_list_to_map(const Rcpp::List &input_list)
{
    std::map<std::string, std::map<std::string, std::string>> result;
    Rcpp::CharacterVector outer_names = input_list.names();
    size_t n_input_list = input_list.size();
    for (size_t i = 0; i < n_input_list; ++i)
    {
        std::string outer_key = Rcpp::as<std::string>(outer_names[i]);
        Rcpp::List inner_list = Rcpp::as<Rcpp::List>(input_list[i]);
        size_t n_inner_list = inner_list.size();

        Rcpp::CharacterVector inner_names = inner_list.names();
        std::map<std::string, std::string> inner_map;

        for (size_t j = 0; j < n_inner_list; ++j)
        {
            std::string inner_key = Rcpp::as<std::string>(inner_names[j]);
            std::string value = Rcpp::as<std::string>(inner_list[j]);

            inner_map[inner_key] = value;
        }

        result[outer_key] = inner_map;
    }
    return result;
}

inline std::map<std::string, double>
constants_to_map(const Rcpp::List &input_list)
{
    std::map<std::string, double> out;
    Rcpp::CharacterVector constants_names = input_list.names();
    size_t n_constant = constants_names.size();

    for (size_t i = 0; i < n_constant; ++i)
    {
        std::string key = Rcpp::as<std::string>(constants_names[i]);
        double value = input_list[key];
        out[key] = value;
    }
    return out;
}

template <typename T, typename R>
inline std::vector<std::vector<R>> r_mat_to_std_mat(const T &input)
{
    if (input.nrow() == 0 || input.ncol() == 0)
    {
        return {}; // Return empty vector if input is empty
    }

    size_t nrows = input.nrow();
    size_t ncols = input.ncol();

    std::vector<std::vector<R>> output(nrows, std::vector<R>(ncols));

    for (size_t i = 0; i < nrows; ++i)
    {
        for (size_t j = 0; j < ncols; ++j)
        {
            output[i][j] = static_cast<R>(input(i, j));
        }
    }

    return output;
}

/* ------------- design light type casting ------------- */
/*--R to std (only used in ggdmcLikelihood and lbaModel for
 * new_design_light)--
 new_design_light_rt_model
 */
inline std::vector<std::vector<std::vector<bool>>>
R_ucube_to_std_ucube(const Rcpp::LogicalVector &input)
{
    // Extract and validate dimensions
    Rcpp::IntegerVector dims = input.attr("dim");
    if (dims.size() != 3)
    {
        Rcpp::stop("Input must be a 3D array.");
    }

    size_t n_cell = dims[0];
    size_t n_param = dims[1];
    size_t n_acc = dims[2];
    size_t n_input = input.size(); // Casting it to size_t
    if (n_input != n_cell * n_param * n_acc)
    {
        Rcpp::stop("Input size does not match specified dimensions.");
    }

    // Initialize 3D std::vector
    std::vector<std::vector<std::vector<bool>>> out(
        n_cell,
        std::vector<std::vector<bool>>(n_param, std::vector<bool>(n_acc)));

    // Fill it by matching R's column-major order
    for (size_t acc = 0; acc < n_acc; ++acc)
    {
        for (size_t param = 0; param < n_param; ++param)
        {
            for (size_t cell = 0; cell < n_cell; ++cell)
            {
                size_t index = cell + n_cell * param + n_cell * n_param * acc;
                out[cell][param][acc] = input[index];
            }
        }
    }

    return out;
}
