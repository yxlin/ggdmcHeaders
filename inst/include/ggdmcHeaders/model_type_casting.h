#pragma once
#include <RcppArmadillo.h>

/* ------------- Type casting ------------- */
inline Rcpp::LogicalVector
std_ucube_to_R_ucube(const std::vector<std::vector<std::vector<bool>>> &input)
{
    // the reverse is in lba_type_casting.h
    if (input.empty() || input[0].empty() || input[0][0].empty())
    {
        return Rcpp::LogicalVector();
    }

    size_t n_cell = input.size();
    size_t n_param = input[0].size();
    size_t n_acc = input[0][0].size();

    Rcpp::LogicalVector out(n_cell * n_param * n_acc);
    for (size_t acc = 0; acc < n_acc; ++acc)
    {
        for (size_t param = 0; param < n_param; ++param)
        {
            for (size_t cell = 0; cell < n_cell; ++cell)
            {
                // Correct R-style (column-major) indexing
                size_t index = cell + n_cell * param + n_cell * n_param * acc;
                out[index] = input[cell][param][acc];
            }
        }
    }
    out.attr("dim") = Rcpp::IntegerVector::create(n_cell, n_param, n_acc);
    return out;
}

inline Rcpp::IntegerMatrix
std_umat_to_R_int_mat(const std::vector<std::vector<unsigned int>> &input)
{
    if (input.empty())
    {
        return Rcpp::IntegerMatrix(0, 0);
    }
    size_t nrows = input.size();
    size_t ncols = input[0].size();

    Rcpp::IntegerMatrix out(nrows, ncols);
    for (size_t i = 0; i < nrows; ++i)
    {
        if (input[i].size() != ncols)
        {
            Rcpp::stop("All inner vectors must have the same length for matrix "
                       "conversion.");
        }
        out(i, Rcpp::_) = Rcpp::IntegerVector(input[i].begin(), input[i].end());
    }
    return out;
}

///////////////////////////////////////////////////
/* ------------- Design helpers  ------------- */
///////////////////////////////////////////////////
inline std::shared_ptr<design::design_class> new_design(const Rcpp::S4 &model_r)
{
    // model2R has called design.h and common_type_casting.h
    //
    // 1. The user enterd elements
    Rcpp::List parameter_map_r = model_r.slot("parameter_map");
    Rcpp::CharacterVector accumulators_r = model_r.slot("accumulators");
    Rcpp::List factors_r = model_r.slot("factors");
    Rcpp::List match_map_r = model_r.slot("match_map");
    Rcpp::List constants_r = model_r.slot("constants");

    // Convert to C++ types
    std::map<std::string, std::vector<std::string>> parameter_map =
        list_to_map<std::string>(parameter_map_r);

    std::vector<std::string> accumulators =
        Rcpp::as<std::vector<std::string>>(accumulators_r);

    std::map<std::string, std::vector<std::string>> factors =
        list_to_map<std::string>(factors_r);

    std::map<std::string, std::map<std::string, std::string>> match_map =
        nested_list_to_map(match_map_r);

    std::map<std::string, double> constants = constants_to_map(constants_r);

    return std::make_shared<design::design_class>(
        match_map, parameter_map, factors, accumulators, constants);
}
