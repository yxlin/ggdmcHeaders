#pragma once
#include "type_aliases.h"
#include <RcppArmadillo.h>

inline bool is_null_or_empty(const Rcpp::RObject &x)
{
    if (x.isNULL())
        return true;
    if (TYPEOF(x) == NILSXP)
        return true;
    if (Rf_length(x) == 0)
        return true;
    return false;
}

inline std::string get_model_type(const Rcpp::S4 &model)
{
    Rcpp::CharacterVector cls = model.attr("class");
    std::string cname = Rcpp::as<std::string>(cls[0]);

    if (cname == "model_lba")
        return "lba";
    if (cname == "model_fastdm")
        return "fastdm";
    if (cname == "model_hyper")
        return "hyper";
    if (cname == "model_cdm")
        return "cdm";

    Rcpp::stop("Unsupported model class: " + cname);
}

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

// Convertors should handle NULL by returning empty structures.
inline std::map<std::string, double> constants_to_map(const Rcpp::RObject &x)
{
    if (is_null_or_empty(x))
        return {};
    Rcpp::List lst(x);
    std::map<std::string, double> out;
    Rcpp::CharacterVector nms = lst.names();
    for (int i = 0; i < lst.size(); ++i)
    {
        std::string k = Rcpp::as<std::string>(nms[i]);
        double v = Rcpp::as<double>(lst[i]);
        out.emplace(std::move(k), v);
    }
    return out;
}

// NULL-OK: character vector (or list of scalars) -> std::vector<std::string>
inline strVec accumulators_to_vec(const Rcpp::RObject &input)
{

    if (input.isNULL())
    {
        return {};
    }

    if (TYPEOF(input) == STRSXP)
    {
        Rcpp::CharacterVector cv(input);
        return Rcpp::as<strVec>(cv); // length 0 -> {}
    }
    if (TYPEOF(input) == VECSXP)
    {
        Rcpp::List lst(input);
        strVec out;
        out.reserve(lst.size());

        for (R_xlen_t i = 0; i < lst.size(); ++i)
        {
            Rcpp::RObject elem = lst[i];
            if (elem.isNULL() || Rf_length(elem) == 0)
            {
                Rcpp::stop(
                    "`accumulators[[%ld]]` is NULL/empty; expected a string.",
                    (long)i + 1);
            }
            std::string s = Rcpp::as<std::string>(elem);
            out.emplace_back(std::move(s));
        }
        return out;
    }
    Rcpp::stop("`accumulators` must be NULL, a character vector, or a list of "
               "strings.");
    return {}; // not reached
}

// NULL-OK: named nested list -> map<string, map<string,string>>
inline std::map<std::string, std::map<std::string, std::string>>
match_map_to_map(const Rcpp::RObject &input)
{
    if (input.isNULL())
        return {};
    if (TYPEOF(input) != VECSXP)
        Rcpp::stop("`match_map` must be NULL or a named list.");
    return nested_list_to_map(Rcpp::List(input)); // your existing converter
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

inline arma::mat r_mat_to_arma_mat_any(const Rcpp::RObject &obj)
{
    if (!Rf_isMatrix(obj))
        Rcpp::stop("Expected a matrix (base R matrix).");

    const SEXPTYPE tp = TYPEOF(obj);
    const int nr = Rf_nrows(obj);
    const int nc = Rf_ncols(obj);

    if (nr == 0 || nc == 0)
        return arma::mat(); // empty OK

    if (tp == REALSXP)
    {
        Rcpp::NumericMatrix nm(obj);
        arma::mat out(nr, nc);
        std::copy(nm.begin(), nm.end(), out.memptr()); // col-major→col-major
        return out;
    }
    else if (tp == INTSXP || tp == LGLSXP)
    {
        Rcpp::IntegerMatrix im(obj);
        arma::mat out(nr, nc);
        for (int j = 0; j < nc; ++j)
            for (int i = 0; i < nr; ++i)
            {
                int v = im(i, j);
                out(i, j) =
                    (v == NA_INTEGER) ? NA_REAL : static_cast<double>(v);
            }
        return out;
    }

    Rcpp::stop("Expected a numeric or integer matrix.");
}

inline std::vector<std::vector<double>>
list_to_std_mat(const Rcpp::List &data_r)
{
    std::vector<std::string> cell_names = data_r.names();
    std::size_t ncell = cell_names.size();

    std::vector<std::vector<double>> data_cpp;
    data_cpp.reserve(ncell);

    for (std::size_t cell_idx = 0; cell_idx < ncell; ++cell_idx)
    {
        // Each entry is assumed to be a NumericVector
        Rcpp::NumericVector rt_r = data_r[cell_idx]; // could RT or response 1/0
        std::vector<double> rt(rt_r.begin(), rt_r.end());
        data_cpp.push_back(std::move(rt));
    }
    return data_cpp;
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
