#pragma once
// #include "common_type_casting.h"
#include "likelihood.h" // incluce design_light.h
#include <RcppArmadillo.h>

///////////////////////////////////////////////////
/* ------------- Design helpers  ------------- */
///////////////////////////////////////////////////
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

inline Rcpp::LogicalMatrix
std_matrix_to_R_matrix(const std::vector<std::vector<bool>> &m)
{
    const std::size_t nrow = m.size();
    const std::size_t ncol = nrow ? m[0].size() : 0;

    // Rectangularity check
    for (std::size_t i = 0; i < nrow; ++i)
    {
        if (m[i].size() != ncol)
        {
            Rcpp::stop("All inner vectors must have the same length for matrix "
                       "conversion.");
        }
    }

    Rcpp::LogicalMatrix out(nrow, ncol);

    // Fill column-major for R
    for (std::size_t j = 0; j < ncol; ++j)
    {
        for (std::size_t i = 0; i < nrow; ++i)
        {
            out(i, j) = m[i][j] ? 1 : 0;
        }
    }
    return out;
}

///////////////////////////////////////////////////
/* ------------- Design helpers  ------------- */
///////////////////////////////////////////////////
// Main implementation
inline std::shared_ptr<design::design_class>
create_design_impl(const Rcpp::S4 &model_r)
{

    // Rcpp::Rcout << "In create_design_impl: likelihood_type_casting.h\n";

    Rcpp::List parameter_map_r = model_r.slot("parameter_map");
    auto parameter_map = list_to_map<std::string>(parameter_map_r);

    std::map<std::string, double> constants =
        constants_to_map(model_r.slot("constants"));

    auto cell_names = Rcpp::as<strVec>(model_r.slot("cell_names"));
    auto model_boolean = R_ucube_to_std_ucube(model_r.slot("model_boolean"));
    std::string ms = get_model_type(model_r);

    auto parameter_x_condition_names =
        Rcpp::as<strVec>(model_r.slot("parameter_x_condition_names"));

    Rcpp::RObject accumulators_obj = model_r.slot("accumulators");
    strVec accumulators = accumulators_to_vec(accumulators_obj);

    return std::make_shared<design::design_class>(
        parameter_map, accumulators, cell_names, parameter_x_condition_names,
        constants, model_boolean, ms);
}

// Overload for direct model input
inline std::shared_ptr<design::design_class> new_design(const Rcpp::S4 &model_r)
{
    return create_design_impl(model_r);
}

// Overload for DMI input (extracts model first)
inline std::shared_ptr<design::design_class>
new_design_light(const Rcpp::S4 &dmi)
{
    // Rcpp::Rcout << "In new_design_light: likelihood_type_casting.h\n";

    Rcpp::S4 model_r = dmi.slot("model");
    return create_design_impl(model_r);
}

///////////////////////////////////////////////////
/* ------------- Likelihood  ------------- */
///////////////////////////////////////////////////
inline void cdm_stop_with_usage(const std::string &what,
                                const std::string &detail = "")
{
    Rcpp::stop("CDM requires '%s' in `dmi`.\n"
               "%s\n"
               "Did you create `dmi` with all CDM arguments? Example:\n"
               "  sub_dmis <- BuildDMI(dat$responses, model,\n"
               "                       q_matrix = Q,\n"
               "                       prior_pi = pi_uniform,\n"
               "                       rule = \"DINA\")",
               what.c_str(), detail.empty() ? "" : detail.c_str());
}

struct CDMInputs
{
    arma::mat q_matrix;
    std::vector<double> prior_pi;
    std::string rule;
};

// Validates presence, non-NULL, types, and 2^K length for prior
inline CDMInputs require_cdm_inputs(const Rcpp::S4 &dmi)
{
    // --- q_matrix ---
    if (!dmi.hasSlot("q_matrix"))
        cdm_stop_with_usage("q_matrix", "Slot 'q_matrix' is missing.");
    Rcpp::RObject q_obj = dmi.slot("q_matrix");
    if (q_obj.isNULL())
        cdm_stop_with_usage("q_matrix", "Slot 'q_matrix' is NULL.");
    if (!(TYPEOF(q_obj) == REALSXP || TYPEOF(q_obj) == INTSXP) ||
        !Rf_isMatrix(q_obj))
        cdm_stop_with_usage(
            "q_matrix", "Slot 'q_matrix' must be a numeric/integer matrix.");

    arma::mat q_mat = r_mat_to_arma_mat(q_obj);

    // --- prior_pi ---
    if (!dmi.hasSlot("prior_pi"))
        cdm_stop_with_usage("prior_pi", "Slot 'prior_pi' is missing.");
    Rcpp::RObject pi_obj = dmi.slot("prior_pi");
    if (pi_obj.isNULL())
        cdm_stop_with_usage("prior_pi", "Slot 'prior_pi' is NULL.");
    if (!((TYPEOF(pi_obj) == REALSXP || TYPEOF(pi_obj) == INTSXP) &&
          !Rf_isMatrix(pi_obj)))
        cdm_stop_with_usage(
            "prior_pi",
            "Slot 'prior_pi' must be a numeric vector (length typically 2^K).");

    std::vector<double> pi_prior = Rcpp::as<std::vector<double>>(pi_obj);

    // Sanity: 2^K matches length(prior_pi)
    {
        std::size_t K = static_cast<std::size_t>(q_mat.n_cols);
        std::size_t expected =
            (K >= (sizeof(std::size_t) * 8)) ? 0u : (1ull << K);
        if (expected > 0 && pi_prior.size() != expected)
        {
            Rcpp::stop("CDM prior length mismatch: got prior_pi length = %d, "
                       "but expected 2^K = %d (K = %d).\n"
                       "Check your `prior_pi` and `q_matrix`.\n"
                       "Example of a uniform prior: rep(1/(2^K), 2^K).",
                       static_cast<int>(pi_prior.size()),
                       static_cast<int>(expected), static_cast<int>(K));
        }
    }

    // --- rule ---
    if (!dmi.hasSlot("rule"))
        cdm_stop_with_usage("rule", "Slot 'rule' is missing.");
    Rcpp::RObject rule_obj = dmi.slot("rule");
    if (rule_obj.isNULL())
        cdm_stop_with_usage("rule", "Slot 'rule' is NULL.");
    if (TYPEOF(rule_obj) != STRSXP || Rf_length(rule_obj) != 1)
        cdm_stop_with_usage(
            "rule",
            "Slot 'rule' must be a single string: \"DINA\" or \"DINO\".");

    std::string rule = Rcpp::as<std::string>(rule_obj);
    if (!(rule == "DINA" || rule == "DINO"))
        cdm_stop_with_usage("rule",
                            "Slot 'rule' must be \"DINA\" or \"DINO\".");

    return CDMInputs{std::move(q_mat), std::move(pi_prior), std::move(rule)};
}

inline std::vector<std::string> list_names_to_std(const Rcpp::List &lst)
{
    Rcpp::RObject nm = lst.names();
    if (nm.isNULL())
        return {};
    return Rcpp::as<std::vector<std::string>>(nm);
}

std::vector<std::shared_ptr<likelihood::likelihood_class>>
new_likelihoods(const Rcpp::List &dmis)
{
    const unsigned int n_subject = dmis.size(); // or n_school
    std::vector<std::shared_ptr<likelihood::likelihood_class>> out(n_subject);

    for (size_t subject_idx = 0; subject_idx < n_subject; ++subject_idx)
    {
        Rcpp::S4 dmi = dmis[subject_idx];
        Rcpp::List data_r = dmi.slot("data");
        Rcpp::S4 model_r = dmi.slot("model");

        auto d_ptr = new_design_light(dmi);
        std::string ms = get_model_type(model_r);
        auto cell_names = list_names_to_std(data_r);
        auto data_cpp = list_to_std_mat(data_r);

        if (ms == "cdm")
        {
            CDMInputs in = require_cdm_inputs(dmi);
            out[subject_idx] = std::make_shared<likelihood::likelihood_class>(
                d_ptr, data_cpp, cell_names, ms, in.q_matrix, in.prior_pi,
                in.rule);
        }
        else if (ms == "lba" || ms == "fastdm")
        {
            // Requires is_positive_drift
            if (!dmi.hasSlot("is_positive_drift"))
                Rcpp::stop("DMI missing required slot: is_positive_drift");
            auto is_positive_drift =
                Rcpp::as<std::vector<bool>>(dmi.slot("is_positive_drift"));
            out[subject_idx] = std::make_shared<likelihood::likelihood_class>(
                d_ptr, data_cpp, cell_names, ms, is_positive_drift);
        }
        else if (ms == "hyper")
        {
            Rcpp::stop("Hyper model not supported in new_likelihoods().");
        }
        else
        {
            Rcpp::stop("Unknown model type: %s", ms);
        }
    }
    return out;
}

/**
 * Creates a likelihood object from DMI (Data Model Interface)
 *
 * @param dmi An S4 object containing model and data specifications
 * @param p_prior Optional shared pointer to prior object (default nullptr)
 * @return Shared pointer to likelihood_class object
 */
std::shared_ptr<likelihood::likelihood_class>
new_likelihood(const Rcpp::S4 &dmi,
               const std::shared_ptr<prior::prior_class> &p_prior = nullptr)
{
    Rcpp::RObject data_obj = dmi.slot("data");
    Rcpp::S4 model_r = dmi.slot("model");
    std::string ms = get_model_type(model_r);

    arma::mat theta_data;
    Rcpp::List data_r;
    std::vector<std::string> cell_names;

    if (ms == "hyper")
    {
        if (TYPEOF(data_obj) != REALSXP)
        {
            Rcpp::stop("For hyper model, `dmi@data` must be a numeric matrix.");
        }
        theta_data = Rcpp::as<arma::mat>(data_obj);
    }

    if (TYPEOF(data_obj) != VECSXP)
    {
        Rcpp::stop("For %s model, `dmi@data` must be a list of cells.", ms);
    }
    data_r = Rcpp::as<Rcpp::List>(data_obj);
    cell_names = list_names_to_std(data_r);

    auto d_ptr = new_design_light(dmi);

    if (p_prior)
    {
        Rcpp::stop("p_prior can only be used with hyper model, got %s", ms);
    }

    if (ms == "cdm")
    {
        auto data_cpp = list_to_std_mat(data_r);
        CDMInputs in = require_cdm_inputs(dmi);

        return std::make_shared<likelihood::likelihood_class>(
            d_ptr, data_cpp, cell_names, ms, in.q_matrix, in.prior_pi, in.rule);
    }
    else if (ms == "lba" || ms == "fastdm")
    {
        auto data_cpp = list_to_std_mat(data_r);
        if (!dmi.hasSlot("is_positive_drift"))
            Rcpp::stop("DMI missing required slot: is_positive_drift");
        auto is_positive_drift =
            Rcpp::as<std::vector<bool>>(dmi.slot("is_positive_drift"));

        return std::make_shared<likelihood::likelihood_class>(
            d_ptr, data_cpp, cell_names, ms, is_positive_drift);
    }
    Rcpp::stop("Unknown model type: %s", ms);
}
