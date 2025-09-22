#pragma once
#include "common_type_casting.h"
#include "likelihood_type_casting.h"
#include <RcppArmadillo.h>
#include <optional>

inline Rcpp::IntegerMatrix wrap_imat(const arma::imat &M)
{
    Rcpp::IntegerMatrix out(M.n_rows, M.n_cols);
    std::copy(M.begin(), M.end(),
              out.begin()); // Armadillo is column-major; R is too
    return out;
}

inline Rcpp::LogicalMatrix
wrap_uchar_as_logical(const arma::Mat<unsigned char> &A)
{
    Rcpp::LogicalMatrix out(A.n_rows, A.n_cols);
    auto *p = out.begin();
    for (arma::uword i = 0; i < A.n_elem; ++i)
        p[i] = A[i] != 0;
    return out;
}

struct SimulationResultsCDM
{
    // One entry per condition, in the same order as simulate_conditions()
    std::vector<std::string> conditions; // length = n_condition
    std::vector<arma::imat> Y;           // each is N × J (responses 0/1)
    std::vector<arma::Mat<unsigned char>> Alpha; // each is N × K (mastery 0/1)

    // Optional metadata (can help downstream checks)
    unsigned int J = 0; // n_item
    unsigned int K = 0; // n_skill

    void reserve(size_t n_cond)
    {
        conditions.reserve(n_cond);
        Y.reserve(n_cond);
        Alpha.reserve(n_cond);
    }

    void add_condition(const std::string &cond_name,
                       const arma::imat &Y_block,                   // N × J
                       const arma::Mat<unsigned char> &Alpha_block, // N × K
                       unsigned int n_item, unsigned int n_skill)
    {
        if (J == 0 && K == 0)
        {
            J = n_item;
            K = n_skill;
        }
        // (Optional) sanity checks:
        if (Y_block.n_cols != J)
            throw std::runtime_error("Y.n_cols != J");
        if (Alpha_block.n_cols != K)
            throw std::runtime_error("Alpha.n_cols != K");
        if (Y_block.n_rows != Alpha_block.n_rows)
            throw std::runtime_error("Y.n_rows != Alpha.n_rows (N mismatch)");

        conditions.push_back(cond_name);
        Y.push_back(Y_block);
        Alpha.push_back(Alpha_block);
    }

    // Convenience: return an R list (if you return to R)
    Rcpp::List as_list() const
    {
        Rcpp::List Y_out(Y.size());
        Rcpp::List A_out(Alpha.size());

        for (size_t i = 0; i < Y.size(); ++i)
        {
            // n_student x n_item integer matrix
            // n_student × n_skill logical matrix
            Y_out[i] = wrap_imat(Y[i]);
            A_out[i] = wrap_uchar_as_logical(Alpha[i]);
        }

        Y_out.attr("names") = conditions;
        A_out.attr("names") = conditions;

        return Rcpp::List::create(Rcpp::Named("condition") = conditions,
                                  Rcpp::Named("Y") = Y_out,     // list of N×J
                                  Rcpp::Named("alpha") = A_out, // list of N×K
                                  Rcpp::Named("J") = static_cast<int>(J),
                                  Rcpp::Named("K") = static_cast<int>(K));
    }
};

struct SimulationResults
{
    std::vector<double> reaction_times;
    std::vector<unsigned int> responses;
    std::vector<std::optional<std::string>> conditions;

    // Constructor with pre-allocation
    SimulationResults(size_t total_trials)
    {
        reaction_times.reserve(total_trials);
        responses.reserve(total_trials);
        conditions.reserve(total_trials);
    }

    // Add method for efficient storage
    void add_trial(unsigned int response, double rt,
                   const std::optional<std::string> &condition = std::nullopt)
    {
        responses.push_back(response);
        reaction_times.push_back(rt);
        conditions.push_back(condition);
    }
};

std::tuple<std::string, std::string, bool>
parse_cell_name(const std::string &cell_name)
{
    // Special case: no factors, generic cell
    if (cell_name == "Cell")
    {
        return {"Cell", "", true};
        // or {"", "", true} depending on how you want to treat it
    }

    const size_t last_dot_pos = cell_name.find_last_of('.');

    if (last_dot_pos == std::string::npos)
    {
        Rcpp::warning("Irregular cell name (no dot found): " + cell_name);
        return {"", "", false};
    }

    return {cell_name.substr(0, last_dot_pos),
            cell_name.substr(last_dot_pos + 1), true};
}

///////////////////////////////////////////////////
/* ------------- Likelihood  ------------- */
/* TODO: merge it with other new_likelihood*/
///////////////////////////////////////////////////
std::shared_ptr<likelihood::likelihood_class>
new_likelihood_for_simulation(const Rcpp::S4 &rt_model_r)
{
    // auto d_ptr = new_design_light_rt_model(rt_model_r);
    auto d_ptr = new_design_light(rt_model_r);
    Rcpp::S4 model_r = rt_model_r.slot("model");
    std::string model_str = get_model_type(model_r);

    auto is_positive_drift =
        Rcpp::as<std::vector<bool>>(rt_model_r.slot("is_positive_drift"));

    return std::make_shared<likelihood::likelihood_class>(d_ptr, model_str,
                                                          is_positive_drift);
}

///////////////////////////////////////////////////
/* Simulation output  ------------- */
///////////////////////////////////////////////////
Rcpp::DataFrame new_DataFrame(const SimulationResults &results)
{
    size_t n = results.reaction_times.size();
    Rcpp::CharacterVector cond(n);

    for (size_t i = 0; i < n; ++i)
    {
        if (results.conditions[i].has_value())
        {
            cond[i] = results.conditions[i].value();
        }
        else
        {
            cond[i] =
                Rcpp::CharacterVector::get_na(); // NA if optional is empty
        }
    }

    return Rcpp::DataFrame::create(Rcpp::Named("RT") = results.reaction_times,
                                   Rcpp::Named("R") = results.responses,
                                   Rcpp::Named("Condition") = cond,
                                   Rcpp::Named("stringsAsFactors") = false);
}
