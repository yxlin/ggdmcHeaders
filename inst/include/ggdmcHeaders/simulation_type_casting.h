#pragma once
#include "common_type_casting.h"
#include "likelihood.h"
#include <RcppArmadillo.h>
#include <optional>
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
/* ------------- Design helpers  ------------- */
/* TODO: merge it with other new_design*/
///////////////////////////////////////////////////
// Use in lbaModel only
std::shared_ptr<design::design_class>
new_design_light_rt_model(const Rcpp::S4 &rt_model)
{
    // 1. The user enterd elements
    Rcpp::S4 model_r = rt_model.slot("model");
    Rcpp::List parameter_map_r = model_r.slot("parameter_map");
    Rcpp::CharacterVector accumulators_r = model_r.slot("accumulators");
    Rcpp::List factors_r = model_r.slot("factors");
    Rcpp::List match_map_r = model_r.slot("match_map");
    Rcpp::List constants_r = model_r.slot("constants");
    std::string model_str = model_r.slot("type");

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

    // Light additional arguments
    auto model_boolean = R_ucube_to_std_ucube(model_r.slot("model_boolean"));

    auto cell_names =
        Rcpp::as<std::vector<std::string>>(model_r.slot("cell_names"));
    auto parameter_x_condition_names = Rcpp::as<std::vector<std::string>>(
        model_r.slot("parameter_x_condition_names"));

    auto out = std::make_shared<design::design_class>(
        parameter_map, accumulators, cell_names, parameter_x_condition_names,
        constants, model_boolean, model_str);

    return out;
}

///////////////////////////////////////////////////
/* ------------- Likelihood  ------------- */
/* TODO: merge it with other new_likelihood*/
///////////////////////////////////////////////////
std::shared_ptr<likelihood::likelihood_class>
new_likelihood_for_simulation(const Rcpp::S4 &rt_model_r)
{
    auto d_ptr = new_design_light_rt_model(rt_model_r);
    Rcpp::S4 model_r = rt_model_r.slot("model");
    std::string model_str = model_r.slot("type");

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
