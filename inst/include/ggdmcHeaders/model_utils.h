#pragma once

#include <RcppArmadillo.h> // for Rcout only
#include <algorithm>       // sort
// #include <armadillo>
#include <iomanip>
// #include <iostream>
#include <map>
// #include <string>
#include <utility>
// #include <vector>

#include <sstream> // std::stringstream and getline?
#include <tuple>
#include <unordered_map>

using bool3D = std::vector<std::vector<std::vector<bool>>>;
using bool2D = std::vector<std::vector<bool>>; // redefined
using bool1D = std::vector<bool>;

using uint2D = std::vector<std::vector<unsigned int>>; // redefined
using uint1D = std::vector<unsigned int>;

using double2D = std::vector<std::vector<double>>;

using strVec = std::vector<std::string>;
using strMap = std::map<std::string, std::string>;

///////////////////////////////////////////////////////////////
/* ----------split_parameter_condition (para_list) ----------*/
///////////////////////////////////////////////////////////////
inline strVec split(const std::string &input, char delimiter)
{
    strVec parts;
    std::stringstream ss(input);
    std::string part;
    while (std::getline(ss, part, delimiter))
    {
        parts.push_back(part);
    }
    return parts;
}

// Flexible split function
inline std::tuple<std::string, std::string, std::string, std::string,
                  std::string, std::string, std::string, std::string>
flexible_split(const std::string &param)
{
    strVec parts = split(param, '.');

    size_t num_parts = parts.size();

    if (num_parts == 1)
    {
        return make_tuple(parts[0], "", "", "", "", "", "", "");
    }
    else if (num_parts == 2)
    {
        return make_tuple(parts[0], parts[1], "", "", "", "", "", "");
    }
    else if (num_parts == 3)
    {
        return make_tuple(parts[0], parts[1], parts[2], "", "", "", "", "");
    }
    else if (num_parts == 4)
    {
        return make_tuple(parts[0], parts[1], parts[2], parts[3], "", "", "",
                          "");
    }
    else if (num_parts == 5)
    {
        return make_tuple(parts[0], parts[1], parts[2], parts[3], parts[4], "",
                          "", "");
    }
    else if (num_parts == 6)
    {
        return make_tuple(parts[0], parts[1], parts[2], parts[3], parts[4],
                          parts[5], "", "");
    }
    else if (num_parts == 7)
    {
        return make_tuple(parts[0], parts[1], parts[2], parts[3], parts[4],
                          parts[5], parts[6], "");
    }
    else if (num_parts == 8)
    {
        return make_tuple(parts[0], parts[1], parts[2], parts[3], parts[4],
                          parts[5], parts[6], parts[7]);
    }
    else
    { // Handle more than 8 parts by combining the rest
        Rcpp::Rcout << "The complex factor structure has exceeded the limit."
                    << std::endl;
        std::string combined;
        for (size_t i = 7; i < parts.size(); ++i)
        {
            combined += parts[i];
            if (i != parts.size() - 1)
            {
                combined += ".";
            }
        }
        return make_tuple(parts[0], parts[1], parts[2], parts[3], parts[4],
                          parts[5], parts[6], combined);
    }
}

inline std::vector<strVec> split_parameter_condition(const strVec &parameters)
{
    std::vector<strVec> para_level;

    for (const auto &param : parameters)
    {
        auto para_level_tuple = flexible_split(param);
        strVec non_empty_parts;

        // Check each part of the tuple and add non-empty parts to the vector
        if (!std::get<0>(para_level_tuple).empty())
            non_empty_parts.push_back(std::get<0>(para_level_tuple));

        if (!std::get<1>(para_level_tuple).empty())
            non_empty_parts.push_back(std::get<1>(para_level_tuple));

        if (!std::get<2>(para_level_tuple).empty())
            non_empty_parts.push_back(std::get<2>(para_level_tuple));

        if (!std::get<3>(para_level_tuple).empty())
            non_empty_parts.push_back(std::get<3>(para_level_tuple));

        if (!std::get<4>(para_level_tuple).empty())
            non_empty_parts.push_back(std::get<4>(para_level_tuple));

        if (!std::get<5>(para_level_tuple).empty())
            non_empty_parts.push_back(std::get<5>(para_level_tuple));

        if (!std::get<6>(para_level_tuple).empty())
            non_empty_parts.push_back(std::get<6>(para_level_tuple));

        if (!std::get<7>(para_level_tuple).empty())
            non_empty_parts.push_back(std::get<7>(para_level_tuple));

        para_level.push_back(non_empty_parts);
    }

    return para_level;
}

///////////////////////////////////////////////////////////////
/* ---------- add_M functions ----------*/
///////////////////////////////////////////////////////////////
// Create the parameter_x_condition names
inline strVec add_M(const std::map<std::string, strVec> &parameter_map,
                    const std::map<std::string, strVec> &factors)
{
    strVec all_parameters;

    for (const auto &[param, associations] : parameter_map)
    {
        // Sort associations: "S" first, then alphabetical order
        strVec sorted_associations(associations.begin(), associations.end());

        std::sort(sorted_associations.begin(), sorted_associations.end(),
                  [](const std::string &a, const std::string &b)
                  {
                      if (a == "S")
                          return true; // "S" comes first
                      if (b == "S")
                          return false; // "S" comes first
                      return a < b;     // Otherwise, sort alphabetically
                  });

        bool has_M = (find(associations.begin(), associations.end(), "M") !=
                      associations.end());

        // Collect all factor levels for this parameter
        std::vector<strVec> factor_levels;
        // associations is "S", "G" or "D", "H", "M", or others
        for (const auto &sorted_association : sorted_associations)
        {
            if (factors.find(sorted_association) != factors.end())
            {
                factor_levels.push_back(factors.at(sorted_association));
            }
        }

        // Sort factor levels: "S" first, then alphabetical order
        // auto sort_order = [](const vector<string> &a, const vector<string>
        // &b)
        // {
        //     bool a_has_S = find(a.begin(), a.end(), "S") != a.end();
        //     bool b_has_S = find(b.begin(), b.end(), "S") != b.end();
        //     if (a_has_S != b_has_S)
        //         return a_has_S;
        //     return a < b;
        // };
        // sort(factor_levels.begin(), factor_levels.end(), sort_order);

        if (factor_levels.empty())
        {
            // No factors, just handle "M" if present
            if (has_M)
            {
                all_parameters.push_back(param + ".true");
                all_parameters.push_back(param + ".false");
            }
            else
            {
                all_parameters.push_back(param);
            }
        }
        else
        {
            // Generate combinations of factor levels
            strVec combined_levels = {""}; // Start with an empty level
            for (const auto &levels : factor_levels)
            {
                strVec new_combined_levels;
                for (const auto &level : levels)
                {
                    for (const auto &combined : combined_levels)
                    {
                        new_combined_levels.push_back(
                            combined + (combined.empty() ? "" : ".") + level);
                    }
                }
                combined_levels = new_combined_levels;
            }

            // Append ".true" and ".false" if "M" is present
            for (const auto &level : combined_levels)
            {
                if (has_M)
                {
                    all_parameters.push_back(param + "." + level + ".true");
                    all_parameters.push_back(param + "." + level + ".false");
                }
                else
                {
                    all_parameters.push_back(param + "." + level);
                }
            }
        }
    }

    std::sort(all_parameters.begin(), all_parameters.end());

    return all_parameters;
}

inline strVec
umap_add_M(const std::unordered_map<std::string, strVec> &parameter_map,
           const std::unordered_map<std::string, strVec> &factors)
{
    strVec out;

    for (const auto &[param, associations] : parameter_map)
    {
        // Sort associations: "S" first, then alphabetical order
        strVec sorted_associations(associations.begin(), associations.end());

        std::sort(sorted_associations.begin(), sorted_associations.end(),
                  [](const std::string &a, const std::string &b)
                  {
                      if (a == "S")
                          return true; // "S" comes first
                      if (b == "S")
                          return false; // "S" comes first
                      return a < b;     // Otherwise, sort alphabetically
                  });

        bool has_M = (std::find(associations.begin(), associations.end(),
                                "M") != associations.end());

        // Collect all factor levels for this parameter
        std::vector<strVec> factor_levels;
        // associations is "S", "G" or "D", "H", "M", or others
        for (const auto &sorted_association : sorted_associations)
        {
            if (factors.find(sorted_association) != factors.end())
            {
                factor_levels.push_back(factors.at(sorted_association));
            }
        }

        if (factor_levels.empty())
        {
            // No factors, just handle "M" if present
            if (has_M)
            {
                out.push_back(param + ".true");
                out.push_back(param + ".false");
            }
            else
            {
                out.push_back(param);
            }
        }
        else
        {
            // Generate combinations of factor levels
            strVec combined_levels = {""}; // Start with an empty level
            for (const auto &levels : factor_levels)
            {
                strVec new_combined_levels;
                for (const auto &level : levels)
                {
                    for (const auto &combined : combined_levels)
                    {
                        new_combined_levels.push_back(
                            combined + (combined.empty() ? "" : ".") + level);
                    }
                }
                combined_levels = new_combined_levels;
            }

            // Append ".true" and ".false" if "M" is present
            for (const auto &level : combined_levels)
            {
                if (has_M)
                {
                    out.push_back(param + "." + level + ".true");
                    out.push_back(param + "." + level + ".false");
                }
                else
                {
                    out.push_back(param + "." + level);
                }
            }
        }
    }

    return out;
}

///////////////////////////////////////////////////////////////
/* ---------- build_cell_names ----------*/
///////////////////////////////////////////////////////////////
inline void validate_SR(const std::map<std::string, strVec> &factors,
                        const strVec &accumulators)
{
    // Ensure the default factor "S" is present in the factors
    if (factors.find("S") == factors.end())
    {
        throw std::runtime_error("The 'factors' argument must include the 'S' "
                                 "factor (model_utils.cpp.)");
    }

    // Ensure the levels of "S" match the length of the accumulator
    if (factors.at("S").size() != accumulators.size())
    {
        throw std::runtime_error("The number of levels for 'S' must match the "
                                 "number of accumulators (model_utils.cpp.)");
    }
}

// Helper function to compute the Cartesian product of vectors
inline std::vector<strVec> cartesian_product(const std::vector<strVec> &input)
{
    std::vector<strVec> result = {{}};
    for (const auto &vec : input)
    {
        std::vector<strVec> temp;
        for (const auto &res : result)
        {
            for (const auto &item : vec)
            {
                temp.push_back(res);
                temp.back().push_back(item);
            }
        }
        result = temp;
    }
    return result;
}

inline std::pair<strVec, strVec>
build_cell_names(const std::map<std::string, strVec> &parameter_map,
                 const std::map<std::string, strVec> &factors,
                 const strVec &accumulators)
{
    // map<string, vector<string>> parameter_map0 = {
    //     {"A", {"1"}}, {"B", {"1"}}, {"t0", {"1"}}, {"mean_v",
    //     {"POLITICAL_VIEW", "M"}},
    //     {"sd_v", {"1"}}, {"st0", {"1"}}};
    // map<string, vector<string>> factors0 = {
    //     {"S", {"s1", "s2"}}, {"POLITICAL_VIEW", {"liberal", "conservative"}}
    //     };
    // vector<string> accumulators0 = {"r1", "r2"};

    validate_SR(factors, accumulators);

    // Extract the factors and their levels from the parameterMap
    std::map<std::string, strVec> factorLevels;

    // Although the for loop follows the order the user enter in the
    // parameter_map, e.g., std::map will should do the alphabetical ordering.
    for (const auto &[param, factorList] : parameter_map)
    {
        // factorList could be c("D", "M") in the mean_v = c("D", "M")
        // for instance.
        for (const auto &factor : factorList)
        {
            if (factors.find(factor) != factors.end())
            {
                // factors.at(factor) is "sti_1", "sti_2", "sti_3", "sti_4"
                // in (e.g.):
                //   factors <- list(
                //   S = c("sti_1", "sti_2", "sti_3", "sti_4"),
                //   D = c("d1", "d2"),
                //   E = c("e1", "e2"))
                // , so this creates a dictionary.
                factorLevels[factor] = factors.at(factor);
            }
        }
    }

    // Always include the default factor "S" in the combinations
    if (factorLevels.find("S") == factorLevels.end())
    {
        factorLevels["S"] = factors.at("S");
    }

    // Ensure the default factor 'S' is always first
    strVec sortedFactors = {"S"};
    for (const auto &[factor, _] : factorLevels)
    {
        if (factor != "S")
        {
            sortedFactors.push_back(factor);
        }
    }
    sort(sortedFactors.begin() + 1, sortedFactors.end());

    // Generate all possible combinations of factor levels
    std::vector<strVec> factorCombinations;
    std::vector<strVec> levels;
    for (const auto &factor : sortedFactors)
    {
        levels.push_back(factorLevels.at(factor));
    }
    factorCombinations = cartesian_product(levels);

    // Generate all possible combinations of factor levels and accumulator
    strVec combinations;
    for (const auto &factorCombo : factorCombinations)
    {
        for (const auto &accumulator : accumulators)
        {
            std::string combo;
            for (const auto &level : factorCombo)
            {
                combo += level + ".";
            }
            combo += accumulator;
            combinations.push_back(combo);
        }
    }

    // Sort combinations alphabetically
    std::sort(combinations.begin(), combinations.end());

    return {combinations, sortedFactors};
}

///////////////////////////////////////////////////////////////
/* ---------- parameter_x_condition.hpp ----------*/
///////////////////////////////////////////////////////////////
/* ----------parameter_x_condition association ----------*/
inline std::vector<bool> is_core_parameter_x_condition(
    const std::map<std::string, strVec> &parameter_map,
    const std::map<std::string, strVec> &factors)
{
    std::vector<bool> out(parameter_map.size(), false);
    int i = 0;

    for (const auto &param : parameter_map)
    {
        for (const auto &factor : factors)
        {
            if (std::find(param.second.begin(), param.second.end(),
                          factor.first) != param.second.end())
            {
                // Ensure out[i] is set to true if at least one factor matches
                out[i] = true;
                break; // Exit early if any match is found
            }
        }
        i++;
    }

    return out;
}

inline std::vector<bool> is_parameter_condition_associated(
    const std::map<std::string, strVec> &parameter_map,
    const strVec &parameter_x_condition,
    const std::map<std::string, strVec> &factors)
{
    // Step 1: Get the core parameter associations
    std::vector<bool> is_asso =
        is_core_parameter_x_condition(parameter_map, factors);

    // Step 2: Initialize the output vector with false values
    std::vector<bool> out(parameter_x_condition.size(), false);

    // Step 3: Iterate through each parameter in parameter_x_condition
    for (size_t i = 0; i < parameter_x_condition.size(); ++i)
    {
        const std::string &param = parameter_x_condition[i];

        // Step 4: Iterate through each core parameter in parameter_map
        size_t j = 0;
        for (const auto &core_param_pair : parameter_map)
        {
            const std::string &core_param = core_param_pair.first;

            // Step 5: Check if the core parameter is associated and is part of
            // the current parameter
            // Rcpp::Rcout << "Current param and core_param = " << param << ", "
            //             << core_param << ", j = " << j << ", " << is_asso[j]
            //             << "\n";
            if (is_asso[j])
            {
                if (param == core_param)
                {
                    out[i] = true;
                    break;
                }
                else if (param.find(core_param + ".") == 0)
                {
                    // Rcpp::Rcout << param.find(core_param + ".") << "\t";
                    out[i] = true;
                    break;
                }
            }
            j++;
        }
        // Rcpp::Rcout << "\n";
    }
    // Rcpp::Rcout << "\n\n";

    return out;
}

/* ---------- build_model_boolean tools ----------*/
// Function to extract the first component from a string split by a delimiter
inline std::string get_stimulus_level(const std::string &input,
                                      char delimiter = '.')
{
    strVec components;
    std::stringstream ss(input);
    std::string item;

    while (std::getline(ss, item, delimiter))
    {
        components.push_back(item);
    }

    // Return the first component if the vector is not empty
    if (!components.empty())
    {
        return components[0];
    }

    // Return an empty string if no components were found
    return "";
}

// Function to create a map of factor cells
inline strMap get_factor_cells(const std::string &cell_name,
                               const strVec &factor_names)
{
    strMap out;
    strVec factor_levels = split(cell_name, '.');

    // The last string after "." in a cell name is always response
    // code (e.g., "resp_1"), so the -1 is to remove it.
    // Ensure we don't exceed the size of the components vector
    if (factor_names.size() != (factor_levels.size() - 1))
    {
        throw std::runtime_error(
            "The number of factors must be the same as the "
            "number of cell name component - 1");
    }

    for (size_t i = 0; i < factor_names.size(); ++i)
    {
        // The last component in the factor level is the resposne code.
        // Because the cell_name include a level of all the factors, this
        // cover all factor.
        out[factor_names[i]] = factor_levels[i];
    }
    // The map will make the S comes after
    return out;
}

inline bool is_at_the_same_level_with_M(const strVec &factor_keys,
                                        const strMap &factor_cells,
                                        const strVec &parameter_and_levels,
                                        size_t para_idx,
                                        const std::string &separator = "")
{
    if (factor_keys.size() + 1 != parameter_and_levels.size())
    { // The factor key assocated with a core parameter.
        throw std::runtime_error(
            "The factor_keys size must be one shorter than "
            "the parameter_and_levels size. Check "
            "the cell names, the parameter_x_condition "
            "namaes or any typos in your factor names in "
            "'parameter_map' and 'factors'? (M cases).");
    }

    std::string cell_name_level_combo, para_cell_level_combo;
    for (size_t i = 0; i < factor_keys.size(); ++i)
    {
        if (factor_keys[i] == "M")
            continue;
        // Because the S key always comes first, even the map (factor_cell)
        // put S at the end, its factor level will come first still.
        cell_name_level_combo += factor_cells.at(factor_keys[i]) + separator;

        if (parameter_and_levels[i + 1] != "true" &&
            parameter_and_levels[i + 1] != "false")
        {
            para_cell_level_combo += parameter_and_levels[i + 1] + separator;
        }
    }
    // First element in the parameter_and_levels is the core_parameter name, so
    // I have to skip it. And the last element in the parameter_and_levels is
    // true or false if the factor key is "M", so I have to skip it too.

    return cell_name_level_combo == para_cell_level_combo;
}

inline bool is_at_the_same_level_no_M(const strVec &factor_keys,
                                      const strMap &factor_cells,
                                      const strVec &parameter_and_levels,
                                      size_t para_idx,
                                      const std::string &separator = "")
{
    if (factor_keys.size() + 1 != parameter_and_levels.size())
    {
        throw std::runtime_error(
            "The factor_keys size must be one shorter than "
            "the parameter_and_levels size. Check "
            "the cell names, the parameter_x_condition "
            "namaes or any typos in your factor names in "
            "'parameter_map' and 'factors'? ( no M cases).");
    }

    std::string cell_name_level_combo, para_cell_level_combo;

    for (size_t i = 0; i < factor_keys.size(); ++i)
    {
        cell_name_level_combo += factor_cells.at(factor_keys[i]) + separator;
        para_cell_level_combo += parameter_and_levels[i + 1] + separator;
    }
    // First element in the parameter_and_levels is the core_parameter name, so
    // I have to skip it. And the last element in the parameter_and_levels is
    // true or false if the factor key is "M", so I have to skip it too.
    return cell_name_level_combo == para_cell_level_combo;
}

inline bool is_this_accumulator(const std::map<std::string, strMap> &match_map,
                                const std::string &accumulator,
                                const std::string &stimulus_cell,
                                const std::string &param_type,
                                const std::string &factor_key = "M")
{
    if (match_map.find(factor_key) == match_map.end())
    {
        throw std::runtime_error("Factor key '" + factor_key +
                                 "' not found in match_map.");
    }

    auto it = match_map.at(factor_key).find(stimulus_cell);
    if (it == match_map.at(factor_key).end())
    {
        throw std::runtime_error("Stimulus cell '" + stimulus_cell +
                                 "' not found for factor key '" + factor_key +
                                 "' in the match_map.");
    }

    std::string target_accumulator = it->second;
    return (param_type == "true" && accumulator == target_accumulator) ||
           (param_type == "false" && accumulator != target_accumulator);
}

inline void
handle_non_asso_parameter(bool3D &model_boolean, const std::string &parameter,
                          const std::string &accumulator,
                          const std::string &stimulus_level,
                          const std::map<std::string, strMap> &match_map,
                          size_t cell_idx, size_t para_idx, size_t accu_idx)
{
    if (parameter.find(".true") != std::string::npos ||
        parameter.find(".false") != std::string::npos)
    {
        std::string param_type = parameter.substr(parameter.rfind(".") + 1);

        bool is_the_right_acc = is_this_accumulator(match_map, accumulator,
                                                    stimulus_level, param_type);
        if (is_the_right_acc)
        {
            model_boolean[cell_idx][para_idx][accu_idx] = true;
        }
    }
    else
    {
        model_boolean[cell_idx][para_idx][accu_idx] = true;
    }
}

inline void handle_non_asso_parameter(
    arma::ucube &model_boolean, const std::string &parameter,
    const std::string &accumulator, const std::string &stimulus_level,
    const std::map<std::string, strMap> &match_map, size_t cell_idx,
    size_t para_idx, size_t accu_idx)
{
    if (parameter.find(".true") != std::string::npos ||
        parameter.find(".false") != std::string::npos)
    {
        std::string param_type = parameter.substr(parameter.rfind(".") + 1);

        bool is_the_right_acc = is_this_accumulator(match_map, accumulator,
                                                    stimulus_level, param_type);
        if (is_the_right_acc)
        {
            model_boolean(cell_idx, para_idx, accu_idx) = true;
        }
    }
    else
    {
        model_boolean(cell_idx, para_idx, accu_idx) = true;
    }
}

// Function to handle the case when "M" is found in factor_keys
inline void handle_m_factor_case(
    bool3D &model_boolean, const strVec &factor_keys,
    const strMap &factor_cells, const strVec &parameter_and_levels,
    const std::string &parameter,
    const std::map<std::string, strMap> &match_map,
    const std::string &accumulator, const std::string &stimulus_level,
    size_t cell_idx, size_t para_idx, size_t accu_idx)
{
    if (factor_keys.size() >= 2)
    {
        // factor_cell still under the map sorting control. S will come last
        // sometimes. RESOLVED.
        bool is_same_level = is_at_the_same_level_with_M(
            factor_keys, factor_cells, parameter_and_levels, para_idx);

        std::string is_true_or_false =
            parameter.substr(parameter.rfind(".") + 1);
        bool is_right_acc = is_this_accumulator(
            match_map, accumulator, stimulus_level, is_true_or_false);

        if (is_same_level && is_right_acc)
        {
            model_boolean[cell_idx][para_idx][accu_idx] = true;
        }
    }
}

inline void handle_m_factor_case(
    arma::ucube &model_boolean, const strVec &factor_keys,
    const strMap &factor_cells, const strVec &parameter_and_levels,
    const std::string &parameter,
    const std::map<std::string, strMap> &match_map,
    const std::string &accumulator, const std::string &stimulus_level,
    size_t cell_idx, size_t para_idx, size_t accu_idx)
{
    if (factor_keys.size() >= 2)
    {
        // factor_cell still under the map sorting control. S will come last
        // sometimes. RESOLVED.
        bool is_same_level = is_at_the_same_level_with_M(
            factor_keys, factor_cells, parameter_and_levels, para_idx);

        std::string is_true_or_false =
            parameter.substr(parameter.rfind(".") + 1);
        bool is_right_acc = is_this_accumulator(
            match_map, accumulator, stimulus_level, is_true_or_false);

        if (is_same_level && is_right_acc)
        {
            model_boolean(cell_idx, para_idx, accu_idx) = true;
        }
    }
}

// Function to handle the case when "M" is not found in factor_keys
inline void handle_no_m_factor_case(bool3D &model_boolean,
                                    const strVec &factor_keys,
                                    const strMap &factor_cells,
                                    const strVec &parameter_and_levels,
                                    size_t cell_idx, size_t para_idx,
                                    size_t accu_idx)
{
    if (factor_keys.size() == 1)
    {
        std::string which_level = factor_cells.at(factor_keys[0]);
        std::string para_x_cell = parameter_and_levels[1];

        if (which_level == para_x_cell)
        {
            model_boolean[cell_idx][para_idx][accu_idx] = true;
        }
    }
    else
    {
        bool is_same_level = is_at_the_same_level_no_M(
            factor_keys, factor_cells, parameter_and_levels, para_idx);

        if (is_same_level)
        {
            model_boolean[cell_idx][para_idx][accu_idx] = true;
        }
    }
}

inline void handle_no_m_factor_case(arma::ucube &model_boolean,
                                    const strVec &factor_keys,
                                    const strMap &factor_cells,
                                    const strVec &parameter_and_levels,
                                    size_t cell_idx, size_t para_idx,
                                    size_t accu_idx)
{
    if (factor_keys.size() == 1)
    {
        std::string which_level = factor_cells.at(factor_keys[0]);
        std::string para_x_cell = parameter_and_levels[1];

        if (which_level == para_x_cell)
        {
            model_boolean(cell_idx, para_idx, accu_idx) = true;
        }
    }
    else
    {
        bool is_same_level = is_at_the_same_level_no_M(
            factor_keys, factor_cells, parameter_and_levels, para_idx);

        if (is_same_level)
        {
            model_boolean(cell_idx, para_idx, accu_idx) = true;
        }
    }
}

inline bool3D
build_model_boolean(const std::map<std::string, strVec> &parameter_map,
                    const std::map<std::string, strVec> &factors,
                    const strVec &accumulators,
                    const std::map<std::string, strMap> &match_map)
{
    auto [cell_names, factor_names] =
        build_cell_names(parameter_map, factors, accumulators);

    auto parameter_x_condition_names = add_M(parameter_map, factors);
    auto para_list = split_parameter_condition(parameter_x_condition_names);

    size_t n_cell = cell_names.size();                       // row
    size_t n_parameter = parameter_x_condition_names.size(); // col
    size_t n_accumulator = accumulators.size();              // slice

    bool3D model_boolean(n_cell,
                         bool2D(n_parameter, bool1D(n_accumulator, false)));

    bool1D is_asso = is_parameter_condition_associated(
        parameter_map, parameter_x_condition_names, factors);

    for (size_t accu_idx = 0; accu_idx < n_accumulator; ++accu_idx)
    {
        for (size_t cell_idx = 0; cell_idx < n_cell; ++cell_idx)
        {
            std::string stimulus_level =
                get_stimulus_level(cell_names[cell_idx]);
            strMap factor_cells =
                get_factor_cells(cell_names[cell_idx], factor_names);

            for (size_t para_idx = 0; para_idx < n_parameter; ++para_idx)
            {
                // parameter_and_levels is a vector<string>
                const auto &parameter_and_levels = para_list[para_idx];
                // core_parameter is like "A", "B" or "mean_v" etc.
                std::string core_parameter = parameter_and_levels[0];

                // parameter is like "A.s1", "B.d1" or "mean_v.d1.true" etc.
                std::string parameter = parameter_x_condition_names[para_idx];
                // D, M, H or D, H, M
                strVec user_entered_factor_keys =
                    parameter_map.at(core_parameter);

                // Sorting the factor keys here, so we may allow the user to
                // enter factor key in any order.
                strVec factor_keys(user_entered_factor_keys.begin(),
                                   user_entered_factor_keys.end());

                std::sort(factor_keys.begin(), factor_keys.end(),
                          [](const std::string &a, const std::string &b)
                          {
                              if (a == "S")
                                  return true; // "S" comes first
                              if (b == "S")
                                  return false; // "S" comes first
                              if (a == "M")
                                  return false; // "M" comes last
                              if (b == "M")
                                  return true; // "M" comes last
                              return a < b;    // Otherwise, sort alphabetically
                          });

                if (is_asso[para_idx])
                {
                    // No association with the match factor
                    if (std::find(factor_keys.begin(), factor_keys.end(),
                                  "M") != factor_keys.end())
                    {
                        handle_m_factor_case(
                            model_boolean, factor_keys, factor_cells,
                            parameter_and_levels, parameter, match_map,
                            accumulators[accu_idx], stimulus_level, cell_idx,
                            para_idx, accu_idx);
                    }
                    else
                    {
                        handle_no_m_factor_case(
                            model_boolean, factor_keys, factor_cells,
                            parameter_and_levels, cell_idx, para_idx, accu_idx);
                    }
                }
                else
                {
                    handle_non_asso_parameter(model_boolean, parameter,
                                              accumulators[accu_idx],
                                              stimulus_level, match_map,
                                              cell_idx, para_idx, accu_idx);
                }
            }
        }
    }
    return model_boolean;
}

inline arma::ucube
build_model_boolean_arma(const std::map<std::string, strVec> &parameter_map,
                         const std::map<std::string, strVec> &factors,
                         const strVec &accumulators,
                         const std::map<std::string, strMap> &match_map)
{
    auto [cell_names, factor_names] =
        build_cell_names(parameter_map, factors, accumulators);

    auto parameter_x_condition_names = add_M(parameter_map, factors);
    auto para_list = split_parameter_condition(parameter_x_condition_names);

    size_t n_cell = cell_names.size();                       // row
    size_t n_parameter = parameter_x_condition_names.size(); // col
    size_t n_accumulator = accumulators.size();              // slice

    // for(const auto)
    // for (const auto &item : cell_names)

    // Rcpp::Rcout << "Inside build_model_boolean_arma\n";

    arma::ucube model_boolean =
        arma::ucube(n_cell, n_parameter, n_accumulator, arma::fill::zeros);

    // Rcpp::Rcout << "model_boolean_arma\n";
    // model_boolean.print();

    bool1D is_asso = is_parameter_condition_associated(
        parameter_map, parameter_x_condition_names, factors);

    for (size_t accu_idx = 0; accu_idx < n_accumulator; ++accu_idx)
    {
        for (size_t cell_idx = 0; cell_idx < n_cell; ++cell_idx)
        {
            std::string stimulus_level =
                get_stimulus_level(cell_names[cell_idx]);
            strMap factor_cells =
                get_factor_cells(cell_names[cell_idx], factor_names);

            for (size_t para_idx = 0; para_idx < n_parameter; ++para_idx)
            {
                // parameter_and_levels is a vector<string>
                const auto &parameter_and_levels = para_list[para_idx];
                // core_parameter is like "A", "B" or "mean_v" etc.
                std::string core_parameter = parameter_and_levels[0];

                // parameter is like "A.s1", "B.d1" or "mean_v.d1.true" etc.
                std::string parameter = parameter_x_condition_names[para_idx];
                // D, M, H or D, H, M
                strVec user_entered_factor_keys =
                    parameter_map.at(core_parameter);

                // Sorting the factor keys here, so we may allow the user to
                // enter factor key in any order.
                strVec factor_keys(user_entered_factor_keys.begin(),
                                   user_entered_factor_keys.end());

                std::sort(factor_keys.begin(), factor_keys.end(),
                          [](const std::string &a, const std::string &b)
                          {
                              if (a == "S")
                                  return true; // "S" comes first
                              if (b == "S")
                                  return false; // "S" comes first
                              if (a == "M")
                                  return false; // "M" comes last
                              if (b == "M")
                                  return true; // "M" comes last
                              return a < b;    // Otherwise, sort alphabetically
                          });

                if (is_asso[para_idx])
                {
                    // No association with the match factor
                    if (std::find(factor_keys.begin(), factor_keys.end(),
                                  "M") != factor_keys.end())
                    {
                        handle_m_factor_case(
                            model_boolean, factor_keys, factor_cells,
                            parameter_and_levels, parameter, match_map,
                            accumulators[accu_idx], stimulus_level, cell_idx,
                            para_idx, accu_idx);
                    }
                    else
                    {
                        handle_no_m_factor_case(
                            model_boolean, factor_keys, factor_cells,
                            parameter_and_levels, cell_idx, para_idx, accu_idx);
                    }
                }
                else
                {
                    handle_non_asso_parameter(model_boolean, parameter,
                                              accumulators[accu_idx],
                                              stimulus_level, match_map,
                                              cell_idx, para_idx, accu_idx);
                }
            }
        }
    }
    return model_boolean;
}

inline uint2D get_node_1_index(const strVec &cell_names,
                               const strVec &accumulators)
{
    size_t n_accumulator = accumulators.size();
    size_t n_cell = cell_names.size();

    uint2D out(n_cell, uint1D(n_accumulator, 0));

    for (size_t i = 0; i < n_cell; ++i)
    {
        strVec cell_level = split(cell_names[i], '.');
        std::string response = cell_level.back();

        auto it = std::find(accumulators.begin(), accumulators.end(), response);
        if (it == accumulators.end())
        {
            throw std::runtime_error(
                "Error: Response cell not found in accumulators");
        }
        size_t acc_index = distance(accumulators.begin(), it);

        if (acc_index == 0)
        {
            for (size_t j = 0; j < n_accumulator; ++j)
            {
                out[i][j] = j;
            }
        }
        else
        {
            std::vector<unsigned int> row(n_accumulator);
            // Fill with 0, 1, ..., naccumulator-1
            std::iota(row.begin(), row.end(), 0);
            std::swap(row[0], row[acc_index]);
            std::sort(row.begin() + 1, row.end());
            for (size_t j = 0; j < n_accumulator; ++j)
            {
                out[i][j] = row[j];
            }
        }
    }
    return out;
}