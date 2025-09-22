#pragma once
#include <RcppArmadillo.h>
#include <algorithm>
#include <limits>
#include <map>
#include <memory>
#include <numeric> // iota
#include <optional>
#include <sstream> // std::stringstream and getline?
#include <stdexcept>
#include <vector>

using bool3D = std::vector<std::vector<std::vector<bool>>>;
using uint2D = std::vector<std::vector<unsigned int>>;
using double2D = std::vector<std::vector<double>>;
using strVec = std::vector<std::string>;
using uint1D = std::vector<unsigned int>;

namespace design
{
// alias for readability
// using strVec = std::vector<std::string>;
// using MapStrVec = std::map<std::string, strVec>;
// using MapStrStr = std::map<std::string, std::map<std::string, std::string>>;
// using MapStrDbl = std::map<std::string, double>;

class design_class
{
  private:
    void find_core_parameters(
        std::map<std::string, std::vector<std::string>> parameter_map)
    {
        m_core_parameter_names.reserve(m_n_core_parameter);
        for (const auto &param : parameter_map)
        {
            m_core_parameter_names.push_back(param.first);
        }
    }
    void find_free_parameters()
    { // This function inferred two member variables:
        m_is_free_parameter.resize(m_n_parameter_x_condition, false);

        for (size_t i = 0; i < m_parameter_x_condition_names.size(); ++i)
        {
            // Going over all parameters and check if it is not listed
            // in the key position of the constant_names.
            const std::string &str = m_parameter_x_condition_names[i];

            // If NOT found in m_constants => it's a free parameter
            if (m_constants.find(str) == m_constants.end())
            {
                m_free_parameter_names.push_back(str);
                m_is_free_parameter[i] = true;
            }
        }

        m_n_free_parameter = m_free_parameter_names.size();
    }
    void find_constant_parameters()
    {
        // This function inferred two member variables:
        // 1. m_constant_names
        // 2. m_constant_value

        m_constant_names.clear();
        m_constant_values.clear();
        // Push names and values in the SAME order to keep indices aligned
        m_constant_names.reserve(m_constants.size());
        m_constant_values.reserve(m_constants.size());

        for (const auto &str : m_parameter_x_condition_names)
        {
            auto it = m_constants.find(str);
            if (it != m_constants.end())
            {
                m_constant_names.push_back(it->first);
                m_constant_values.push_back(it->second);
            }
        }
    }

    size_t find_free_parameter_index(const std::string &target_name)
    {
        auto it = std::find(m_free_parameter_names.begin(),
                            m_free_parameter_names.end(), target_name);
        if (it != m_free_parameter_names.end())
        {
            return std::distance(m_free_parameter_names.begin(), it);
        }
        throw std::runtime_error(
            "Parameter name not found in the parameter vector: " + target_name);
    }
    size_t find_constant_index(const std::string &target_name)
    {
        auto it = std::find(m_constant_names.begin(), m_constant_names.end(),
                            target_name);
        if (it != m_constant_names.end())
        {
            return std::distance(m_constant_names.begin(), it);
        }
        throw std::runtime_error("Parameter name not found (in constants): " +
                                 target_name);
    }
    void map_parameter2condition(size_t cell_idx, size_t accu_idx)
    {
        size_t core_param_i = 0;

        for (size_t param_idx = 0; param_idx < m_n_parameter_x_condition;
             param_idx++)
        {
            if (!m_model_boolean[cell_idx][param_idx][accu_idx])
                continue;

            const std::string &target_name =
                m_parameter_x_condition_names[param_idx];

            bool is_free = m_is_free_parameter[param_idx];

            size_t idx = is_free ? find_free_parameter_index(target_name)
                                 : find_constant_index(target_name);

            m_tmp_param_map[accu_idx][cell_idx][core_param_i][0] = idx;
            m_tmp_param_map[accu_idx][cell_idx][core_param_i][1] = is_free;
            core_param_i++;
        }
    }

  public:
    // 1. process parameter_map
    // m_n_parameter_x_condition includes fixed parameters
    size_t m_n_core_parameter;
    std::vector<std::string> m_core_parameter_names;

    // 2. copy over accumulator names, cell_names and
    // parameter_x_condition_names
    std::vector<std::string> m_accumulator_names;
    std::vector<std::string> m_cell_names;
    std::vector<std::string> m_parameter_x_condition_names;
    size_t m_n_accumulator, m_n_cell, m_n_parameter_x_condition;

    // 3. prepare constants for the parameter map
    std::map<std::string, double> m_constants;
    std::vector<std::string> m_constant_names;
    std::vector<double> m_constant_values;

    std::vector<std::vector<std::vector<bool>>> m_model_boolean;
    std::vector<std::vector<unsigned int>> m_node_1_index;

    // Created members
    std::vector<std::string> m_free_parameter_names; // no fixed parameters
    size_t m_n_free_parameter;

    std::string m_model_str;

    // Constructor
    design_class(std::map<std::string, std::vector<std::string>> &parameter_map,
                 std::vector<std::string> &accumulator_names,
                 std::vector<std::string> &cell_names,
                 std::vector<std::string> &parameter_x_condition_names,
                 std::map<std::string, double> &constants,
                 std::vector<std::vector<std::vector<bool>>> &model_boolean,
                 std::string model_str = "")
        : m_accumulator_names(std::move(accumulator_names)),
          m_cell_names(std::move(cell_names)),
          m_parameter_x_condition_names(std::move(parameter_x_condition_names)),
          m_constants(std::move(constants)),
          m_model_boolean(std::move(model_boolean)),
          m_model_str(std::move(model_str))
    {
        // 1.
        m_n_core_parameter = parameter_map.size();
        find_core_parameters(parameter_map);

        // 2.
        m_n_accumulator = m_accumulator_names.size();
        if (m_n_accumulator == 0)
        {
            // Rcpp::Rcout
            //     << "Non-accumulator model, setting n_accumulator to 1\n";
            m_n_accumulator = 1;
        }

        m_n_cell = m_cell_names.size();
        m_n_parameter_x_condition = m_parameter_x_condition_names.size();

        find_free_parameters();

        if (m_constants.size() != 0)
        {
            // Rcpp::Rcout << "m_constant size = " << m_constants.size() <<
            // "\n";
            find_constant_parameters();
        }

        // 3. We should know what model type the user wants.
        if (m_model_str == "lba")
        {
            m_node_1_index =
                get_node_1_index(m_cell_names, m_accumulator_names);
        }

        if (m_model_str == "lba" || m_model_str == "fastdm" ||
            m_model_str == "cdm")
        {
            m_tmp_param_map.resize(m_n_accumulator);
            m_param_map.resize(m_n_accumulator);
            allocate_parameters();
            transform();
            prepare_parameter_matrix();
        }
        else if (m_model_str == "hyper")
        {
            // Rcpp::Rcout << "hyper model type " << m_model_str
            //             << " does not prepare parameter matrix;\n";
        }
        else
        {
            Rcpp::Rcout << "model type " << m_model_str << "\n";

            throw std::runtime_error(
                "Unknown model type detected in design_light.h");
        }
    }

    ~design_class()
    {
    }

    /* ---------------- Members and functions -------------------------- */
    // Data structures to store the pre-computed mapping
    std::vector<std::vector<uint2D>> m_tmp_param_map;
    std::vector<std::vector<uint2D>> m_param_map;
    std::vector<bool> m_is_free_parameter;

    /* ----- Prepare the info for the parameter matrix -----*/
    // The outer layer of the std::vector is a n_cell vector.
    // The middle layer of the std::vector is a n_core_parameter vector
    // The inner most layer of the std::vector is a n_accumulator vector.
    std::vector<double2D> m_parameter_matrix;

    void prepare_parameter_matrix()
    {
        // Rcpp::Rcout << "prepare parameter matrix\n\n";
        // Resize the parameter matrix to match the dimensions
        m_parameter_matrix.resize(m_n_cell);

        for (size_t cell_idx = 0; cell_idx < m_n_cell; ++cell_idx)
        {
            m_parameter_matrix[cell_idx].resize(m_n_core_parameter);

            for (size_t para_idx = 0; para_idx < m_n_core_parameter; ++para_idx)
            {
                m_parameter_matrix[cell_idx][para_idx].resize(m_n_accumulator);
                for (size_t accu_idx = 0; accu_idx < m_n_accumulator;
                     ++accu_idx)
                {
                    m_parameter_matrix[cell_idx][para_idx][accu_idx] =
                        std::numeric_limits<double>::quiet_NaN();
                }

                for (size_t accu_idx = 0; accu_idx < m_n_accumulator;
                     ++accu_idx)
                {
                    unsigned int index =
                        m_param_map[accu_idx][cell_idx][para_idx][0];
                    unsigned int is_free_param =
                        m_param_map[accu_idx][cell_idx][para_idx][1];

                    if (!is_free_param)
                    {
                        m_parameter_matrix[cell_idx][para_idx][accu_idx] =
                            m_constant_values[index];
                    }
                }
            }
        }
    }
    void set_parameter_values(size_t cell_idx,
                              const std::vector<double> &parameters)
    {
        for (size_t para_idx = 0; para_idx < m_n_core_parameter; ++para_idx)
        {
            for (size_t accu_idx = 0; accu_idx < m_n_accumulator; ++accu_idx)
            {
                // The 2nd column indicates if this is a free or a fixed
                // parameter.
                if (m_param_map[accu_idx][cell_idx][para_idx][1])
                {
                    // The 1st column stores the index on the core parameter
                    // (NA) vector, where the estimated parameter value should
                    // be placed.
                    m_parameter_matrix[cell_idx][para_idx][accu_idx] =
                        parameters[m_param_map[accu_idx][cell_idx][para_idx]
                                              [0]];
                }

                /* LBA specific b = A + B; A, b, mean_v, sd_v, st0, t0 */
                // parameter_matrix is a 6 x 2 matrix. This part should be
                // refactored to a separate function.
                if (m_core_parameter_names[para_idx] == "B")
                {

                    m_parameter_matrix[cell_idx][para_idx][accu_idx] +=
                        m_parameter_matrix[cell_idx][para_idx - 1][accu_idx];
                }
            }
        }
    }
    void allocate_parameters()
    {
        // Rcpp::Rcout << "allocate parameters: " << m_model_str << "\n";

        for (size_t accu_idx = 0; accu_idx < m_n_accumulator; ++accu_idx)
        {
            m_tmp_param_map[accu_idx].resize(m_n_cell);

            for (size_t cell_idx = 0; cell_idx < m_n_cell; ++cell_idx)
            {
                m_tmp_param_map[accu_idx][cell_idx].resize(m_n_core_parameter);
                for (size_t para_idx = 0; para_idx < m_n_core_parameter;
                     ++para_idx)
                {
                    m_tmp_param_map[accu_idx][cell_idx][para_idx].resize(2);
                }

                map_parameter2condition(cell_idx, accu_idx);
            }
        }
    }
    void transform()
    {
        // Rcpp::Rcout << "transform: " << m_model_str << "\n";

        for (size_t accu_idx = 0; accu_idx < m_n_accumulator; accu_idx++)
        {
            m_param_map[accu_idx].resize(m_n_cell);
            for (size_t cell_idx = 0; cell_idx < m_n_cell; cell_idx++)
            {
                m_param_map[accu_idx][cell_idx].resize(m_n_core_parameter);

                for (size_t para_idx = 0; para_idx < m_n_core_parameter;
                     ++para_idx)
                {
                    m_param_map[accu_idx][cell_idx][para_idx].resize(2);
                }

                if (m_model_str == "lba")
                {
                    size_t node_idx = m_node_1_index[cell_idx][accu_idx];
                    m_param_map[accu_idx][cell_idx] =
                        std::move(m_tmp_param_map[node_idx][cell_idx]);
                }
                else
                {
                    m_param_map[accu_idx][cell_idx] =
                        std::move(m_tmp_param_map[accu_idx][cell_idx]);
                }
            }
        }
    }
};

} // namespace design
