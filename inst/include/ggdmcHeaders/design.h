#pragma once

#include "model_utils.h"
#include "type_aliases.h"

namespace design
{

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
    {
        // This function inferred two member variables:
        // 1. m_parameter_names,
        // 2. m_n_parameter.
        m_is_free_parameter.resize(m_n_parameter_x_condition, false);

        for (size_t i = 0; i < m_parameter_x_condition_names.size(); ++i)
        {
            // Going over all parameters and check if it is not listed
            // in the key position of the constant_names.
            std::string str = m_parameter_x_condition_names[i];

            if (m_constants.find(str) == m_constants.end())
            {
                // Add to the result if not found
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
        for (const auto &item : m_constants)
        {
            m_constant_names.push_back(item.first);
        }

        for (const auto &str : m_parameter_x_condition_names)
        {
            auto it = m_constants.find(str);
            if (it != m_constants.end())
            {
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

    /* ----------------non-accumulator models------------------------- */
    bool m_has_accumulators = true;

    // strVec m_accumulator_names;
    // strVec m_cell_names;
    strVec m_factor_names;
    // strVec m_factor_names_hint; // optional: if your add_M needs something
    // strVec m_parameter_x_condition_names;
    // MapStrDbl m_constants;

    // ---------- shared init ----------
    void init_common(const MapStrVec &parameter_map, const MapStrVec &factors)
    {
        m_n_core_parameter = parameter_map.size();
        find_core_parameters(parameter_map);

        m_parameter_x_condition_names = add_M(parameter_map, factors);
        m_n_parameter_x_condition = m_parameter_x_condition_names.size();
    }

    // ---------- no-acc path ----------
    void init_no_accumulators(const MapStrVec &parameter_map,
                              const MapStrVec &factors)
    {
        auto [cells, facs] =
            build_cell_names_no_accumulator(parameter_map, factors);
        m_cell_names = std::move(cells);
        m_factor_names = std::move(facs);

        m_n_cell = m_cell_names.size();
        m_n_accumulator = 1;

        m_model_boolean =
            build_model_boolean_noaccu_nomatch(parameter_map, factors);
    }
    // ---------- finish-up common tail ----------
    void finalize()
    {
        find_free_parameters();
        find_constant_parameters();

        /*---------- parameter map to be used in likelihood_class----------*/
        m_tmp_param_map.resize(m_n_accumulator);
        m_param_map.resize(m_n_accumulator);

        // allocate_parameters();
        // prepare_parameter_matrix();
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

    // 3.
    std::map<std::string, double> m_constants;
    std::vector<std::string> m_constant_names;
    std::vector<double> m_constant_values;

    std::vector<std::vector<std::vector<bool>>> m_model_boolean;
    std::vector<std::vector<unsigned int>> m_node_1_index;

    // Created members
    std::vector<std::string> m_free_parameter_names; // no fixed parameters
    size_t m_n_free_parameter;

    /* ---------------- Members and functions -------------------------- */
    std::vector<std::vector<uint2D>> m_tmp_param_map;
    std::vector<std::vector<uint2D>> m_param_map;
    std::vector<bool> m_is_free_parameter;
    std::vector<double2D> m_parameter_matrix;

    arma::ucube m_model_boolean_arma;
    arma::cube m_parameter_matrix_arma;

    std::string m_model_str;

    // Work with model_utils; TODO: make this constructor obsoleted
    design_class(
        std::map<std::string, std::map<std::string, std::string>> match_map,
        std::map<std::string, std::vector<std::string>> parameter_map,
        std::map<std::string, std::vector<std::string>> factors,
        std::vector<std::string> accumulator_names,
        std::map<std::string, double> constants)
    {
        Rcpp::Rcout << "Constructor WITH accumulators\n";

        m_n_core_parameter = parameter_map.size();
        m_accumulator_names = accumulator_names;
        m_constants = constants;

        /*-------------- Inferred members --------------*/
        find_core_parameters(parameter_map);

        // 1. (1) tmp_cell_names = std::vector<std::string>
        //    (2) tmp_factor_names = std::vector<std::string>
        auto [tmp_cell_names, tmp_factor_names] =
            build_cell_names(parameter_map, factors, accumulator_names);
        m_cell_names = tmp_cell_names;
        // m_factor_names = tmp_factor_names;

        // 2. strVec add_M(const std::map<std::string, strVec> &parameter_map,
        // const std::map<std::string, strVec> &factors);
        // m_parameter_x_condition_names = strVec
        m_parameter_x_condition_names = add_M(parameter_map, factors);

        m_n_cell = m_cell_names.size();
        m_n_parameter_x_condition = m_parameter_x_condition_names.size();
        m_n_accumulator = m_accumulator_names.size();

        // 3. using uint2D = m_node_1_index
        m_node_1_index = get_node_1_index(m_cell_names, m_accumulator_names);

        // 4. uint3D = m_model_boolean
        m_model_boolean = build_model_boolean(parameter_map, factors,
                                              accumulator_names, match_map);

        // Debugging
        m_model_boolean_arma = build_model_boolean_arma(
            parameter_map, factors, accumulator_names, match_map);

        find_free_parameters();
        find_constant_parameters();

        /*---------- parameter map to be used in likelihood_class----------*/
        m_tmp_param_map.resize(m_n_accumulator);
        m_param_map.resize(m_n_accumulator);

        allocate_parameters();
        transform();
        prepare_parameter_matrix();
    }

    // Temp constructor (same as the one in design_light.h) to
    // bridge design.h to merge with design_light.h
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
            m_n_accumulator = 1;
        }

        m_n_cell = m_cell_names.size();
        m_n_parameter_x_condition = m_parameter_x_condition_names.size();

        find_free_parameters();

        if (m_constants.size() != 0)
        {
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

    // ---- Constructor WITHOUT accumulators ----; retired this constructor,
    // too.
    design_class(const MapStrVec &parameter_map, const MapStrVec &factors,
                 const MapStrDbl &constants)
        : m_has_accumulators(false), m_constants(constants)
    {
        // Rcpp::Rcout << "Constructor WITHOUT accumulators\n";
        // m_constants = constants;

        init_common(parameter_map, factors);
        init_no_accumulators(parameter_map, factors);
        finalize();
    }
    /* ------------------ Public functions ------------------ */
    void allocate_parameters()
    {
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

    void prepare_parameter_matrix()
    {
        // Resize the parameter matrix to match the dimensions
        m_parameter_matrix.resize(m_n_cell);
        m_parameter_matrix_arma = arma::cube(
            m_n_cell, m_n_core_parameter, m_n_accumulator, arma::fill::zeros);

        for (size_t cell_idx = 0; cell_idx < m_n_cell; ++cell_idx)
        {
            m_parameter_matrix[cell_idx].resize(m_n_core_parameter);

            for (size_t para_idx = 0; para_idx < m_n_core_parameter; ++para_idx)
            {
                m_parameter_matrix[cell_idx][para_idx].resize(m_n_accumulator);

                for (size_t accu_idx = 0; accu_idx < m_n_accumulator;
                     ++accu_idx)
                {
                    // As a defult, set it to NaN
                    m_parameter_matrix[cell_idx][para_idx][accu_idx] =
                        std::numeric_limits<double>::quiet_NaN();

                    unsigned int index =
                        m_param_map[accu_idx][cell_idx][para_idx][0];
                    unsigned int is_free_param =
                        m_param_map[accu_idx][cell_idx][para_idx][1];

                    // The user set it to a fixed value, so we place the user
                    // indicated value into the slot.
                    if (!is_free_param)
                    {
                        m_parameter_matrix[cell_idx][para_idx][accu_idx] =
                            m_constant_values[index];

                        m_parameter_matrix_arma(cell_idx, para_idx, accu_idx) =
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
                    m_parameter_matrix_arma(cell_idx, para_idx, accu_idx) =
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
                    m_parameter_matrix_arma(cell_idx, para_idx, accu_idx) +=
                        m_parameter_matrix[cell_idx][para_idx - 1][accu_idx];
                }
            }
        }
    }

    /* ---------------Model specific methods--------------- */
    void transform()
    {
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
                    // 1st index used node_idx
                    m_param_map[accu_idx][cell_idx] =
                        std::move(m_tmp_param_map[node_idx][cell_idx]);
                }
                else
                {
                    // 1st index used acc_idx
                    m_param_map[accu_idx][cell_idx] =
                        std::move(m_tmp_param_map[accu_idx][cell_idx]);
                }
            }
        }
    }

    /* ---------------Print methods--------------- */
    void print_parameter_map(const std::string str = "")
    {
        Rcpp::Rcout << str;
        for (size_t j = 0; j < m_n_cell; ++j) // Iterate over cells
        {
            Rcpp::Rcout << "Cell, " << m_cell_names[j] << ":\n";

            for (size_t i = 0; i < m_n_accumulator;
                 ++i) // Iterate over accumulators
            {
                uint2D &param_matrix = m_param_map[i][j];
                size_t rows = param_matrix.size();

                // Print accumulator header
                Rcpp::Rcout << "Acc " << i << ": ";

                if (rows == 0)
                {
                    Rcpp::Rcout << "(empty)\n";
                    continue;
                }
                // Print first row (indices)
                for (size_t row = 0; row < rows; ++row)
                {
                    Rcpp::Rcout << param_matrix[row][0] << " ";
                }
                Rcpp::Rcout << "\n       "; // Indentation for second row

                // Print second row (mapped values)
                for (size_t row = 0; row < rows; ++row)
                {
                    Rcpp::Rcout << param_matrix[row][1] << " ";
                }

                Rcpp::Rcout << "\n"; // End line for accumulator
            }
            Rcpp::Rcout
                << std::endl; // Separate different cells with a blank line
        }
    }
    void print_all_parameters(const std::string str = "")
    {
        Rcpp::Rcout << str;
        for (const auto &item : m_parameter_x_condition_names)
        {
            Rcpp::Rcout << item << "\t";
        }
        Rcpp::Rcout << std::endl;
    }
    void print_core_parameters(const std::string str = "")
    {
        Rcpp::Rcout << str;
        for (const auto &item : m_core_parameter_names)
        {
            Rcpp::Rcout << item << "\t";
        }
        Rcpp::Rcout << std::endl;
    }
    void print_free_parameters(const std::string str = "")
    {
        Rcpp::Rcout << str;
        for (const auto &item : m_free_parameter_names)
        {
            Rcpp::Rcout << item << "\t";
        }
        Rcpp::Rcout << std::endl;
    }
    // void print_constants(const std::string str = "")
    // {
    //     if (m_constant_names.size() != m_constant_values.size())
    //     {
    //         Rcpp::Rcout << "Constant vector length " <<
    //         m_constant_values.size()
    //                     << "\n";
    //         throw std::runtime_error("Constant vector (names and values) "
    //                                  "mismatched: name vector length: " +
    //                                  std::to_string(m_constant_names.size()));
    //     }
    //     Rcpp::Rcout << str;
    //     for (size_t i = 0; i < m_constant_names.size(); ++i)
    //     {
    //         Rcpp::Rcout << m_constant_names[i] << ": " <<
    //         m_constant_values[i]
    //                     << "\t";
    //     }
    //     Rcpp::Rcout << std::endl;
    // }

    void print_constants(const std::string &str = "")
    {
        if (m_constant_names.size() != m_constant_values.size())
        {
            Rcpp::Rcout << "❌ Constant vector length mismatch detected!\n"
                        << "  Names size: " << m_constant_names.size()
                        << ", Values size: " << m_constant_values.size()
                        << "\n";

            Rcpp::Rcout << "  --- Names ---\n";
            for (size_t i = 0; i < m_constant_names.size(); ++i)
            {
                Rcpp::Rcout << "  [" << i << "] " << m_constant_names[i]
                            << "\n";
            }

            Rcpp::Rcout << "  --- Values ---\n";
            for (size_t i = 0; i < m_constant_values.size(); ++i)
            {
                Rcpp::Rcout << "  [" << i << "] " << m_constant_values[i]
                            << "\n";
            }

            throw std::runtime_error("Constant vector mismatch: names and "
                                     "values have different lengths.");
        }

        Rcpp::Rcout << str;
        for (size_t i = 0; i < m_constant_names.size(); ++i)
        {
            Rcpp::Rcout << m_constant_names[i] << ": " << m_constant_values[i]
                        << "\t";
        }
        Rcpp::Rcout << std::endl;
    }

    void print_parameter_matrix(const std::string str = "")
    {
        Rcpp::Rcout << "\n" << str << "\n";
        for (size_t cell_idx = 0; cell_idx < m_n_cell; ++cell_idx)
        {
            Rcpp::Rcout << "[Cell, " << m_cell_names[cell_idx]
                        << "]:" << std::endl;

            for (size_t accu_idx = 0; accu_idx < m_n_accumulator; ++accu_idx)
            {
                Rcpp::Rcout << m_accumulator_names[accu_idx] << ": ";

                for (size_t para_idx = 0; para_idx < m_n_core_parameter;
                     ++para_idx)
                {
                    Rcpp::Rcout
                        << m_parameter_matrix[cell_idx][para_idx][accu_idx]
                        << "\t";
                }
                Rcpp::Rcout << std::endl; // Move to the next line after
                                          // printing all parameters
            }
            Rcpp::Rcout << std::endl; // Blank line between cells
        }
    }
};

} // namespace design
