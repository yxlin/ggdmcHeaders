#pragma once
#include "design_light.h"
#include "lba.h"
#include "simulation_type_casting.h"
#include <string>
#include <vector>
void simulate_each_condition(
    const std::shared_ptr<design::design_class> &design,
    lba::lba_class &lba_obj, const std::vector<double> &parameters,
    const std::vector<bool> &is_positive_drift, unsigned int n_trial_per_cell,
    bool use_inverse_method, SimulationResults &results, bool debug)
{
    std::string prev_condition;
    std::string node_1 = design->m_accumulator_names[0];
    if (debug)
    {
        Rcpp::Rcout << "Node 1 = " << node_1 << std::endl;
    }

    for (size_t cell_idx = 0; cell_idx < design->m_n_cell; ++cell_idx)
    {
        const auto &cell_name = design->m_cell_names[cell_idx];
        const auto [condition, accumulator, valid] = parse_cell_name(cell_name);

        if (accumulator != node_1)
        {
            if (debug)
            {
                Rcpp::Rcout << "Verify if this is node 1...\n";
                Rcpp::Rcout << "Current and previous conditions: " << condition
                            << ", " << prev_condition << std::endl;
                Rcpp::Rcout << "Current node is: " << accumulator
                            << ", which is not the node 1, " << node_1 << "\n"
                            << std::endl;
            }
            continue;
        }
        if (!valid || (cell_idx > 0 && condition == prev_condition))
        {
            if (debug)
            {
                Rcpp::Rcout << "Verify if this condition has done before:"
                            << std::endl;
                Rcpp::Rcout << "Current and previous conditions:" << condition
                            << ", " << prev_condition << std::endl
                            << std::endl;
            }

            continue;
        }

        prev_condition = condition;

        design->set_parameter_values(cell_idx, parameters);

        lba_obj.set_parameters(design->m_parameter_matrix[cell_idx],
                               is_positive_drift);

        bool is_valid = lba_obj.validate_parameters();
        if (!is_valid)
        {
            lba_obj.print_parameters("simulate_lba_trials");
            throw std::runtime_error("Invalid paraemters");
        }
        if (debug)
        {
            Rcpp::Rcout << "LBA model implements node 1 (current condition = "
                        << cell_name << ") to simulate all "
                        << design->m_n_accumulator << " possible responses\n";
            lba_obj.print_parameters();
        }

        std::vector<std::pair<unsigned int, double>> cell_data(
            n_trial_per_cell);
        if (use_inverse_method)
        {
            lba_obj.r_inverse(cell_data);
        }
        else
        {
            lba_obj.r(cell_data);
        }

        for (const auto &trial : cell_data)
        {
            results.add_trial(trial.first, trial.second, condition);
        }
    }
}
