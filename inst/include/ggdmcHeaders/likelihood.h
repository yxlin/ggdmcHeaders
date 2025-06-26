#pragma once

#include "ddm.h"
#include "design_light.h"
#include "lba.h"
#include "prior.h"

namespace likelihood
{
class likelihood_class
{
  private:
    enum ModelType
    {
        LBA = 1,
        FAST_DM = 2,
        HYPER = 3,
        DEFAULT // DEFAULT will implicitly be 0
    };

    ModelType resolve_string(std::string type)
    {
        if (type == "lba")
            return LBA; // 1
        if (type == "fastdm")
            return FAST_DM; // 2
        if (type == "hyper")
            return HYPER; // 3
        return DEFAULT;   // 0
    }

    void update_empty_cells()
    {
        // It does not matter whether data_cell_names are sorted or not
        // here, but the sequence of the rt and cell_boolean from the data
        // are still critical.
        unsigned int n_cell = m_model->m_cell_names.size();
        m_is_empty_cell.resize(n_cell);

        for (size_t i = 0; i < n_cell; ++i)
        {
            m_is_empty_cell[i] = true; // Assume the cell is empty initially

            for (const auto &cell_name : m_data_cell_names)
            {
                if (m_model->m_cell_names[i] == cell_name)
                {
                    m_is_empty_cell[i] = false;
                    break;
                }
            }
        }
    }

    // Index for non_missing_rt
    void fill_in_rt_vector()
    {
        size_t rt_index = 0;
        m_rt.resize(m_model->m_n_cell);
        for (size_t i = 0; i < m_model->m_n_cell; ++i)
        {
            if (m_is_empty_cell[i])
            {
                m_rt[i] = {};
            }
            else
            {
                m_rt[i] = m_data_rt[rt_index++];
            }
        }
    }

    void lba_likelihood(const std::vector<double> &parameters, bool debug)
    {
        m_density.resize(m_model->m_n_cell);

        static lba::lba_class lba_obj;
        bool is_valid = false;

        for (size_t cell_idx = 0; cell_idx < m_model->m_n_cell; ++cell_idx)
        {
            if (m_is_empty_cell[cell_idx])
                continue;

            m_model->set_parameter_values(cell_idx, parameters);
            lba_obj.set_parameters(m_model->m_parameter_matrix[cell_idx],
                                   m_is_positive_drift);

            is_valid = lba_obj.validate_parameters(debug);

            if (debug)
            {
                const auto &cell_name = m_model->m_cell_names[cell_idx];
                lba_obj.print_parameters(cell_name);
            }

            if (is_valid)
            {
                // Rcpp::Rcout << "valid\n";
                m_density[cell_idx] = lba_obj.dlba(m_rt[cell_idx]);
            }
            else
            {
                // Rcpp::Rcout << "invalid\n";
                m_density[cell_idx].assign(m_rt[cell_idx].size(), 1e-10);
            }
        }
    }

    void print_parameter_matrix(const std::vector<std::vector<double>> &input)
    {
        // REMOVE it after debugging
        const int width = 10;
        const int precision = 3;

        for (size_t i = 0; i < input.size(); ++i)
        {
            std::vector<double> row = input[i];
            Rcpp::Rcout << "Row " << i << ": ";
            for (size_t j = 0; j < row.size(); ++j)
            {
                Rcpp::Rcout << std::setw(width) << std::fixed
                            << std::setprecision(precision) << row[j];
            }
            Rcpp::Rcout << "\n";
        }
    }

    void ddm_likelihood(const std::vector<double> &parameters, bool debug)
    {
        m_density.resize(m_model->m_n_cell);
        static ddm::ddm_class ddm_obj;
        bool is_valid = false;

        for (size_t cell_idx = 0; cell_idx < m_model->m_n_cell; ++cell_idx)
        {
            if (m_is_empty_cell[cell_idx])
                continue;

            m_model->set_parameter_values(cell_idx, parameters);
            ddm_obj.set_parameters(m_model->m_parameter_matrix[cell_idx],
                                   !m_is_positive_drift[cell_idx]);

            is_valid = ddm_obj.validate_parameters(debug);

            if (debug)
            {
                const auto &cell_name = m_model->m_cell_names[cell_idx];
                ddm_obj.print(cell_name);
            }

            if (is_valid)
            {
                m_density[cell_idx] = ddm_obj.dddm(m_rt[cell_idx]);
            }
            else
            {
                m_density[cell_idx].assign(m_rt[cell_idx].size(), 1e-10);
            }
        }
    }

  public:
    std::shared_ptr<design::design_class> m_model;

    // - m_data_rt is from the empirical data, which at some occassions will
    // have no data for some cell/condition.
    // - m_rt preserves an emply slot for the cell/condition without data.
    // Therefore, the loop going over each cell knows if a cell has no data.
    std::vector<std::vector<double>> m_data_rt, m_rt, m_density;
    std::vector<std::string> m_data_cell_names;
    std::string m_model_str;

    // is_positive_drift -> for tnorm in lba.h or is_lower in ddm.h
    std::vector<bool> m_is_empty_cell, m_is_positive_drift;

    // Hyper only members
    arma::mat m_theta_data;
    std::shared_ptr<prior::prior_class> m_p_prior;

    // LBA constructor (with is_positive_drift)
    likelihood_class(std::shared_ptr<design::design_class> model,
                     std::vector<std::vector<double>> data_rt,
                     std::vector<std::string> data_cell_names,
                     std::string model_str, std::vector<bool> is_positive_drift)
        : m_model(std::move(model)), m_data_rt(std::move(data_rt)),
          m_data_cell_names(std::move(data_cell_names)),
          m_model_str(std::move(model_str)),
          m_is_positive_drift(std::move(is_positive_drift))
    {
        update_empty_cells();
        fill_in_rt_vector();
        // Rcpp::Rcout << "likelihood constructor m_is_positive_drift\n";
    }

    // LBA simulation constructor. No data_rt, nor data_cell_names.
    // with is_positive_drift
    likelihood_class(std::shared_ptr<design::design_class> model,
                     std::string model_str, std::vector<bool> is_positive_drift)
        : m_model(std::move(model)), m_model_str(std::move(model_str)),
          m_is_positive_drift(std::move(is_positive_drift))
    {
        // Rcpp::Rcout << "Likelihood simulation constructor\n";
    }

    // Hyper constructor
    ////////////////////////////
    likelihood_class(std::shared_ptr<design::design_class> model,
                     std::shared_ptr<prior::prior_class> p_prior,
                     arma::mat theta_data, std::string model_str)
        : m_model(std::move(model)), m_p_prior(std::move(p_prior))

    {
        // Rcpp::Rcout << "Likelihood hyper constructor\n";
        m_theta_data = theta_data.t();
        m_model_str = model_str;
    }

    ~likelihood_class(){};

    // Use by the ggdmcLikelihood
    template <typename T>
    void likelihood(const T &parameters, bool debug = false)
    {
        if constexpr (std::is_same_v<T, arma::vec>)
        {
            std::vector<double> parameters_std =
                arma::conv_to<std::vector<double>>::from(parameters);
            likelihood_impl(parameters_std, debug);
        }
        else // Assumed to be std::vector<double>
        {
            likelihood_impl(parameters, debug);
        }
    }

    void likelihood_impl(const std::vector<double> &parameters,
                         bool debug = false)
    {
        // ggdmcLikelihood enters from here
        if (m_model_str == "lba")
        { // Fill in m_density with likelihoods;
            // Rcpp::Rcout << "lba_likelihood ggdmcLikelihood\n";

            lba_likelihood(parameters, debug);
        }
        else if (m_model_str == "fastdm")
        {
            ddm_likelihood(parameters, debug);
        }
        else
        {
            throw std::runtime_error("Undefined model likelihood functions");
        }
    }

    void sumloghlike(const arma::vec &parameters, double &out)
    {
        m_p_prior->m_p0 = parameters.head(parameters.n_elem / 2);
        m_p_prior->m_p1 = parameters.tail(parameters.n_elem / 2);

        // m_theta_data has been transposed in the constructor
        size_t n_subject = m_theta_data.n_cols;

        for (size_t i = 0; i < n_subject; ++i)
        {
            out += m_p_prior->sumlogprior(m_theta_data.col(i));
        }
    }

    // Use by the DE sampler
    double sumloglike(arma::vec parameters, bool debug = false)
    {
        // DE sampler enters from here
        std::vector<double> parameters_std =
            arma::conv_to<std::vector<double>>::from(parameters);

        double out = 0;
        switch (resolve_string(m_model_str))
        {
        case LBA:
            // Rcpp::Rcout << "lba_likelihood sumloglike\n";
            lba_likelihood(parameters_std, debug); // calculate m_density;
            for (const auto &inner_vec : m_density)
            {
                for (double val : inner_vec)
                {
                    out += std::log(val);
                }
            }

            break;
        case FAST_DM:
            // Rcpp::Rcout << "before ddm_likelihood\n";
            ddm_likelihood(parameters_std, debug); // calculate m_density;
            for (const auto &inner_vec : m_density)
            {
                for (double val : inner_vec)
                {
                    // DBL_MIN from <cfloat>
                    // Safeguard the DDM underflow problem that causes negative
                    // density
                    out += std::log(std::max(val, DBL_MIN));
                }
            }

            break;
        case HYPER:
            sumloghlike(parameters, out);
            break;
        case DEFAULT:
            throw std::runtime_error("Undefined model type\n");
            break;
        }
        return out;
    }

    void print_rt(const std::string &title = "Response Times (m_rt)",
                  bool print_summary = false)
    {
        // Save original formatting
        std::ios::fmtflags old_flags = Rcpp::Rcout.flags();
        std::streamsize old_precision = Rcpp::Rcout.precision();

        // Set fixed formatting with 3 decimal places
        Rcpp::Rcout << std::fixed << std::setprecision(3);

        // Print header
        Rcpp::Rcout << "\n" << title << " [" << m_rt.size() << " cells]\n";
        Rcpp::Rcout << std::string(60, '-') << "\n";

        // Print each accumulator's RTs
        for (size_t cell_idx = 0; cell_idx < m_rt.size(); ++cell_idx)
        {
            Rcpp::Rcout << "Cell " << m_model->m_cell_names[cell_idx] << ": ";

            // Print up to first 5 RTs (or all if fewer than 5)
            size_t print_count = std::min<size_t>(5, m_rt[cell_idx].size());
            for (size_t i = 0; i < print_count; ++i)
            {
                Rcpp::Rcout << std::setw(8) << m_rt[cell_idx][i];
            }

            // Add ellipsis if more RTs exist
            if (m_rt[cell_idx].size() > 5)
            {
                size_t remaining = m_rt[cell_idx].size() - 5;
                Rcpp::Rcout << " ... [+" << std::setw(2) << std::setfill('0')
                            << remaining << std::setfill(' ') << " more]";
            }

            // Summary stats
            if (!m_rt[cell_idx].empty() && print_summary)
            {
                auto [min_it, max_it] = std::minmax_element(
                    m_rt[cell_idx].begin(), m_rt[cell_idx].end());
                double mean = std::accumulate(m_rt[cell_idx].begin(),
                                              m_rt[cell_idx].end(), 0.0) /
                              m_rt[cell_idx].size();

                Rcpp::Rcout << " | Min: " << std::setw(6) << *min_it
                            << " Max: " << std::setw(6) << *max_it
                            << " Mean: " << std::setw(6) << mean;
            }

            Rcpp::Rcout << "\n";
        }

        // Restore original formatting
        Rcpp::Rcout.flags(old_flags);
        Rcpp::Rcout.precision(old_precision);
        Rcpp::Rcout << std::string(60, '-') << "\n\n";
    }
    void print_is_empty_cell(const std::string &title = "")
    {
        Rcpp::Rcout << title;
        for (const auto &item : m_is_empty_cell)
        {
            Rcpp::Rcout << item << " ";
        }
        Rcpp::Rcout << std::endl;
    }
};
} // namespace likelihood

using LPtr = std::shared_ptr<likelihood::likelihood_class>;
