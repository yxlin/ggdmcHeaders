#pragma once

#include "tnorm.h" // <- must include here for ggdmcPrior to work
#include <cstddef>
#include <iomanip> // For formatting output
#include <memory>

namespace prior
{
class prior_class
{
  private:
    void validate_input_sizes() const
    {
        if (m_stdp1.size() != m_nparameter || m_lower.size() != m_nparameter ||
            m_upper.size() != m_nparameter ||
            m_dist_code.size() != m_nparameter ||
            m_log_p.size() != m_nparameter)
        {
            throw std::runtime_error("All arguments to the prior constructor "
                                     "must have the same length.");
        }
    }

    double dcauchy_trunc(double x, double location, double scale, double lower,
                         double upper, bool log_p = false)
    {
        // Check for invalid scale or bounds
        if (scale <= 0 || lower >= upper)
        {
            return R_NaN; // Return NaN for invalid inputs
        }

        // Standard Cauchy density
        double cauchy_density = Rf_dcauchy(x, location, scale, 0);

        // CDF at bounds
        double cdf_lower = Rf_pcauchy(lower, location, scale, 1, 0);
        double cdf_upper = Rf_pcauchy(upper, location, scale, 1, 0);

        // Normalized truncated density
        double trunc_density = 0.0;
        if (x >= lower && x <= upper)
        {
            trunc_density = cauchy_density / (cdf_upper - cdf_lower);
        }

        // Return log-density if log_p = 1
        if (log_p)
        {
            return (trunc_density > 0.0) ? log(trunc_density) : R_NegInf;
        }
        else
        {
            return trunc_density;
        }
    }

    double rcauchy(double location, double scale)
    {
        return location +
               scale * Rf_rcauchy(0.0, 1.0); // R's rcauchy() uses standard
                                             // Cauchy (location=0, scale=1)
    }

    void set_bounds(double &bound, double default_value)
    {
        bound = std::isnan(bound) ? default_value : bound;
    }
    void handleTnorm(size_t i, size_t &tnorm_counter)
    {
        set_bounds(m_lower[i], R_NegInf);
        set_bounds(m_upper[i], R_PosInf);

        m_tnorm_objects[tnorm_counter].set_parameters(
            m_p0[i], m_p1[i], m_lower[i], m_upper[i], m_log_p[i]);

        tnorm_counter++;
    }
    void handleBetaLU(size_t i)
    {
        set_bounds(m_lower[i], 0);
        set_bounds(m_upper[i], 1);

        m_beta_range = m_upper[i] - m_lower[i];
        m_log_beta_range = std::log(m_upper[i] - m_lower[i]);
    }
    void handleGammaLOrLnormL(size_t i)
    {
        set_bounds(m_lower[i], 0);
    }
    void handleCauchy(size_t i)
    {
        set_bounds(m_lower[i], R_NegInf);
        set_bounds(m_upper[i], R_PosInf);
    }
    void handleUniform(size_t i)
    {
        set_bounds(m_lower[i], R_NegInf);
        set_bounds(m_upper[i], R_PosInf);

        if (m_p0[i] == m_p1[i])
        {
            throw std::runtime_error(
                "lower_bound == upper_bound in the uniform PDF.");
        }
    }
    void handleNorm(size_t i)
    {
    }
    void handleDefault()
    {
        Rcpp::Rcout
            << "You have entered an undefined distribution. Please submit "
               "your request to the "
               "developer.\n";
    }

    void initialise_distributions()
    {
        size_t tnorm_counter = 0;

        for (size_t i = 0; i < m_nparameter; i++)
        {
            switch (m_dist_code[i])
            {
            case TNORM:
                handleTnorm(i, tnorm_counter);
                break;
            case BETA_LU:
                handleBetaLU(i);
                break;
            case GAMMA_L:
            case LNORM_L:
                handleGammaLOrLnormL(i);
                break;
            case CAUCHY:
                handleCauchy(i);
                break;
            case UNIF:
                handleUniform(i);
                break;
            case NORM:
                handleNorm(i);
                break;
            default:
                handleDefault();
                break;
            }
        }
    }
    void initialise_tnorm_objects()
    {
        for (size_t i = 0; i < m_nparameter; ++i)
        {
            if (m_dist_code[i] == 1) // dist type code = 1 mean tnorm
            {
                m_tnorm_index.push_back(i);
            }
        }
        m_tnorm_objects.resize(m_tnorm_index.size());

        for (size_t i = 0; i < m_tnorm_index.size(); ++i)
        {

            size_t j = m_tnorm_index[i];
            set_bounds(m_lower[j], R_NegInf);
            set_bounds(m_upper[j], R_PosInf);
            m_tnorm_objects[i] = tnorm::tnorm_class(
                m_stdp0[j], m_stdp1[j], m_lower[j], m_upper[j], m_log_p[j]);
        }
    }

  public:
    std::vector<double> m_stdp0, m_stdp1, m_lower, m_upper;
    std::vector<unsigned int> m_dist_code;
    std::vector<bool> m_log_p;
    unsigned int m_nparameter;

    arma::vec m_p0, m_p1;
    std::vector<unsigned int> m_tnorm_index;
    std::vector<tnorm::tnorm_class> m_tnorm_objects;

    double m_beta_range, m_log_beta_range;

    enum DistributionType
    {
        TNORM = 1,   // OK
        BETA_LU = 2, // OK
        GAMMA_L = 3, // OK
        LNORM_L = 4, // OK
        CAUCHY = 5,  // OK
        UNIF = 6,    // Error
        NORM = 7,
        DEFAULT // DEFAULT will implicitly be 0
    };

    DistributionType resolve_option(int type)
    {
        if (type < 1 || type > 6)
            return DEFAULT;

        switch (type)
        {
        case 1:
            return TNORM;
        case 2:
            return BETA_LU;
        case 3:
            return GAMMA_L;
        case 4:
            return LNORM_L;
        case 5:
            return CAUCHY;
        case 6:
            return UNIF;
        case 7:
            return NORM;
        }

        return DEFAULT; // Should never reach here
    }

    prior_class(std::vector<double> p0, std::vector<double> p1,
                std::vector<double> lower, std::vector<double> upper,
                std::vector<unsigned int> dist_code, std::vector<bool> log_p)
        : m_stdp0(std::move(p0)), m_stdp1(std::move(p1)),
          m_lower(std::move(lower)), m_upper(std::move(upper)),
          m_dist_code(std::move(dist_code)), m_log_p(std::move(log_p)),
          m_nparameter(m_stdp0.size())
    {
        validate_input_sizes();
        arma::vec tmp_p0(m_stdp0.data(), m_stdp0.size(), true);
        arma::vec tmp_p1(m_stdp1.data(), m_stdp1.size(), true);
        m_p0 = tmp_p0;
        m_p1 = tmp_p1;

        initialise_tnorm_objects();
        initialise_distributions();
    }

    ~prior_class()
    {
    }

    // The main functions
    void dprior(std::vector<double> &parameters, std::vector<double> &out)
    {
        size_t tnorm_counter = 0;

        // Lambda to compute x for shifted distributions
        auto compute_x = [&](size_t i)
        {
            return std::isfinite(m_lower[i]) ? parameters[i] - m_lower[i]
                                             : parameters[i];
        };

        for (size_t i = 0; i < m_nparameter; i++)
        {
            switch (m_dist_code[i])
            {
            case TNORM:
            {

                m_tnorm_objects[tnorm_counter].set_parameters(m_stdp0[i],
                                                              m_stdp1[i]);
                out[i] = m_tnorm_objects[tnorm_counter].d(parameters[i]);

                tnorm_counter++; // go to the next free parameter, if there is
                                 // any.
                break;
            }
            case BETA_LU:
            {
                double x = (parameters[i] - m_lower[i]) / m_beta_range;
                // check if Rf_dbeta is faster
                double density = R_NegInf;

                if (m_stdp0[i] >= 0 && m_stdp1[i] >= 0)
                {
                    density = Rf_dbeta(x, m_p0[i], m_p1[i], m_log_p[i]);
                }

                out[i] = m_log_p[i] ? density - m_log_beta_range
                                    : density / m_beta_range;
                break;
            }
            case GAMMA_L:
            case LNORM_L:
            {
                // R::dlnorm may be faster
                out[i] = (m_dist_code[i] == GAMMA_L)
                             ? Rf_dgamma(compute_x(i), m_stdp0[i], m_stdp1[i],
                                         m_log_p[i])
                             : Rf_dlnorm(compute_x(i), m_stdp0[i], m_stdp1[i],
                                         m_log_p[i]);
                break;
            }
            case CAUCHY:
            {
                // Not tested
                out[i] = dcauchy_trunc(parameters[i], m_stdp0[i], m_stdp1[i],
                                       m_lower[i], m_upper[i], m_log_p[i]);
                break;
            }
            case UNIF:
            {
                // Boost dunif ran faster than R::dunif

                double val =
                    Rf_dunif(parameters[i], m_stdp0[i], m_stdp1[i], m_log_p[i]);
                out[i] = std::isnan(val) ? -1e10 : val;
                break;
            }
            case NORM:
            {
                out[i] = Rf_dnorm4(parameters[i], m_stdp0[i], m_stdp1[i],
                                   m_log_p[i]);
                break;
            }

            default:
            {
                Rcpp::Rcout << "An undefined distribution (dprior).\n";
                out[i] = NA_REAL;
                break;
            }
            }
        }
    }
    void dprior(arma::vec &parameters, arma::vec &out)
    {
        size_t tnorm_counter = 0;

        // Lambda to compute x for shifted distributions
        auto compute_x = [&](size_t i)
        {
            return std::isfinite(m_lower[i]) ? parameters[i] - m_lower[i]
                                             : parameters[i];
        };

        for (size_t i = 0; i < m_nparameter; i++)
        {
            switch (m_dist_code[i])
            {
            case TNORM:
            {

                m_tnorm_objects[tnorm_counter].set_parameters(m_p0[i], m_p1[i]);
                out[i] = m_tnorm_objects[tnorm_counter].d(parameters[i]);

                tnorm_counter++; // go to the next free parameter, if there is
                                 // any.
                break;
            }
            case BETA_LU:
            {
                double x = (parameters[i] - m_lower[i]) / m_beta_range;

                double density = R_NegInf;
                if (m_p0[i] >= 0 && m_p1[i] >= 0)
                {
                    density = Rf_dbeta(x, m_p0[i], m_p1[i], m_log_p[i]);
                }

                out[i] = m_log_p[i] ? density - m_log_beta_range
                                    : density / m_beta_range;
                break;
            }
            case GAMMA_L:
            case LNORM_L:
            {
                // R::dlnorm may be faster
                out[i] =
                    (m_dist_code[i] == GAMMA_L)
                        ? Rf_dgamma(compute_x(i), m_p0[i], m_p1[i], m_log_p[i])
                        : Rf_dlnorm(compute_x(i), m_p0[i], m_p1[i], m_log_p[i]);
                break;
            }
            case CAUCHY:
            {
                // Not tested
                out[i] = dcauchy_trunc(parameters[i], m_p0[i], m_p1[i],
                                       m_lower[i], m_upper[i], m_log_p[i]);
                break;
            }
            case UNIF:
            {
                // Boost dunif ran faster than R::dunif
                // out[i] = stats::dunif(parameters[i], m_p0[i], m_p1[i],
                // m_log_p[i]);
                double val =
                    Rf_dunif(parameters[i], m_p0[i], m_p1[i], m_log_p[i]);
                out[i] = std::isnan(val) ? -1e10 : val;
                break;
            }
            case NORM:
            {
                out[i] = Rf_dnorm4(parameters[i], m_p0[i], m_p1[i], m_log_p[i]);
                break;
            }

            default:
            {
                Rcpp::Rcout << "An undefined distribution (dprior).\n";
                out[i] = NA_REAL;
                break;
            }
            }
        }
    }

    arma::vec rprior()
    {
        arma::vec out(m_nparameter);

        size_t tnorm_counter = 0;

        for (size_t i = 0; i < m_nparameter; i++)
        {
            switch (m_dist_code[i])
            {
            case TNORM:
            {
                out[i] = m_tnorm_objects[tnorm_counter++].r();
                break;
            }
            case BETA_LU:
            {
                out[i] = m_lower[i] + Rf_rbeta(m_stdp0[i], m_stdp1[i]) *
                                          (m_upper[i] - m_lower[i]);
                break;
            }
            case GAMMA_L:
            case LNORM_L:
            {
                out[i] = (m_dist_code[i] == GAMMA_L)
                             ? Rf_rgamma(m_stdp0[i], m_stdp1[i]) + m_lower[i]
                             : Rf_rlnorm(m_stdp0[i], m_stdp1[i]) + m_lower[i];

                break;
            }
            case CAUCHY:
            {
                Rcpp::Rcout << "The rcauchy is untested.\n";
                out[i] = rcauchy(m_stdp0[i], m_stdp1[i]);
                break;
            }
            case UNIF:
            {
                out[i] = Rf_runif(m_stdp0[i], m_stdp1[i]);
                break;
            }
            case NORM:
            {
                out[i] = Rf_rnorm(m_stdp0[i], m_stdp1[i]);
                break;
            }
            default:
            {
                Rcpp::Rcout << "An undefined distribution.\n";
                out[i] = NA_REAL;
                break;
            }
            }
        }
        return out;
    }
    double sumlogprior(arma::vec parameters)
    {

        arma::vec out(m_nparameter);
        dprior(parameters, out);
        return arma::accu(out);
    }
    // Use only in the print function
    std::vector<std::string> convert_dist_code2string() const
    {
        std::vector<std::string> dist_str(m_nparameter);
        for (size_t i = 0; i < m_nparameter; ++i)
        {
            if (m_dist_code[i] == 1)
            {
                dist_str[i] = "tnorm";
            }
            else if (m_dist_code[i] == 2)
            {
                dist_str[i] = "beta_lu";
            }
            else if (m_dist_code[i] == 3)
            {
                dist_str[i] = "gamma_l";
            }
            else if (m_dist_code[i] == 4)
            {
                dist_str[i] = "lnorm_l";
            }
            else if (m_dist_code[i] == 5)
            {
                dist_str[i] = "cauchy";
            }
            else if (m_dist_code[i] == 6)
            {
                dist_str[i] = "uniform";
            }
            else if (m_dist_code[i] == 7)
            {
                dist_str[i] = "norm";
            }
            else
            {
                dist_str[i] = "default (unknown)";
            }
        }
        return dist_str;
    }
    void print(std::vector<std::string> parameter_names) const
    {
        std::vector<std::string> dist_str = convert_dist_code2string();

        // Find the maximum length of parameter names to adjust the column width
        size_t max_name_length = 0;
        for (const auto &name : parameter_names)
        {
            if (name.length() > max_name_length)
            {
                max_name_length = name.length();
            }
        }

        // Ensure the "Index" column is wide enough to fit the longest name
        size_t index_column_width = std::max(
            max_name_length, size_t(5)); // At least 5 characters for "Index"

        Rcpp::Rcout << std::fixed << std::setprecision(2);
        Rcpp::Rcout
            << "-------------------------------------------------------------"
               "-------\n";
        Rcpp::Rcout
            << "| Parameter |   p0   |   p1   | Lower  | Upper  | Log_P | "
               "Dist  |\n";
        Rcpp::Rcout
            << "-------------------------------------------------------------"
               "-------\n";

        for (size_t i = 0; i < m_nparameter; ++i)
        {
            Rcpp::Rcout << "| " << std::setw(index_column_width)
                        << parameter_names[i] << " | " << std::setw(6)
                        << m_stdp0[i] << " | " << std::setw(6) << m_stdp1[i]
                        << " | " << std::setw(6) << m_lower[i] << " | "
                        << std::setw(6) << m_upper[i] << " | " << std::setw(5)
                        << (m_log_p[i] ? "True" : "False") << " | "
                        << std::setw(4) << dist_str[i] << " |\n";
        }
        Rcpp::Rcout
            << "-------------------------------------------------------------"
               "-------\n";
    }
};
///////////////////////////////////////////////////
/* Prior ---------------------------------------*/
///////////////////////////////////////////////////
inline std::shared_ptr<prior_class> new_prior(const Rcpp::List &p_prior_r)
{
    std::vector<std::string> parameter_names = p_prior_r.names();
    unsigned int nparameter = p_prior_r.size();

    std::vector<double> p0(nparameter), p1(nparameter);
    std::vector<double> lower(nparameter), upper(nparameter);
    std::vector<bool> log_p(nparameter);
    std::vector<unsigned int> dist_type(nparameter);

    for (size_t i = 0; i < nparameter; i++)
    {
        if (!p_prior_r.containsElementNamed(parameter_names[i].c_str()))
        {
            Rcpp::stop("Parameter name not found in prior list.");
        }
        // Use the name to extract the individual prior distribution, so
        // we must order the name first.
        Rcpp::List parameter = p_prior_r[parameter_names[i]];

        if (parameter.size() < 6)
        {
            Rcpp::stop("Parameter '" + parameter_names[i] + "' has " +
                       std::to_string(parameter.size()) +
                       " element. It should have 6 elements.");
        }

        p0[i] = Rcpp::as<double>(parameter[0]);
        p1[i] = Rcpp::as<double>(parameter[1]);
        lower[i] = Rcpp::as<double>(parameter[2]);
        upper[i] = Rcpp::as<double>(parameter[3]);
        dist_type[i] = Rcpp::as<unsigned int>(parameter[4]);
        log_p[i] = Rcpp::as<bool>(parameter[5]);
    }

    return std::make_shared<prior_class>(p0, p1, lower, upper, dist_type,
                                         log_p);
}
} // namespace prior

using PriorPtr = std::shared_ptr<prior::prior_class>;
