#pragma once

#include "tnorm.h"
#include <iomanip> // For std::setw to align output
#include <memory>  // For std::maked_shared in lba.cpp

namespace lba
{
class lba_class
{
  private:
    unsigned int m_n_accumulator, m_n_point;
    std::vector<double> m_A, m_b, m_mean_v, m_sd_v, m_st0, m_t0;
    std::vector<double> m_denom, m_t0_actual;

    std::vector<bool> m_is_positive_drift;
    double m_min_dt, m_max_dt, m_time_step;

    void create_time_grid()
    {
        m_n_point =
            static_cast<unsigned int>((m_max_dt - m_min_dt) / m_time_step) + 1;
        m_rt_grid.resize(m_n_accumulator, std::vector<double>(m_n_point));

        for (size_t acc_idx = 0; acc_idx < m_n_accumulator; ++acc_idx)
        {
            for (size_t i = 0; i < m_n_point; ++i)
            {
                m_rt_grid[acc_idx][i] =
                    m_t0_actual[acc_idx] + m_min_dt + i * m_time_step;
            }
            m_rt_grid[acc_idx].back() = m_t0_actual[acc_idx] + m_max_dt;
        }
    }
    bool is_negative(double value, const std::string &name, bool is_debug)
    {
        if (value < 0)
        {
            if (is_debug)
                throw std::runtime_error("invalid parameter " + name + " = " +
                                         std::to_string(value));
            return false;
        }
        return true;
    }
    /* Find the closest time point index using binary search */
    double get_density(double rt, unsigned int acc_idx)
    {
        if (rt < m_min_dt || rt > m_max_dt)
            return 0.0;

        // Must generated m_rt_grid first in order to use find_time_index
        auto pdfs = theoretical_dlba();

        auto it = std::lower_bound(m_rt_grid[acc_idx].begin(),
                                   m_rt_grid[acc_idx].end(), rt);

        // Handle edge cases (RT exactly at grid boundaries)
        if (it == m_rt_grid[acc_idx].begin())
        {
            return pdfs[acc_idx][0];
        }
        if (it == m_rt_grid[acc_idx].end())
        {
            return pdfs[acc_idx].back();
        }

        // Get neighboring points
        double t_upper = *it;
        double t_lower = *(it - 1);

        // Get corresponding PDF values
        size_t idx_upper = it - m_rt_grid[acc_idx].begin();
        size_t idx_lower = idx_upper - 1;

        double pdf_upper = pdfs[acc_idx][idx_upper];
        double pdf_lower = pdfs[acc_idx][idx_lower];

        // Linear interpolation: weight = (rt - t_lower) / (t_upper - t_lower)
        double weight = (rt - t_lower) / (t_upper - t_lower);
        return pdf_lower + weight * (pdf_upper - pdf_lower);
    }

  public:
    double m_dt, m_pdf, m_cdf;
    std::vector<std::vector<double>> m_rt_grid;
    void
    set_parameters(const std::vector<std::vector<double>> &parameter_matrix,
                   const std::vector<bool> &is_positive_drift)
    {

        if (m_n_accumulator != parameter_matrix[0].size())
        {
            m_n_accumulator = parameter_matrix[0].size();

            m_denom.resize(m_n_accumulator);
            m_is_positive_drift.resize(m_n_accumulator);
            m_t0_actual.resize(m_n_accumulator);
        }

        m_A = parameter_matrix[0];
        m_b = parameter_matrix[1];
        m_mean_v = parameter_matrix[2];
        m_sd_v = parameter_matrix[3];
        m_st0 = parameter_matrix[4];
        m_t0 = parameter_matrix[5];

        m_is_positive_drift = is_positive_drift;

        for (size_t i = 0; i < m_n_accumulator; ++i)
        {
            m_denom[i] =
                !m_is_positive_drift[i]
                    ? 1.0
                    : std::fmax(Rf_pnorm5(m_mean_v[i] / m_sd_v[i], 0, 1, 1, 0),
                                1e-10);
            m_t0_actual[i] = m_t0[i] + m_st0[i] * Rf_runif(0.0, 1.0);
        }
    }

    bool validate_parameters(bool is_debug = false)
    {
        for (size_t i = 0; i < m_n_accumulator; ++i)
        {
            if (!is_negative(m_A[i], "A", is_debug) || m_A[i] < 0)
                return false;
            if (!is_negative(m_b[i], "b", is_debug))
                return false;
            if (m_b[i] < m_A[i])
            {
                if (is_debug)
                    Rcpp::Rcout << "b must be greater than A. b = " << m_b[i]
                                << ", A = " << m_A[i] << std::endl;
                return false;
            }
            if (!is_negative(m_sd_v[i], "sd_v", is_debug))
                return false;
            if (!is_negative(m_st0[i], "st0", is_debug))
                return false;
            if (!is_negative(m_t0[i], "t0", is_debug))
                return false;
        }

        return true;
    }

    void set_times(const std::vector<double> &time_parameter)
    {
        m_min_dt = time_parameter[0];
        m_max_dt = time_parameter[1];
        m_time_step = time_parameter[2];
    }
    void print_parameters(std::string str = "")
    {
        std::ios oldState(nullptr);
        oldState.copyfmt(Rcpp::Rcout);

        unsigned int sep_size = 8;

        Rcpp::Rcout << str << std::endl;
        // Set fixed-point notation with 2 decimal places
        Rcpp::Rcout << std::fixed << std::setprecision(2);

        Rcpp::Rcout << std::setw(6) << "Accu" << std::setw(sep_size) << "A"
                    << std::setw(sep_size) << "b" << std::setw(sep_size)
                    << "mean_v" << std::setw(sep_size) << "sd_v"
                    << std::setw(sep_size) << "st0" << std::setw(sep_size)
                    << "t0" << std::setw(sep_size) << "t0+st0"
                    << std::setw(sep_size) << "+drift?" << std::setw(sep_size)
                    << "denom" << std::endl;
        for (size_t i = 0; i < m_n_accumulator; ++i)
        {
            // Print values
            Rcpp::Rcout << std::setw(6) << i << std::setw(sep_size) << m_A[i]
                        << std::setw(sep_size) << m_b[i] << std::setw(sep_size)
                        << m_mean_v[i] << std::setw(sep_size) << m_sd_v[i]
                        << std::setw(sep_size) << m_st0[i]
                        << std::setw(sep_size) << m_t0[i] << std::setw(sep_size)
                        << m_t0_actual[i] << std::setw(sep_size)
                        << m_is_positive_drift[i] << std::setw(sep_size)
                        << m_denom[i] << std::endl;
        }
        Rcpp::Rcout << "[min, max, time_step]: [" << m_min_dt << ", "
                    << m_max_dt << ", " << m_time_step << "]" << std::endl;

        // Restore original formatting
        Rcpp::Rcout.copyfmt(oldState);
        Rcpp::Rcout << std::endl;
    }
    lba_class(const std::vector<std::vector<double>> &parameter_matrix =
                  {{1.2, 1.2},
                   {1.5, 1.5},
                   {2.4, 2.0},
                   {1.0, 1.0},
                   {0.0, 0.0},
                   {0.05, 0.05}},
              const std::vector<bool> &is_positive_drift = {true, true},
              const std::vector<double> &time_parameter = {0, 10, 0.01})
    {
        m_n_accumulator = parameter_matrix[0].size();
        m_denom.resize(m_n_accumulator);
        m_t0_actual.resize(m_n_accumulator);

        m_pdf = NA_REAL;
        m_cdf = NA_REAL;

        set_parameters(parameter_matrix, is_positive_drift);
        set_times(time_parameter);
    }

    ~lba_class(){};

    void d(double x)
    {
        m_dt = x - m_t0_actual[0];

        if (m_dt < 0)
        {
            m_pdf = 1e-10;
        }
        else if (m_A[0] < 1e-10)
        {
            // x, 0, 1, false
            double term2 = Rf_dnorm4(m_b[0] / m_dt, m_mean_v[0], m_sd_v[0], 0) /
                           m_denom[0];
            double term3 = (m_b[0] / (m_dt * m_dt)) * term2;
            m_pdf = std::fmax(1e-10, term3);
        }
        else
        {
            double ts = m_dt * m_sd_v[0];   // zs
            double tv = m_dt * m_mean_v[0]; // zu

            //  x, mean = 0.0, sd = 1.0, lower_tail = true, log_p = false;
            double term4 = m_mean_v[0] *
                           (Rf_pnorm5((m_b[0] - tv) / ts, 0, 1, 1, 0) -
                            Rf_pnorm5((m_b[0] - m_A[0] - tv) / ts, 0, 1, 1, 0));

            // x, mean = 0, sd = 1, log_p = FALSE
            double term5 =
                m_sd_v[0] * (Rf_dnorm4((m_b[0] - m_A[0] - tv) / ts, 0, 1, 0) -
                             Rf_dnorm4((m_b[0] - tv) / ts, 0, 1, 0));

            m_pdf = std::fmax(1e-10, (term4 + term5) / (m_A[0] * m_denom[0]));
        }

        m_pdf = std::isnan(m_pdf) ? 1e-10 : m_pdf;
    }
    void d(double x, size_t acc_idx)
    {
        m_dt = x - m_t0_actual[acc_idx];

        if (m_dt < 0)
        {
            m_pdf = 1e-10;
        }
        else if (m_A[acc_idx] < 1e-10)
        {
            double term0 = Rf_dnorm4(m_b[acc_idx] / m_dt, m_mean_v[acc_idx],
                                     m_sd_v[acc_idx], 0) /
                           m_denom[acc_idx];
            m_pdf = std::fmax(1e-10, (m_b[acc_idx] / (m_dt * m_dt)) * term0 /
                                         m_denom[acc_idx]);
        }
        else
        {
            double ts = m_dt * m_sd_v[acc_idx];
            double tv = m_dt * m_mean_v[acc_idx];
            double b_prime = m_b[acc_idx] - tv;
            double bA_prime = b_prime - m_A[acc_idx];

            double term1 =
                m_mean_v[acc_idx] * (Rf_pnorm5(b_prime / ts, 0, 1, 1, 0) -
                                     Rf_pnorm5(bA_prime / ts, 0, 1, 1, 0));
            double term2 =
                m_sd_v[acc_idx] * (Rf_dnorm4(bA_prime / ts, 0, 1, 0) -
                                   Rf_dnorm4(b_prime / ts, 0, 1, 0));

            m_pdf = std::fmax(1e-10, (term1 + term2) /
                                         (m_A[acc_idx] * m_denom[acc_idx]));
        }

        m_pdf = std::isnan(m_pdf) ? 1e-10 : m_pdf;
    }

    void p(double x)
    {
        // Starting from the 2nd accumulator.
        static thread_local double last_x =
            -std::numeric_limits<double>::infinity();
        static thread_local double last_t0 =
            -std::numeric_limits<double>::infinity();
        static thread_local double cached_dt = 0.0;

        for (size_t i = 1; i < m_n_accumulator; ++i)
        {
            // Check if we can reuse previous calculation
            if (x != last_x || m_t0_actual[i] != last_t0)
            {
                m_dt = x - m_t0_actual[i];
                last_x = x;
                last_t0 = m_t0_actual[i];
                cached_dt = m_dt;
            }
            else
            {
                m_dt = cached_dt; // Reuse cached value
            }

            if (m_dt < 0.0)
            {
                m_cdf = 1e-10; // 0.0;
                               // m_cdf = 0;
            }
            else if (m_A[i] < 1e-10)
            {
                m_cdf = Rf_pnorm5(m_b[i] / m_dt, m_mean_v[i], m_sd_v[i], false,
                                  false) /
                        m_denom[i];
                m_cdf = std::clamp(m_cdf, 1e-10, 1.0);
            }
            else
            {
                double ts = m_dt * m_sd_v[i];
                double tv = m_dt * m_mean_v[i];

                double b_prime = m_b[i] - tv;
                double bA_prime = b_prime - m_A[i];

                double term1 =
                    bA_prime * Rf_pnorm5(bA_prime / ts, 0.0, 1.0, 1, 0);
                double term2 =
                    b_prime * Rf_pnorm5(b_prime / ts, 0.0, 1.0, 1, 0);
                double term3 = ts * (Rf_dnorm4(bA_prime / ts, 0.0, 1.0, 0) -
                                     Rf_dnorm4(b_prime / ts, 0.0, 1.0, 0));

                m_cdf = (1.0 + (term1 - term2 + term3) / m_A[i]) / m_denom[i];
                m_cdf = std::clamp(m_cdf, 1e-10, 1.0);
            }

            const double temp = m_pdf * (1.0 - m_cdf);
            m_pdf = std::isnan(temp) ? 1e-10 : temp;

        } // end of for loop over accumulators
    }
    void p(double x, size_t acc_idx)
    {
        // Static variables to cache previous values
        static thread_local double last_x =
            -std::numeric_limits<double>::infinity();
        static thread_local double last_t0 =
            -std::numeric_limits<double>::infinity();
        static thread_local double cached_dt = 0.0;

        // Check if we can reuse previous calculation
        if (x != last_x || m_t0_actual[acc_idx] != last_t0)
        {
            m_dt = x - m_t0_actual[acc_idx];
            last_x = x;
            last_t0 = m_t0_actual[acc_idx];
            cached_dt = m_dt;
        }
        else
        {
            m_dt = cached_dt; // Reuse cached value
        }

        if (m_dt < 0.0)
        {
            m_cdf = 1e-10;
        }
        else if (m_A[acc_idx] < 1e-10)
        {
            m_cdf = Rf_pnorm5(m_b[acc_idx] / m_dt, m_mean_v[acc_idx],
                              m_sd_v[acc_idx], false, false) /
                    m_denom[acc_idx];
            m_cdf = std::clamp(m_cdf, 0.0, 1.0);
        }
        else
        {
            double ts = m_dt * m_sd_v[acc_idx];
            double tv = m_dt * m_mean_v[acc_idx];
            double b_prime = m_b[acc_idx] - tv;
            double bA_prime = b_prime - m_A[acc_idx];

            double term1 = bA_prime * Rf_pnorm5(bA_prime / ts, 0.0, 1.0, 1, 0);
            double term2 = b_prime * Rf_pnorm5(b_prime / ts, 0.0, 1.0, 1, 0);
            double term3 = ts * (Rf_dnorm4(bA_prime / ts, 0.0, 1.0, 0) -
                                 Rf_dnorm4(b_prime / ts, 0.0, 1.0, 0));

            m_cdf = (1.0 + (term1 - term2 + term3) / m_A[acc_idx]) /
                    m_denom[acc_idx];
            m_cdf = std::clamp(m_cdf, 0.0, 1.0);
        }
    }

    void fptcdf(double x)
    {
        m_dt = x - m_t0_actual[0];

        if (m_dt < 0.0)
        {
            m_cdf = 0.0; // 0.0;
        }
        else if (m_A[-0] < 1e-10)
        {
            m_cdf =
                Rf_pnorm5(m_b[0] / m_dt, m_mean_v[0], m_sd_v[0], false, false) /
                m_denom[0];
            m_cdf = std::clamp(m_cdf, 0.0, 1.0);
        }
        else
        {
            double ts = m_dt * m_sd_v[0];
            double tv = m_dt * m_mean_v[0];

            double b_prime = m_b[0] - tv;
            double bA_prime = b_prime - m_A[0];

            double term1 = bA_prime * Rf_pnorm5(bA_prime / ts, 0.0, 1.0, 1, 0);
            double term2 = b_prime * Rf_pnorm5(b_prime / ts, 0.0, 1.0, 1, 0);
            double term3 = ts * (Rf_dnorm4(bA_prime / ts, 0.0, 1.0, 0) -
                                 Rf_dnorm4(b_prime / ts, 0.0, 1.0, 0));

            m_cdf = (1.0 + (term1 - term2 + term3) / m_A[0]) / m_denom[0];
            m_cdf = std::clamp(m_cdf, 0.0, 1.0);
        }
    }

    void r(std::vector<std::pair<unsigned int, double>> &out)
    {
        unsigned int n = out.size();
        std::vector<double> lowers(m_n_accumulator);
        for (size_t i = 0; i < m_n_accumulator; ++i)
        {
            lowers[i] = m_is_positive_drift[i] ? 0.0 : R_NegInf;
        }

        // Pre-generate all random numbers
        std::vector<double> runifs(n * m_n_accumulator);
        for (auto &x : runifs)
            x = Rf_runif(0.0, 1.0);

        for (size_t i = 0; i < n; ++i)
        {
            double min_rt = R_PosInf;
            unsigned int min_index = 0;

            for (size_t accu_idx = 0; accu_idx < m_n_accumulator; ++accu_idx)
            {
                tnorm::tnorm_class tnorm_obj(m_mean_v[accu_idx],
                                             m_sd_v[accu_idx], lowers[accu_idx],
                                             R_PosInf);

                double tmp_rt =
                    m_t0_actual[accu_idx] +
                    (m_b[accu_idx] -
                     m_A[accu_idx] * runifs[i * m_n_accumulator + accu_idx]) /
                        (tnorm_obj.r());

                if (!std::isfinite(tmp_rt))
                {
                    throw std::runtime_error(
                        "LBA's r function generates non-finite response time.");
                }

                if (tmp_rt < min_rt)
                {
                    min_rt = tmp_rt;
                    min_index = accu_idx;
                }
            }
            out[i] = {min_index, min_rt};
        }
    }
    void r_inverse(std::vector<std::pair<unsigned int, double>> &out)
    {
        unsigned int n = out.size();

        // 1. Get CDFs on a fine grid
        auto cdfs = theoretical_plba();

        // Precompute max CDF values and total probability
        std::vector<double> acc_max_cdf(m_n_accumulator);
        double total_prob = 0.0;
        for (size_t acc = 0; acc < m_n_accumulator; ++acc)
        {
            acc_max_cdf[acc] = cdfs[acc].back();
            total_prob += acc_max_cdf[acc];
        }

        // Normalization factor (handle floating-point imprecision)
        const double norm_factor = 1.0 / total_prob;

        for (size_t i = 0; i < n; ++i)
        {
            double u = Rf_runif(0.0, 1.0);

            // Clamp u to [0, 1] and account for normalization
            u = std::min(std::max(u, 0.0), 1.0) * total_prob * norm_factor;

            // 2. Find which accumulator's CDF contains the random number
            size_t chosen_acc = 0;
            double running_sum = 0.0;
            for (; chosen_acc < m_n_accumulator; ++chosen_acc)
            {
                running_sum += acc_max_cdf[chosen_acc];
                if (u <= running_sum)
                    break;
            }

            // allback to last accumulator if numerical issues persist
            if (chosen_acc >= m_n_accumulator)
            {
                chosen_acc = m_n_accumulator - 1;
                u = running_sum; // Use the maximum possible value
            }

            // 3. Find the corresponding RT via inverse transform
            double rt = m_max_dt + m_t0_actual[0]; // Default fallback
            if (!cdfs[chosen_acc].empty())
            {
                auto &cdf = cdfs[chosen_acc];
                auto it = std::lower_bound(
                    cdf.begin(), cdf.end(),
                    u - (running_sum - acc_max_cdf[chosen_acc]));
                size_t idx = std::distance(cdf.begin(), it);

                if (idx == 0)
                {
                    rt = m_rt_grid[chosen_acc][0];
                }
                else if (idx < cdf.size())
                {
                    double p0 = cdf[idx - 1];
                    double p1 = cdf[idx];
                    double t0 = m_rt_grid[chosen_acc][idx - 1];
                    double t1 = m_rt_grid[chosen_acc][idx];
                    rt =
                        t0 + (t1 - t0) *
                                 ((u - (running_sum - acc_max_cdf[chosen_acc]) -
                                   p0) /
                                  (p1 - p0));
                }
            }
            else
            {
                print_parameters("Unknown numerical issues occur");
                Rcpp::Rcout << "chosen_acc = " << chosen_acc
                            << " total_prob = " << total_prob
                            << " runif() = " << u << std::endl;
                rt = m_max_dt + m_t0_actual[0];
                chosen_acc = 0;
            }

            out[i] = {chosen_acc, rt};
        }
    }

    std::vector<double> dlba(const std::vector<double> &rt)
    {
        std::vector<double> out(rt.size());

        for (size_t i = 0; i < rt.size(); ++i)
        {
            d(rt[i]);
            p(rt[i]);
            out[i] = m_pdf;
        }

        return out;
    }

    std::vector<std::vector<double>> dlba_all(const std::vector<double> &rt)
    {
        std::vector<std::vector<double>> out(m_n_accumulator,
                                             std::vector<double>(rt.size()));
        std::vector<double> survivor_products(m_n_accumulator);

        for (size_t i = 0; i < rt.size(); ++i)
        {
            // Calculate PDFs for all accumulators
            for (size_t acc_idx = 0; acc_idx < m_n_accumulator; ++acc_idx)
            {
                d(rt[i], acc_idx);
                out[acc_idx][i] = m_pdf;
            }

            // For non-defective densities, multiply each PDF by the probability
            // accumulators of other being slower
            for (size_t acc_idx = 0; acc_idx < m_n_accumulator; ++acc_idx)
            {
                double survivor_product = 1.0;
                for (size_t other_idx = 0; other_idx < m_n_accumulator;
                     ++other_idx)
                {
                    if (other_idx != acc_idx)
                    {
                        p(rt[i], other_idx);
                        survivor_product *= (1.0 - m_cdf);
                    }
                }
                out[acc_idx][i] *= survivor_product;
            }
        }
        return out;
    }
    std::vector<std::vector<double>> plba_all(const std::vector<double> &rt)
    {
        auto cdfs = theoretical_plba();

        // Check if we need normalization
        double total_max = 0.0;
        for (size_t acc_idx = 0; acc_idx < m_n_accumulator; ++acc_idx)
        {
            // add up the last element of each accumulator.
            total_max += cdfs[acc_idx].back();
        }
        Rcpp::Rcout << "total_max = " << total_max << std::endl;

        // Normalise if the sum exceeds 1 (with small tolerance)
        if (total_max > 1.0001)
        { // Using small tolerance for floating point comparison
            double scale_factor = 1.0 / total_max;
            Rcpp::Rcout << "Normalising using the scale factor = "
                        << scale_factor << std::endl;

            for (size_t acc_idx = 0; acc_idx < m_n_accumulator; ++acc_idx)
            {
                for (auto &val : cdfs[acc_idx])
                {
                    val *= scale_factor;
                }
            }
        }

        // Interpolate to get CDF at requested RTs
        std::vector<std::vector<double>> out(m_n_accumulator,
                                             std::vector<double>(rt.size()));

        for (size_t acc_idx = 0; acc_idx < m_n_accumulator; ++acc_idx)
        {
            for (size_t i = 0; i < rt.size(); ++i)
            {
                double t = rt[i];
                if (t <= m_rt_grid[acc_idx][0])
                {
                    out[acc_idx][i] = 0.0;
                }
                else if (t >= m_rt_grid[acc_idx].back())
                {
                    out[acc_idx][i] = 1.0;
                }
                else
                {
                    // Find the interval where rt[i] lies
                    auto it = std::lower_bound(m_rt_grid[acc_idx].begin(),
                                               m_rt_grid[acc_idx].end(), t);
                    size_t idx = std::distance(m_rt_grid[acc_idx].begin(), it);
                    double t0 = m_rt_grid[acc_idx][idx - 1],
                           t1 = m_rt_grid[acc_idx][idx];
                    double cdf0 = cdfs[acc_idx][idx - 1],
                           cdf1 = cdfs[acc_idx][idx];
                    // Linear interpolation
                    out[acc_idx][i] =
                        cdf0 + (cdf1 - cdf0) * (t - t0) / (t1 - t0);
                }
            }
        }

        return out;
    }

    std::vector<std::vector<double>> theoretical_dlba()
    {
        // Calculate number of points based on time_step
        create_time_grid();

        std::vector<std::vector<double>> out(m_n_accumulator,
                                             std::vector<double>(m_n_point));

        for (size_t i = 0; i < m_n_point; ++i)
        {
            // Calculate PDFs for all accumulators
            for (size_t acc_idx = 0; acc_idx < m_n_accumulator; ++acc_idx)
            {
                d(m_rt_grid[acc_idx][i], acc_idx);
                out[acc_idx][i] = m_pdf;
            }

            // For non-defective densities, multiply each PDF by the probability
            // accumulators of other being slower
            for (size_t acc_idx = 0; acc_idx < m_n_accumulator; ++acc_idx)
            {
                double survivor_product = 1.0;
                for (size_t other_idx = 0; other_idx < m_n_accumulator;
                     ++other_idx)
                {
                    if (other_idx != acc_idx)
                    {
                        p(m_rt_grid[acc_idx][i], other_idx);
                        survivor_product *= (1.0 - m_cdf);
                    }
                }
                out[acc_idx][i] *= survivor_product;
            }
        }

        return out;
    }
    std::vector<std::vector<double>> theoretical_plba()
    {
        auto pdfs = theoretical_dlba();

        std::vector<std::vector<double>> cdfs(
            m_n_accumulator, std::vector<double>(m_n_point, 0.0));

        for (size_t acc_idx = 0; acc_idx < m_n_accumulator; ++acc_idx)
        {
            // Numerically integrate PDFs to get CDFs (trapezoidal rule)
            for (size_t i = 0; i < m_n_point; ++i)
            {
                double dt = m_rt_grid[acc_idx][i] - m_rt_grid[acc_idx][i - 1];
                cdfs[acc_idx][i] =
                    cdfs[acc_idx][i - 1] +
                    0.5 * (pdfs[acc_idx][i - 1] + pdfs[acc_idx][i]) * dt;
            }
        }

        return cdfs;
    }

    // Get density for a specific RT and accumulator
    std::vector<double> dlba_inverse(const std::vector<double> &rts,
                                     const std::vector<unsigned int> &responses)
    {
        std::vector<double> densities(rts.size());

        for (size_t i = 0; i < rts.size(); ++i)
        {
            if (responses[i] < 0 || responses[i] >= m_n_accumulator)
            {
                densities[i] = 0.0;
                continue;
            }

            densities[i] = get_density(rts[i], responses[i]);
        }

        return densities;
    }
};

} // namespace lba
