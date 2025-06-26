#pragma once
#include <RcppArmadillo.h>

namespace tnorm
{
static constexpr double SQRT_2PI = 2.506628274631000502415765284811;

class tnorm_class
{
  public:
    double m_mean, m_sd, m_lower, m_upper; // mean, sd, lower, upper,
    bool m_lower_tail, m_log_p;            // lower tail, log probability
    double m_denom, m_log_denom;

    void validate_parameters(double lower, double upper, double mean, double sd)
    {
        if (upper < lower)
        {
            Rcpp::Rcout << "upper = " << upper << ", lower = " << lower
                        << std::endl;
            throw std::runtime_error("'upper' must be greater than 'lower'.");
        }
        if (sd <= 0 || sd == R_PosInf)
        {
            throw std::runtime_error("The standard deviation must be a finite "
                                     "value, gretaer than 0");
        }
        if (mean == R_NegInf || mean == R_PosInf)
        {
            throw std::runtime_error("The mean must have a finite value.");
        }
    }
    bool validate_parameters(double x, double mean, double sd)
    {
        bool is_valid = true;
        if (sd == R_PosInf || sd < 0 || mean == R_NegInf || mean == R_PosInf)
        {
            is_valid = false;
        }

        if (x < m_lower || x > m_upper)
        {
            is_valid = false;
            // return m_log_p ? R_NegInf : 1e-10;
        }

        return is_valid;
    }

    void set_parameters(double mean, double sd, double lower, double upper,
                        bool log_p)
    {
        this->m_mean = mean;
        this->m_sd = sd;
        this->m_lower = lower;
        this->m_upper = upper;
        this->m_log_p = log_p;
    }
    void set_parameters(double mean, double sd)
    {
        this->m_mean = mean;
        this->m_sd = sd;
        //  x, mean = 0.0, sd = 1.0, lower_tail = true, log_p = false;
        this->m_denom = Rf_pnorm5(m_upper, m_mean, m_sd, 1, 0) -
                        Rf_pnorm5(m_lower, m_mean, m_sd, 1, 0);
        this->m_log_denom = std::log(m_denom);
    }

    void printParameters()
    {
        Rcpp::Rcout << "[mean, sd, lower, upper, log_p]" << this->m_mean << ", "
                    << this->m_sd << ", " << this->m_lower << ", "
                    << this->m_upper << ", " << this->m_log_p << ", "
                    << std::endl;
    }

    /* --- A no-argument, default constructor for the std::vector to work --- */
    tnorm_class() : m_mean(0), m_sd(0), m_lower(0), m_upper(0), m_log_p(false)
    {
        m_lower_tail = false;
        m_denom = 0;
        m_log_denom = -std::numeric_limits<double>::infinity();
    }

    /* ------------------- ptnorm ------------------------------------- */
    tnorm_class(double mean, double sd, double lower, double upper,
                bool lower_tail, bool log_p)
        : m_mean(mean), m_sd(sd), m_lower(lower), m_upper(upper),
          m_lower_tail(lower_tail), m_log_p(log_p)
    {
        validate_parameters(m_lower, m_upper, m_mean, m_sd);
        m_denom = 0;
        m_log_denom = -std::numeric_limits<double>::infinity();
    }

    /* ------------------- dtnorm (no lower_tail-------------------------- */
    tnorm_class(double mean, double sd, double lower, double upper, bool log_p)
        : m_mean(mean), m_sd(sd), m_lower(lower), m_upper(upper), m_log_p(log_p)
    {
        validate_parameters(m_lower, m_upper, m_mean, m_sd);
        // Precompute denominator
        m_denom = Rf_pnorm5(m_upper, m_mean, m_sd, 1, 0) -
                  Rf_pnorm5(m_lower, m_mean, m_sd, 1, 0);
        m_log_denom = std::log(m_denom);
    }

    /* ------------------- rtnorm ------------------------------------- */
    tnorm_class(double mean, double sd, double lower, double upper)
        : m_mean(mean), m_sd(sd), m_lower(lower), m_upper(upper)
    {
        validate_parameters(m_lower, m_upper, m_mean, m_sd);

        m_lower_tail = false;
        m_denom = 0;
        m_log_denom = -std::numeric_limits<double>::infinity();
    }

    double d(double x) const
    {
        if (x < m_lower || x > m_upper)
        {
            return m_log_p ? R_NegInf : 1e-10;
        }

        double numer = Rf_dnorm4(x, m_mean, m_sd, m_log_p);
        return m_log_p ? numer - m_log_denom : numer / m_denom;
    }
    double p(double q)
    {
        double out, qtmp;

        if (m_lower_tail)
        {
            out = (q < m_lower) ? 0 : 1;
        }
        else
        {
            out = (q < m_lower) ? 1 : 0;
        }

        if ((q >= m_lower) && (q <= m_upper))
        {
            double term1 = Rf_pnorm5(m_upper, m_mean, m_sd, 1, 0);
            double term2 = Rf_pnorm5(m_lower, m_mean, m_sd, 1, 0);
            double term3 = Rf_pnorm5(q, m_mean, m_sd, 1, 0);
            qtmp = m_lower_tail ? (term3 - term2) : (term1 - term3);
            out = m_log_p ? (std::log(qtmp) - m_log_denom) : (qtmp / m_denom);
        }

        return out;
    }

    // Vectorised functions
    void p(std::vector<double> &q, std::vector<double> &out)
    {
        out.resize(q.size()); // Ensure out has the same size as x
        std::transform(q.begin(), q.end(), out.begin(),
                       [this](double val) { return p(val); });
    }
    void d(const std::vector<double> &x, std::vector<double> &out)
    {
        out.resize(x.size()); // Ensure out has the same size as x
        std::transform(x.begin(), x.end(), out.begin(),
                       [this](double val) { return d(val); });
    }

    // rtnorm internal
    // Accept-Reject Algorithm 0; Naive method A-R method
    double rtnorm0(const double &l, const double &u)
    {
        bool invalid = true;
        double z;
        while (invalid)
        {
            z = Rf_rnorm(0.0, 1.0);
            if (z <= u && z >= l)
                break;
        }
        return z;
    }

    // Algorithm 1; 'expl'; use when lower > mean; upper = INFINITY; p 122,
    // right
    double rtnorm1(const double &l, const double &u)
    {
        bool invalid = true;
        double z, r, num; // a stands for alphaStar in Robert (1995)
        double a = 0.5 * (std::sqrt(l * l + 4.0) + l);

        while (invalid)
        {
            z = (-1.0 / a) * std::log(Rf_runif(0.0, 1.0)) +
                l; // control lower boundary
            num = Rf_runif(0.0, 1.0);
            r = std::exp(-0.5 * (z - a) * (z - a));
            if (num <= r && z <= u)
                break;
        }

        return z;
    }

    // Algorithm 2; 'expu'; use when upper < mean; lower = -INFINITY.
    double rtnorm2(const double &l, const double &u)
    {
        bool invalid = true;
        double z, r, num;
        double a = 0.5 * (std::sqrt(u * u + 4.0) - u);

        while (invalid)
        {
            z = (-1.0 / a) * std::log(Rf_runif(0.0, 1.0)) -
                u; // control lower boundary
            num = Rf_runif(0.0, 1.0);
            r = std::exp(-0.5 * (z - a) * (z - a));
            if (num <= r && z <= -l)
                break;
        }

        return -z; // note the negative
    }

    // Algorithm 3; page 123. 2.2. Two-sided truncated normal dist.
    double rtnorm3(const double &l, const double &u)
    {
        bool invalid = true;
        double z, r, num; // a stands for alphaStar in Robert (1995)
        double l2 = l * l;
        double u2 = u * u;

        while (invalid)
        {
            z = Rf_runif(l, u);
            if (l > 0)
            {
                r = std::exp(0.5 * (l2 - z * z));
            }
            else if (u < 0)
            {
                r = std::exp(0.5 * (u2 - z * z));
            }
            else
            {
                r = std::exp(-0.5 * z * z);
            }
            num = Rf_runif(0.0, 1.0);
            if (num <= r)
                break;
        }

        return z;
    }

    double r()
    {
        // auto &rng = random_class::get_instance().get_rng();

        // Standardised lower and upper
        double z, stdl, stdl2, stdu, stdu2, eq_a1, eq_a2;
        bool a0, a1, a2;
        stdl = (m_lower - m_mean) /
               m_sd; // stdl = standard lower, stdu = standard upper
        stdu = (m_upper - m_mean) / m_sd;
        stdl2 = stdl * stdl;
        stdu2 = stdu * stdu;

        // Accept-Reject Algorithm 0;
        // Algorithm (1): Use Proposition 2.3 with only lower truncation.
        // upper==INFINITY
        // rejection sampling with exponential proposal. Use if lower > mean
        //
        // Algorithm (2): Use -x ~ N_+ (-mu, -mu^+, sigma^2) on page 123.
        // lower==-INFINITY
        // rejection sampling with exponential proposal. Use if upper < mean.
        //
        // Algorithm (3, else): rejection sampling with uniform proposal.
        // Use if bounds are narrow and central.
        a0 = (stdl < 0 && m_upper == R_PosInf) ||
             (stdl == R_NegInf && stdu > 0) ||
             (std::isfinite(stdl) && std::isfinite(m_upper) && stdl < 0 &&
              stdu > 0 && (stdu - stdl) > SQRT_2PI);

        eq_a1 =
            stdl +
            (2.0 * std::sqrt(M_E) / (stdl + std::sqrt(stdl2 + 4.0))) *
                (std::exp(0.25 * (2.0 * stdl - stdl * std::sqrt(stdl2 + 4.0))));

        a1 = (stdl >= 0) && (stdu > eq_a1);

        eq_a2 =
            -stdu +
            (2.0 * std::sqrt(M_E) / (-stdu + std::sqrt(stdu2 + 4.0))) *
                (std::exp(0.25 * (2.0 * stdu + stdu * std::sqrt(stdu2 + 4.0))));

        a2 = (stdu <= 0) && (-stdl > eq_a2);

        if (a0)
        {
            z = rtnorm0(stdl, stdu);
        }
        else if (a1)
        {
            z = rtnorm1(stdl, stdu);
        }
        else if (a2)
        {
            z = rtnorm2(stdl, stdu);
        }
        else
        {
            z = rtnorm3(stdl, stdu);
        }

        return z * m_sd + m_mean;
    }
};

} // namespace tnorm