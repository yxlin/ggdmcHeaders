#pragma once
#include <RcppArmadillo.h>

namespace ddm
{
namespace
{
enum ParamIndex
{
    PARAM_a = 0,
    PARAM_d,
    PARAM_precision,
    PARAM_s,
    PARAM_st0,
    PARAM_sv,
    PARAM_sz,
    PARAM_t0,
    PARAM_v,
    PARAM_z,
    NUM_PARAMS // Optional: total number of parameters
};
// Core parameters: a      d       precision       s       st0     sv      sz t0
// v      z
const double EPSILON = 1e-6;

class ddm_class
{
  public:
    double m_a;  // Boundary separation
    double m_v;  // Mean of the drift
    double m_t0; // Non-decision time
    double m_d;  // Difference between boundaries of non-decision time

    // the user enters the reltive szr
    double m_szr; // width of zr distribution
    double m_sv;  // standard deviation of v distribution
    double m_st0; // width of t0 distribution
    double m_zr;  // Mean of diffusion starting point relative to boundaries
    double m_s;   // standard deviation; sqrt(diffusion constant)
    double m_precision; // precision
    // m_plus = true uses g_plus (upper); other use g_minus (lower)
    double m_plus;

    double m_t_offset; // actual t0
    double m_z, m_sz;

    double TUNE_DZ;
    double TUNE_DV;
    double TUNE_DT0;

    // If std=c++11 we can use C++ defaults to set as = 1e-6;
    double TUNE_PDE_DT_MIN;
    double TUNE_PDE_DT_MAX;   // ... we can default to = 1e-6;
    double TUNE_PDE_DT_SCALE; // ... we can default to = 0.0;
    double TUNE_INT_T0;
    double TUNE_INT_Z;
    double TUNE_SV_EPSILON;
    double TUNE_SZ_EPSILON;
    double TUNE_ST0_EPSILON;

    double m_t_min, m_t_max, m_time_step;

    ddm_class(const std::vector<double> &parameters = {1.0, 1.5, 0.5, 0.0, 0.0,
                                                       0.0, 0.1, 0.0, 1.0, 2.5},
              bool is_minus = true,
              const std::vector<double> &time_parameter = {-0.6, 0.6, 0.001});

    void set_parameters(const std::vector<double> &parameters,
                        bool is_minus = true);

    void set_parameters(const std::vector<std::vector<double>> &parameters,
                        bool is_minus = true);

    // void set_parameters2(const std::vector<std::vector<double>> &parameters,
    //                      bool is_minus = true);
    bool validate_parameters(bool is_debug = false);
    void print(std::string header = "") const;

    // Density function
    double compute_g_series(double ta2, double zr, bool small_time, int N);
    // double compute_g_factor(double dt, double a, double zr, double v,
    // bool no_var_case);

    double compute_g_factor(double dt, double zr, bool no_var_case);

    std::pair<int, int> get_N(double dt, double ta2, double eps);

    double g_no_var(double dt, double zr);
    // double g_no_var(double dt, double a, double zr, double v);
    // double integral_v(double dt, double a, double zr, double v);
    double integral_v(double dt, double zr);

    double integrate_v_over_zr(double dt);
    // double integrate_v_over_zr(double dt, double zr);
    double integral_z(double dt);

    double integrate_z_over_t(double dt);
    double integral_t0(double dt);
    double g(double rt);

    std::vector<double> dddm(const std::vector<double> &rt);
    void set_times(const std::vector<double> &time_parameter);

  private:
    void set_precision();

    // Precomputed constants
    const double M_PI_SQUARED = M_PI * M_PI;
    double m_a2, m_v2, m_sv2;
};

// Implementation --------------------------
void ddm_class::set_times(const std::vector<double> &time_parameter)
{
    m_t_min = time_parameter[0];
    m_t_max = time_parameter[1];
    m_time_step = time_parameter[2];
}
void ddm_class::set_precision()
{
    // Try to achieve an accuracy of approximately 10^{-p} for the CDF.
    TUNE_PDE_DT_MIN = std::pow(10, -0.400825 * m_precision - 1.422813);
    TUNE_PDE_DT_MAX = std::pow(10, -0.627224 * m_precision + 0.492689);
    TUNE_PDE_DT_SCALE = std::pow(10, -1.012677 * m_precision + 2.261668);

    TUNE_DZ = std::pow(10, -0.5 * m_precision - 0.033403); // CDF
    TUNE_DV = std::pow(10, -1.0 * m_precision + 1.4);
    TUNE_DT0 = std::pow(10, -0.5 * m_precision - 0.323859);

    TUNE_INT_T0 = 0.089045 * std::exp(-1.037580 * m_precision); // PDF
    TUNE_INT_Z = 0.508061 * std::exp(-1.022373 * m_precision);

    TUNE_SV_EPSILON = std::pow(10, -(m_precision + 2.0)); // Used by pdiffusion

    // Used by ddiffusion and pdiffusion
    TUNE_SZ_EPSILON = std::pow(10, -(m_precision + 2.0));
    TUNE_ST0_EPSILON = std::pow(10, -(m_precision + 2.0));
}

ddm_class::ddm_class(const std::vector<double> &parameters, bool is_minus,
                     const std::vector<double> &time_parameter)
{
    set_parameters(parameters, is_minus);
    set_times(time_parameter);
}

void ddm_class::set_parameters(const std::vector<double> &parameters,
                               bool is_minus)
{
    // Use by dfastdm, pfastdm and rfastdm.
    // Store common expressions
    double s = parameters[PARAM_s];
    double scaling = (s == 1) ? 1.0 : (1.0 / s);

    // Used in F_calculator
    m_s = s;
    m_precision = parameters[PARAM_precision];
    m_d = parameters[PARAM_d];

    m_st0 = parameters[PARAM_st0];

    m_a = parameters[PARAM_a] * scaling;
    m_sv = parameters[PARAM_sv] * scaling;
    m_v = is_minus ? parameters[PARAM_v] * scaling
                   : (-parameters[PARAM_v]) * scaling;

    m_t_offset = parameters[PARAM_t0] + 0.5 * m_st0;
    m_t0 = (m_st0 == 0) ? parameters[PARAM_t0] : m_t_offset;

    m_plus = is_minus ? false : true;

    double tmp_zr = parameters[PARAM_z] / m_a;
    m_zr = is_minus ? tmp_zr : 1 - tmp_zr;

    m_sz = parameters[PARAM_sz];
    m_szr = m_sz / m_a;
    set_precision();

    // Precompute trivial constants
    m_a2 = m_a * m_a;
    m_v2 = m_v * m_v;
    m_sv2 = m_sv * m_sv;
}

void ddm_class::set_parameters(
    const std::vector<std::vector<double>> &parameters, bool is_minus)
{
    // Assuming the user enters absolute z and sz.
    // I then convert them to zr and szr, so that the user will
    // not need to know that the range of zr and szr are clamped
    // within 0 and 1.
    //
    // Store common expressions
    double s = parameters[PARAM_s][0];
    double scaling = (s == 1) ? 1.0 : (1.0 / s);

    // Used in F_calculator
    m_s = s;
    m_precision = parameters[PARAM_precision][0];
    m_d = parameters[PARAM_d][0];
    m_st0 = parameters[PARAM_st0][0];

    m_a = parameters[PARAM_a][0] * scaling;
    m_sv = parameters[PARAM_sv][0] * scaling;
    m_v = is_minus ? parameters[PARAM_v][0] * scaling
                   : (-parameters[PARAM_v][0]) * scaling;

    m_t_offset = parameters[PARAM_t0][0] + 0.5 * parameters[PARAM_st0][0];
    m_t0 = (m_st0 == 0) ? parameters[PARAM_t0][0] : m_t_offset;

    m_plus = is_minus ? false : true;

    m_z = parameters[PARAM_z][0];
    m_sz = parameters[PARAM_sz][0];

    double tmp_zr = m_z / m_a;
    m_zr = is_minus ? tmp_zr : 1 - tmp_zr;
    m_szr = m_sz / m_a;

    set_precision();

    // Precompute trivial constants
    m_a2 = m_a * m_a;
    m_v2 = m_v * m_v;
    m_sv2 = m_sv * m_sv;
}

bool ddm_class::validate_parameters(bool is_debug)
{
    bool valid = true;

    if (m_a <= 0)
    {
        valid = false;
        if (is_debug)
        {
            Rcpp::Rcout << "error: invalid parameter a = " << m_a << std::endl;
        }
    }
    if (m_szr < 0 || m_szr > 1)
    {
        valid = false;
        if (is_debug)
        {
            Rcpp::Rcout << "Warning: the variability of the relative starting "
                        << "point is not within 0 and 1, szr = " << m_szr
                        << std::endl;
        }
    }
    if (m_st0 < 0)
    {
        valid = false;
        if (is_debug)
        {
            Rcpp::Rcout << "error: invalid parameter st0 = " << m_st0
                        << std::endl;
        }
    }
    if (m_sv < 0)
    {
        valid = false;
        if (is_debug)
        {
            Rcpp::Rcout << "Warning: the parameter sv = " << m_sv
                        << " is not positive." << std::endl;
        }
    }
    if (m_t_offset < 0)
    {
        valid = false;
        if (is_debug)
        {
            // double t_offset = std::fabs(0.5 * m_d) + 0.5 * m_st0;
            Rcpp::Rcout << "Warning: t_offset is negative, t0 = " << m_t0
                        << ", t_offset = " << m_t_offset << std::endl;
        }
    }
    if (m_zr - 0.5 * m_szr <= 0)
    {
        valid = false;
        if (is_debug)
        {
            Rcpp::Rcout << "Warning: zr - szr/2 = "
                        << std::to_string(m_zr - 0.5 * m_szr)
                        << ", which must be greater than 0" << std::endl;
        }
    }
    if (m_zr + 0.5 * m_szr >= 1)
    {
        valid = false;
        if (is_debug)
        {
            Rcpp::Rcout << "Warning: zr + szr/2 = "
                        << std::to_string(m_zr + 0.5 * m_szr)
                        << ", which must be less than 1" << std::endl;
        }
    }
    if (m_s <= 0)
    {
        valid = false;
        if (is_debug)
        {
            Rcpp::Rcout << "Warning: diffusion constant must be greater "
                           "than  0, s = "
                        << m_s;
        }
    }
    return valid;
}

void ddm_class::print(std::string header) const
{
    std::string boundary = m_plus ? "upper" : "lower";

    // Header
    Rcpp::Rcout << "\n=== " << header << " ===\n";

    // Core parameters (aligned)
    Rcpp::Rcout << std::left << std::setw(6) << "a:" << std::setw(10) << m_a
                << std::setw(6) << "v:" << std::setw(10) << m_v << std::setw(6)
                << "t0:" << std::setw(10) << m_t0 << "\n"
                << std::setw(6) << "d:" << std::setw(10) << m_d << std::setw(6)
                << "zr:" << std::setw(10) << m_zr << std::setw(6)
                << "z:" << std::setw(10) << m_z << std::setw(6)
                << "sz:" << std::setw(10) << m_sz << "\n";

    // Variability parameters (aligned)
    Rcpp::Rcout << std::setw(6) << "szr:" << std::setw(10) << m_szr
                << std::setw(6) << "sv:" << std::setw(10) << m_sv
                << std::setw(6) << "st0:" << std::setw(10) << m_st0 << "\n"
                << std::setw(6) << "s:" << std::setw(10) << m_s << std::setw(6)
                << "prec:" << std::setw(10) << m_precision << std::setw(6)
                << "b:" << std::setw(10) << boundary << "\n";

    // Time parameters (aligned)
    Rcpp::Rcout << std::setw(6) << "t_min:" << std::setw(10) << m_t_min
                << std::setw(6) << "t_max:" << std::setw(10) << m_t_max
                << std::setw(6) << "dt:" << std::setw(10) << m_time_step
                << "\n\n";
}
// DDM theorteical density -----------------------------
double ddm_class::compute_g_series(double ta2, double zr, bool small_time,
                                   int N)
{
    // if (ta2 <= 0)
    // {
    //     throw std::runtime_error("dt / (a * a) be greater than 0.");
    // }

    double out = 0;
    if (small_time)
    {
        // A3 Formula on p1217; Appendix Mathematical Details V, R, & V
        // (2004) See also Feller (1971, p359 & p370 Problem 22); Note zr =
        // z/a. v = 0 & a = 1
        const double ta2_3 = ta2 * ta2 * ta2;
        const double norm_factor = 1.0 / std::sqrt(M_2PI * ta2_3);

        int start = -N / 2;
        int end = N / 2;

        for (int i = start; i <= end; i++)
        {
            const double step = 2.0 * i + zr;
            out += std::exp(-step * step / (2.0 * ta2)) * step;
        }
        return out * norm_factor;
    }
    else
    {
        // Large-time approximation (A4 formula, p1217)
        for (int i = 1; i <= N; i++)
        {
            const double step = i * M_PI;
            out += std::exp(-0.5 * step * step * ta2) * std::sin(step * zr) * i;
        }
        return out * M_PI;
    }
}

double ddm_class::compute_g_factor(double dt, double zr, bool no_var_case)

{
    if (no_var_case)
    {
        // Case when sv = 0 (no variability)
        const double factor =
            std::exp(-m_a * zr * m_v - 0.5 * m_v2 * dt) / m_a2;
        return std::isfinite(factor) ? factor : 0.0;
    }
    else
    {
        // Case with variability (sv != 0)
        // const double zr2 = zr * zr;

        const double denominator = dt * m_sv2 + 1.0;
        const double exponent =
            -0.5 * (m_v2 * dt + 2 * m_v * m_a * zr - m_a2 * zr * zr * m_sv2) /
            denominator;
        const double factor =
            std::exp(exponent) / (m_a2 * std::sqrt(denominator));

        return std::isfinite(factor) ? factor : 0.0;
    }
}

std::pair<int, int> ddm_class::get_N(double dt, double ta2, double eps)
{
    int N_large = static_cast<int>(std::ceil(1.0 / (M_PI * std::sqrt(dt))));
    if (M_PI * ta2 * eps < 1.0)
    {
        double term0 = -2 * std::log(M_PI * ta2 * eps) / (M_PI_SQUARED * ta2);
        double term1 = std::ceil(std::sqrt(term0));
        N_large = std::max(N_large, static_cast<int>(term1));
    }

    int N_small;
    if (2.0 * std::sqrt(M_2PI * ta2) * eps < 1.0)
    {
        double term0 = std::log(2.0 * eps * std::sqrt(M_2PI * ta2));
        double term1 = 2.0 + std::sqrt(-2.0 * ta2 * term0);
        double term2 = std::sqrt(ta2) + 1;
        N_small = static_cast<int>(std::ceil(std::max(term2, term1)));
    }
    else
    {
        N_small = 2;
    }
    return {N_small, N_large};
}

double ddm_class::g_no_var(double dt, double zr)
{
    if (dt <= 0)
    {
        return 0.0;
    }

    double ta2 = dt / m_a2;
    const double factor = compute_g_factor(dt, zr, true);
    // Rcpp::Rcout << "dt, m_a2 " << dt << ", " << m_a2 << ", " << zr << "\n";
    // Rcpp::Rcout << "ta2, factor " << ta2 << ", " << factor << "\n";

    if (factor == 0)
    {
        return 0.0;
    }

    const double eps = EPSILON / factor;
    auto [N_small, N_large] = get_N(dt, ta2, eps);
    const bool use_small_time = (N_small < N_large);
    const int N = use_small_time ? N_small : N_large;
    return factor * compute_g_series(ta2, zr, use_small_time, N);
}

double ddm_class::integral_v(double dt, double zr)
{
    // zr could be m_zr or a range of zr based on m_zr and m_szr
    if (dt <= 0)
    {
        return 0;
    }
    if (m_sv == 0) // take m_sv from the member variable
    {
        // return g_no_var(dt, m_a, zr, m_v);
        return g_no_var(dt, zr);
    }

    const double ta2 = dt / m_a2;
    // Boolean to false
    // const double factor = compute_g_factor(dt, m_a, zr, m_v, false);
    const double factor = compute_g_factor(dt, zr, false);

    if (factor == 0)
    {
        return 0.0;
    }

    const double eps = EPSILON / factor;
    auto [N_small, N_large] = get_N(dt, ta2, eps);
    const bool use_small_time = (N_small < N_large);
    const int N = use_small_time ? N_small : N_large;

    return factor * compute_g_series(ta2, zr, use_small_time, N);
}

double ddm_class::integrate_v_over_zr(double dt)
{
    double lower = m_zr - 0.5 * m_szr;
    double upper = m_zr + 0.5 * m_szr;
    double width = upper - lower;
    double tmp_N = width / TUNE_INT_Z;

    int N = std::max(4, static_cast<int>(tmp_N));
    double step = width / N;

    double result = 0;
    for (double x = lower + 0.5 * step; x < upper; x += step)
    {
        // result += step * integral_v(dt, m_a, x, m_v);
        result += step * integral_v(dt, x); // changing zr value
    }

    return result / m_szr;
}

double ddm_class::integral_z(double dt)
{
    // integral over a uniform, fixed range of zr.
    // return (m_szr < TUNE_SZ_EPSILON) ? integral_v(dt, m_a, m_zr, m_v)
    //                                  : integrate_v_over_zr(dt);
    return (m_szr < TUNE_SZ_EPSILON) ? integral_v(dt, m_zr)
                                     : integrate_v_over_zr(dt);
}

double ddm_class::integrate_z_over_t(double dt)
{
    double lower = dt - 0.5 * m_st0;
    double upper = dt + 0.5 * m_st0;

    double width = upper - lower;
    double tmp_N = width / TUNE_INT_T0;
    int N = std::max(4, static_cast<int>(tmp_N));
    double step = width / N;

    double result = 0;

    for (double x = lower + 0.5 * step; x < upper; x += step)
    {
        result += step * integral_z(x);
    }

    return result / m_st0;
}

double ddm_class::integral_t0(double dt)
{
    // integral over a uniform, fixed range of t, depending on
    // integrate_z_over_t, and integral_z_g, which depends
    // on integrate_v_over_zr and integral_v_g.
    return (m_st0 < TUNE_ST0_EPSILON) ? integral_z(dt) : integrate_z_over_t(dt);
}

double ddm_class::g(double rt)
{
    // Rcpp::Rcout << rt << ", " << m_t_offset << "\n";

    return integral_t0(rt - m_t_offset);
}

std::vector<double> ddm_class::dddm(const std::vector<double> &rt)
{
    std::vector<double> out(rt.size());
    for (size_t i = 0; i < rt.size(); ++i)
    {
        // Rcpp::Rcout << rt[i] << "\n";
        out[i] = g(rt[i]);
    }

    return out;
}

} // namespace
} // namespace ddm