#pragma once
// DO NOT mix use RcppArmadillo.h and armadillo.
// Mixing them together cause memory segfault
#include <RcppArmadillo.h>

class ThetaInput
{
  public:
    unsigned int nmc, thin, nchain, nparameter;
    unsigned int report_length, max_init_attempts;
    bool is_print;                   // Whether to print progress
    std::vector<std::string> pnames; // no fixed parameters

    // Constructor with default values
    ThetaInput(unsigned int nmc = 100, unsigned int thin = 1,
               unsigned int nchain = 4, unsigned int nparameter = 5,
               unsigned int report_length = 10,
               unsigned int max_init_attempts = 1000, bool is_print = true,
               std::vector<std::string> pnames = std::vector<std::string>())
        : nmc(nmc), thin(thin), nchain(nchain), nparameter(nparameter),
          report_length(report_length), max_init_attempts(max_init_attempts),
          is_print(is_print), pnames(pnames)
    {
    }

    ~ThetaInput(){};
};
class SampleInput
{
  public:
    arma::cube theta;
    arma::mat summed_log_prior, log_likelihood;
    std::vector<std::string> pnames;

    SampleInput(const arma::cube &theta, const arma::mat &summed_log_prior,
                const arma::mat &log_likelihood,
                const std::vector<std::string> &pnames)
    {
        if (theta.n_rows == 0 || theta.n_cols == 0 || theta.n_slices == 0 ||
            summed_log_prior.n_rows == 0 || summed_log_prior.n_cols == 0 ||
            log_likelihood.n_rows == 0 || log_likelihood.n_cols == 0)
        {
            throw std::runtime_error(
                "SampleInput: Input cube/matrix has zero dimensions!");
        }

        this->theta = theta;
        this->summed_log_prior = summed_log_prior;
        this->log_likelihood = log_likelihood;
        this->pnames = pnames;
    }
};

class DEInput
{
  public:
    double pop_migration_prob, sub_migration_prob;
    double gamma_precursor, rp;
    bool is_hblocked, is_pblocked;
    unsigned int nparameter, nchain;
    bool pop_debug, sub_debug;

    // Constructor with default values
    DEInput(double pop_migration_prob = 0.00, double sub_migration_prob = 0.00,
            double gamma_precursor = 2.38, double rp = 0.001,
            bool is_hblocked = false, bool is_pblocked = false,
            unsigned int nparameter = 6, unsigned int nchain = 4,
            bool pop_debug = false, bool sub_debug = false)
        : pop_migration_prob(pop_migration_prob),
          sub_migration_prob(sub_migration_prob),
          gamma_precursor(gamma_precursor), rp(rp), is_hblocked(is_hblocked),
          is_pblocked(is_pblocked), nparameter(nparameter), nchain(nchain),
          pop_debug(pop_debug), sub_debug(sub_debug)
    {
    }

    ~DEInput(){};
};

inline DEInput new_DEInput(const Rcpp::S4 &de_input_r)
{
    // Depend on theta_helpers.h
    double pop_migration_prob = de_input_r.slot("pop_migration_prob");
    double sub_migration_prob = de_input_r.slot("sub_migration_prob");
    double gamma_precursor = de_input_r.slot("gamma_precursor");
    double rp = de_input_r.slot("rp");
    bool is_hblocked = de_input_r.slot("is_hblocked");
    bool is_pblocked = de_input_r.slot("is_pblocked");

    unsigned int nparameter = de_input_r.slot("nparameter");
    unsigned int nchain = de_input_r.slot("nchain");

    bool pop_debug = de_input_r.slot("pop_debug");
    bool sub_debug = de_input_r.slot("sub_debug");

    return DEInput(pop_migration_prob, sub_migration_prob, gamma_precursor, rp,
                   is_hblocked, is_pblocked, nparameter, nchain, pop_debug,
                   sub_debug);
}

inline ThetaInput new_ThetaInput(const Rcpp::S4 &theta_input_r)
{
    unsigned int nmc = theta_input_r.slot("nmc");
    unsigned int nchain = theta_input_r.slot("nchain");
    unsigned int thin = theta_input_r.slot("thin");
    unsigned int nparameter = theta_input_r.slot("nparameter");
    std::vector<std::string> pnames = theta_input_r.slot("pnames");
    unsigned int report_length = theta_input_r.slot("report_length");
    unsigned int max_init_attempts = theta_input_r.slot("max_init_attempts");
    bool is_report_progress = theta_input_r.slot("is_print");

    return ThetaInput(nmc, thin, nchain, nparameter, report_length,
                      max_init_attempts, is_report_progress, pnames);
}

inline SampleInput new_SampleInput(const Rcpp::S4 &samples)
{
    Rcpp::NumericVector theta_r = samples.slot("theta");
    Rcpp::NumericMatrix summed_log_prior_r = samples.slot("summed_log_prior");
    Rcpp::NumericMatrix log_likelihood_r = samples.slot("log_likelihoods");
    Rcpp::CharacterVector pnames_r = samples.slot("pnames");

    arma::cube theta = Rcpp::as<arma::cube>(theta_r);
    arma::mat summed_log_prior = Rcpp::as<arma::mat>(summed_log_prior_r);
    arma::mat log_likelihood = Rcpp::as<arma::mat>(log_likelihood_r);

    std::vector<std::string> pnames =
        Rcpp::as<std::vector<std::string>>(pnames_r);

    return SampleInput(theta, summed_log_prior, log_likelihood, pnames);
}
