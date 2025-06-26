#pragma once
#include "prior.h"
#include <RcppArmadillo.h>
///////////////////////////////////////////////////
/* Prior ---------------------------------------*/
///////////////////////////////////////////////////
std::shared_ptr<prior::prior_class> new_prior(const Rcpp::List &p_prior_r)
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

    return std::make_shared<prior::prior_class>(p0, p1, lower, upper, dist_type,
                                                log_p);
}
