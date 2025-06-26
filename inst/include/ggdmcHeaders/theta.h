#pragma once

#include "prior.h"
#include "theta_helpers.h"

class theta_phi
{
  protected:
    void yank_theta_input(const ThetaInput &theta_input)
    {
        m_nmc = theta_input.nmc;
        m_nchain = theta_input.nchain;
        m_thin = theta_input.thin;
        m_nparameter = theta_input.nparameter;
        m_report_length = theta_input.report_length;
        m_max_init_attempts = theta_input.max_init_attempts;
        m_is_print = theta_input.is_print;
    }
    void yank_sample_input(const SampleInput &sample_input)
    {
        m_previous_theta = sample_input.theta;
        m_previous_lp = sample_input.summed_log_prior;
        m_previous_ll = sample_input.log_likelihood;
        m_pnames = sample_input.pnames;
    }
    void fill_in_NAs()
    {

        if (m_previous_nmc == 1)
        {
            m_lp = m_previous_lp;
            m_ll = m_previous_ll;
            m_theta = m_previous_theta;
        }
        else
        {
            m_lp = arma::mat(m_nchain, m_nmc);
            m_ll = arma::mat(m_nchain, m_nmc);
            m_theta = arma::cube(m_nparameter, m_nchain, m_nmc);
            m_lp.fill(R_NegInf);
            m_ll.fill(R_NegInf);
            m_theta.fill(NA_REAL);
            m_theta.slice(0) = m_previous_theta.slice(m_previous_nmc - 1);
            m_lp.col(0) = m_previous_lp.col(m_previous_nmc - 1);
            m_ll.col(0) = m_previous_ll.col(m_previous_nmc - 1);
        }
    }

  public:
    unsigned int m_nmc, m_nchain, m_thin, m_nparameter;
    unsigned int m_report_length, m_max_init_attempts;
    bool m_is_print;

    /*-----------Theta, lp and ll used in DE-------------*/
    unsigned int m_start, m_store_i, m_nsample, m_previous_nmc;
    arma::cube m_theta, m_previous_theta;
    arma::mat m_lp, m_ll, m_used_theta, m_previous_lp, m_previous_ll;
    arma::vec m_used_lp, m_used_ll;
    std::vector<std::string> m_pnames;

    void store(unsigned int MC_index)
    {
        // Check if the current iteration should be stored based on thinning
        if (MC_index % m_thin == 0)
        {
            m_store_i++; // Increment the storage index
            if (m_store_i < m_nmc)
            {
                m_lp.col(m_store_i) = m_used_lp; // nchain x nmc
                m_ll.col(m_store_i) = m_used_ll;
                m_theta.slice(m_store_i) = m_used_theta; // npar x nchain x nmc
            }
        }
    }

    void print_progress(unsigned int MC_index)
    {
        if (MC_index % m_thin == 0)
        {

            if (m_is_print && (m_store_i + 1) % m_report_length == 0)
            {
                Rcpp::Rcout << m_store_i + 1 << " "
                            << std::flush; // Flush the buffer;
            }

            if (m_store_i > m_nmc)
            {
                std::ostringstream warning_msg;
                warning_msg
                    << "(m_store_i, m_nmc): (" << m_store_i << " "
                    << m_nmc + ") "
                    << "Warning: Exceeded storage dimensions at MC iteration "
                    << MC_index;

                Rf_warning("%s", warning_msg.str().c_str());
            }
        }
    }
    
};

class theta_class : public theta_phi
{
  public:
    /* Public theta_class functions */
    // Restart theta
    theta_class(const ThetaInput &theta_input, const SampleInput &sample_input)
    {
        yank_theta_input(theta_input);
        yank_sample_input(sample_input);

        // unsigned int previous_nmc;
        if (m_previous_theta.has_nan())
        {
            m_previous_nmc = 1;
        }
        else
        {
            m_previous_nmc = m_previous_theta.n_slices;
        }

        unsigned int previous_nchain = m_previous_theta.n_cols;

        if (m_nchain != previous_nchain)
        {
            throw std::runtime_error("Changing chain number is not supported.");
        }

        m_start = 1;
        m_store_i = m_start - 1;

        fill_in_NAs();

        ////////////////////////////////////////////////////
        m_nsample = 1 + (m_nmc - m_start) * m_thin;
        m_used_theta = m_theta.slice(m_store_i); // npar x nchain;
        m_used_lp = m_lp.col(m_store_i);         // nchains x 1
        m_used_ll = m_ll.col(m_store_i);         // nchains x 1
    }
};

using ThetaPtr = std::shared_ptr<theta_class>;

// class phi_class : public theta_phi
// {
//   public:
//     /* Phi-specific members */
//     std::vector<std::shared_ptr<theta_class>> m_subj_theta;
//     prior::PriorPtr m_p_prior, m_h_prior;
//     unsigned int m_nsubject;

//     /* Phi-specific constructor */
//     phi_class(const ThetaInput &theta_input, const SampleInput &sample_input,
//               std::vector<std::shared_ptr<theta_class>> subj_theta,
//               prior::PriorPtr p_prior, prior::PriorPtr h_prior)
//         : m_subj_theta(std::move(subj_theta)), m_p_prior(std::move(p_prior)),
//           m_h_prior(std::move(h_prior))
//     {
//         yank_theta_input(theta_input);
//         yank_sample_input(sample_input);
//         m_nsubject = m_subj_theta.size();

//         // unsigned int previous_nmc;
//         if (m_previous_theta.has_nan())
//         {
//             m_previous_nmc = 1;
//         }
//         else
//         {
//             m_previous_nmc = m_previous_theta.n_slices;
//         }

//         unsigned int previous_nchain = m_previous_theta.n_cols;

//         if (m_nchain != previous_nchain)
//         {
//             throw std::runtime_error("Changing chain number is not supported.");
//         }

//         m_start = 1;
//         m_store_i = m_start - 1;

//         fill_in_NAs();

//         ////////////////////////////////////////////////////
//         m_nsample = 1 + (m_nmc - m_start) * m_thin;
//         m_used_theta = m_theta.slice(m_store_i); // npar x nchain;
//         m_used_lp = m_lp.col(m_store_i);         // nchains x 1
//         m_used_ll = m_ll.col(m_store_i);         // nchains x 1
//     }

//     /* Phi-specific methods */
//     double sumloghlike(const arma::vec &parameters, size_t chain_idx)
//     {
//         m_p_prior->m_p0 = parameters.head(parameters.n_elem / 2);
//         m_p_prior->m_p1 = parameters.tail(parameters.n_elem / 2);

//         double out = 0;
//         for (size_t j = 0; j < m_nsubject; ++j)
//         {
//             out += m_p_prior->sumlogprior(
//                 m_subj_theta[j]->m_used_theta.col(chain_idx));
//         }
//         return out;
//     };
// };

// using PhiPtr = std::shared_ptr<phi_class>;