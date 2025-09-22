#pragma once
#include <RcppArmadillo.h>

namespace cdm
{
enum class Rule
{
    DINA,
    DINO
};

enum class CDMParamLayout
{
    BlockGuessSlip, // [g1 ... gN, s1 ... sN]
    Interleaved     // [g1, s1, g2, s2, ...]
};

class cdm_class
{
  private:
    // --- Model spec ---
    arma::imat m_Q;    // J x K (0/1)
    arma::vec m_guess; // J
    arma::vec m_slip;  // J
    Rule m_rule = Rule::DINA;

    // --- Latent generator (optional) ---
    arma::vec m_mean;  // K
    arma::mat m_Sigma; // K x K

    // --- Optional fixed mastery A (N x K, 0/1) ---
    bool m_has_mvn;
    bool m_has_alpha;

    // --- Derived / caches ---
    arma::imat m_profiles; // A: L x K
    arma::imat m_eta;      // L x J
    // std::vector<double> m_pi; // L (prior over classes)
    arma::vec m_pi;

    std::vector<double> m_prior_pi;

    // ---------- helpers ----------
    void build_profiles()
    {
        // Build all 2^K attribute profiles A (L x K, rows are 0/1 vectors)
        // L = 2^K; row l uses the bits of (l) in binary (0-based).
        if (m_nskill <= 0)
            Rcpp::stop("K (m_nskill) must be positive.");

        const unsigned int MAX_SHIFT =
            static_cast<unsigned int>(8 * sizeof(std::size_t) - 1);

        if (m_nskill > MAX_SHIFT)
            Rcpp::stop("K=%u too large for bitshift.", m_nskill);

        if (m_nskill > 25) // your heuristic warning
            Rcpp::warning("K=%u -> 2^K is large.", m_nskill);

        m_profiles.set_size(m_nprofile, m_nskill);
        m_profiles.zeros();

        // Row l encodes K-bit mastery vector (little-endian: skill k at bit k)
        for (std::size_t l = 0; l < m_nprofile; ++l)
        {
            for (std::size_t k = 0; k < m_nskill; ++k)
            {
                m_profiles(static_cast<arma::uword>(l),
                           static_cast<arma::uword>(k)) =
                    static_cast<int>((l >> k) & 1u);
            }
        }
    }

    arma::vec row_logsumexp_weighted_(const arma::mat &X)
    {

        arma::vec w = arma::conv_to<arma::vec>::from(m_prior_pi);

        if ((int)X.n_cols != (int)w.n_elem)
            Rcpp::stop("row_logsumexp_weighted: dim mismatch.");
        arma::vec m = arma::max(X, 1);  // N x 1
        arma::mat Z = X.each_col() - m; // N x L
        arma::mat EZ = arma::exp(Z);    // N x L
        arma::vec s = EZ * w;           // N
        return m + arma::log(s);
    }

    void get_eta()
    {
        // A: L x K mastery profiles; Q: J x K requirement matrix (0/1)
        arma::imat A = arma::conv_to<arma::imat>::from(m_profiles);
        arma::imat Q = arma::conv_to<arma::imat>::from(m_Q);

        // AQ: L x J = A (LxK) * Q^T (KxJ)
        arma::imat AQ = arma::conv_to<arma::imat>::from(
            arma::conv_to<arma::mat>::from(A) *
            arma::conv_to<arma::mat>::from(Q).t());

        m_eta.set_size(AQ.n_rows, AQ.n_cols); // L x J

        if (m_rule == Rule::DINA)
        {
            // Need: J x 1 (number of required attributes per item)
            arma::ivec need = arma::sum(Q, 1); // Jx1
            // Broadcast need^T to L x J and compare
            m_eta = arma::conv_to<arma::imat>::from(
                AQ == arma::repmat(need.t(), AQ.n_rows, 1));
        }
        else
        { // DINO
            // Rcpp::Rcout << "DINO rule in get_eta\n";
            m_eta = arma::conv_to<arma::imat>::from(AQ >= 1);
        }
    }

    arma::mat compute_eta_probability(double eps = 1e-12) const
    {
        // p = eta*(1-slip) + (1-eta)*guess
        arma::mat eta = arma::conv_to<arma::mat>::from(m_eta); // L x J
        arma::rowvec one_minus_slip = (1.0 - m_slip).t();      // 1 x J
        arma::rowvec guess_row = m_guess.t();                  // 1 x J

        arma::mat term1 = eta;
        term1.each_row() %= one_minus_slip;
        arma::mat term2 = 1.0 - eta;
        term2.each_row() %= guess_row;

        return arma::clamp(term1 + term2, eps, 1.0 - eps); // L x J
    }
    static bool
    extract_acc_column(const std::vector<std::vector<double>> &parameters,
                       std::size_t acc_idx, std::vector<double> &out,
                       std::string &err)
    {
        out.clear();
        out.reserve(parameters.size());
        for (std::size_t r = 0; r < parameters.size(); ++r)
        {
            if (parameters[r].size() <= acc_idx)
            {
                err = "parameters[" + std::to_string(r) + "] has only " +
                      std::to_string(parameters[r].size()) +
                      " accumulators, need index " + std::to_string(acc_idx);
                return false;
            }
            out.push_back(parameters[r][acc_idx]);
        }
        return true;
    }

    void map_guess_slip(const std::vector<double> &acc_vec,
                        CDMParamLayout layout)
    {

        if (layout == CDMParamLayout::BlockGuessSlip)
        {
            // acc_vec = [g1..gN, s1..sN]
            std::copy_n(acc_vec.begin(), m_nitem, m_guess.begin());
            std::copy_n(acc_vec.begin() + m_nitem, m_nitem, m_slip.begin());
        }
        else
        {
            // acc_vec = [g1, s1, g2, s2, ...]
            for (std::size_t i = 0; i < m_nitem; ++i)
            {
                m_guess[i] = acc_vec[2 * i];
                m_slip[i] = acc_vec[2 * i + 1];
            }
        }
    }

    // --- helper: validate & convert Q (std::vector<vector<double>>) ->
    // arma::imat
    inline arma::imat
    q_std_to_arma_imat(const std::vector<std::vector<double>> &Q,
                       std::string *err = nullptr)
    {
        if (Q.empty() || Q[0].empty())
        {
            if (err)
                *err = "Q dims invalid (got 0 x 0).";
            return arma::imat(); // empty
        }
        const std::size_t n_rows = Q.size();
        const std::size_t n_cols = Q[0].size();

        // rectangular check
        for (std::size_t i = 1; i < n_rows; ++i)
        {
            if (Q[i].size() != n_cols)
            {
                if (err)
                    *err = "Q is not rectangular: row " + std::to_string(i) +
                           " has size " + std::to_string(Q[i].size()) +
                           " vs first row " + std::to_string(n_cols) + ".";
                return arma::imat();
            }
        }

        arma::imat out(n_rows, n_cols, arma::fill::zeros);
        // accept near-0/near-1 values and coerce to 0/1
        for (std::size_t i = 0; i < n_rows; ++i)
        {
            for (std::size_t j = 0; j < n_cols; ++j)
            {
                const double x = Q[i][j];
                int v;
                if (std::fabs(x) < 1e-12)
                    v = 0;
                else if (std::fabs(x - 1.0) < 1e-12)
                    v = 1;
                else
                {
                    if (err)
                        *err = "Q must be binary; got Q[" + std::to_string(i) +
                               "," + std::to_string(j) +
                               "]=" + std::to_string(x) + ".";
                    return arma::imat();
                }
                out(i, j) = v;
            }
        }
        return out;
    }

    // ---------- simulation helpers ----------
    arma::mat sample_mvn(unsigned int n_student)
    {
        // const unsigned int K = mean.n_elem;
        arma::mat L = arma::chol(m_Sigma, "lower"); // throws if not SPD
        arma::mat Z = arma::randn<arma::mat>(n_student, m_nskill); // N x K
        arma::mat X = Z * L.t();                                   // N x K
        X.each_row() += m_mean.t();
        return X;
    }

    static arma::Mat<unsigned char> thresh_(const arma::mat &X)
    {
        arma::umat M = (X > 0.0);
        return arma::conv_to<arma::Mat<unsigned char>>::from(M);
    }

  public:
    arma::imat m_simulation_output;
    arma::Mat<unsigned char> m_alpha; // N x K (0/1) as bytes
    arma::mat m_X;
    unsigned int m_nitem, m_nskill, m_nprofile;

    cdm_class(const std::vector<double> &parameters = {0.21, 0.21, 0.21, 0.21,
                                                       0.21, 0.15, 0.15, 0.15,
                                                       0.15, 0.15},
              const std::vector<std::vector<double>> &Q = {
                  {1.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}, {1.0, 0.0}, {0.0, 1.0}})
    {
        m_nitem = Q.size();
        m_nskill = (Q.empty() ? 0u : static_cast<unsigned>(Q[0].size()));
        m_nprofile =
            (m_nskill >= (sizeof(unsigned) * 8) ? 0u : (1u << m_nskill));

        m_prior_pi = std::vector<double>(
            m_nprofile,
            1.0 / (m_nprofile ? static_cast<double>(m_nprofile) : 1.0));
    }

    ~cdm_class() {};

    // ----- Setters -----
    void set_parameters(const std::vector<std::vector<double>> &parameters,
                        const arma::mat &Q,
                        CDMParamLayout layout = CDMParamLayout::BlockGuessSlip,
                        std::size_t acc_idx = 0, bool debug = false)
    {
        // set_default_parameters(Q);
        build_profiles();

        // 2) Extract one accumulator column (usually 0)
        std::string err;
        std::vector<double> parameter_vec;
        if (!extract_acc_column(parameters, acc_idx, parameter_vec, err))
        {
            Rcpp::stop("set_parameters(): %s", err.c_str());
        }

        // 3) Validate *shape* vs expected (no value checks here)
        const std::size_t expected = (layout == CDMParamLayout::BlockGuessSlip)
                                         ? 2 * m_nitem
                                         : 2 * m_nitem;

        if (parameter_vec.size() != expected)
        {
            Rcpp::stop(
                "Parameter length mismatch: got %zu for accumulator %zu, "
                "expected %zu (= 2*nitem=%zu) with layout=%s.",
                parameter_vec.size(), acc_idx, expected, m_nitem,
                (layout == CDMParamLayout::BlockGuessSlip ? "BlockGuessSlip"
                                                          : "Interleaved"));
        }

        // 4) Map to internal guess/slip
        map_guess_slip(parameter_vec, layout);
        get_eta();

        if (debug)
        {
            m_guess.t().print("guess (set_parameters)");
            m_slip.t().print("slip (set_parameters)");
        }
        // TODO: store m_X and m_alpha here
    }

    bool validate_parameters(bool debug = false) const
    {
        // const std::size_t nitem = static_cast<std::size_t>(m_Q.n_rows);

        if (m_guess.n_elem != m_nitem || m_slip.n_elem != m_nitem)
        {
            if (debug)
            {
                Rcpp::Rcout << "Size mismatch: guess/slip vs Q rows. "
                            << "guess=" << m_guess.n_elem
                            << ", slip=" << m_slip.n_elem
                            << ", Q.n_rows=" << m_Q.n_rows << "\n";
            }
            return false;
        }

        for (std::size_t j = 0; j < m_nitem; ++j)
        {
            const double g = m_guess[j];
            const double s = m_slip[j];
            if (!(g >= 0.0 && g <= 1.0 && s >= 0.0 && s <= 1.0))
            {
                if (debug)
                {
                    Rcpp::Rcout << "guess/slip out of [0,1] at item " << j
                                << " (g=" << g << ", s=" << s << ")\n";
                }
                return false;
            }
        }
        return true;
    }

    void print_parameters(const std::string &cell_name = "") const
    {
        if (!cell_name.empty())
            Rcpp::Rcout << "[CDM] Cell: " << cell_name << "\n";
        m_guess.t().print("guess");
        m_slip.t().print("slip");
        // Rcpp::Rcout << std::endl;
    }

    void set_rule(const std::string &rule)
    {
        if (rule == "DINA")
            m_rule = Rule::DINA;
        else if (rule == "DINO")
            m_rule = Rule::DINO;
        else
            Rcpp::stop("Unknown rule. Use 'DINA' or 'DINO'.");
    }

    std::vector<double> dcdm(const std::vector<double> &rt)
    {
        // std::vector<double> out(rt.size());
        // arma::mat p = dina_prob_(); // L x J
        arma::mat p = compute_eta_probability();

        // Prepare log p and log(1-p) (J x L) by transposing later
        arma::mat logp = arma::log(p);         // L x J
        arma::mat log1mp = arma::log(1.0 - p); // L x J

        std::size_t N = rt.size() / m_nitem;
        arma::mat W(N, m_nitem, arma::fill::ones);
        arma::mat Y(N, m_nitem, arma::fill::zeros);

        for (std::size_t i = 0; i < N; ++i)
        {
            for (std::size_t j = 0; j < m_nitem; ++j)
            {
                double v = rt[i * m_nitem + j]; // flatten: row-major index

                if (std::isnan(v)) // check missing
                {
                    W(i, j) = 0.0;
                    Y(i, j) = 0.0; // ignored
                }
                else
                {
                    if (v != 0.0 && v != 1.0)
                        Rcpp::stop("Y must be 0/1/NA.");
                    Y(i, j) = v;
                }
            }
        }

        arma::mat term0 = (Y % W) * logp.t();
        arma::mat term1 = ((1.0 - Y) % W) * log1mp.t();
        arma::mat li_mat = term0 + term1; // N x L

        arma::vec ll_i = row_logsumexp_weighted_(li_mat); // length N
        arma::vec likelihood_i = arma::exp(ll_i);
        std::vector<double> out(likelihood_i.begin(), likelihood_i.end());

        return out;
    }

    // ----- Simulation ---------
    // void set_mvn(const arma::vec &mean, const arma::mat &Sigma)
    // {
    //     if ((int)Sigma.n_rows != (int)Sigma.n_cols ||
    //         (int)Sigma.n_rows != (int)mean.n_elem)
    //         Rcpp::stop("Sigma must be KxK and match mean length.");
    //     m_mean = mean;
    //     m_Sigma = Sigma;
    //     m_has_mvn = true;
    // }

    void set_default_parameters(const arma::mat &Q)
    {
        // 1) Structural checks on Q
        if (Q.n_rows <= 0 || Q.n_cols <= 0)
        {
            Rcpp::stop("Q dims invalid (got %lu x %lu).",
                       static_cast<unsigned long>(Q.n_rows),
                       static_cast<unsigned long>(Q.n_cols));
        }
        m_Q = arma::conv_to<arma::imat>::from(Q);
        m_nitem = Q.n_rows;
        m_nskill = Q.n_cols;

        m_nprofile =
            (m_nskill >= (sizeof(unsigned) * 8) ? 0u : (1u << m_nskill));

        m_prior_pi = std::vector<double>(
            m_nprofile,
            1.0 / (m_nprofile ? static_cast<double>(m_nprofile) : 1.0));

        // mean = zero vector of length K
        m_mean = arma::zeros<arma::vec>(m_nskill);

        // Sigma = identity matrix of size K
        m_Sigma = arma::eye<arma::mat>(m_nskill, m_nskill);
        m_has_mvn = false;
        m_has_alpha = false;

        m_guess.set_size(m_nitem);
        m_slip.set_size(m_nitem);
    }

    void rcdm(unsigned int n_student)
    {
        if (!m_has_alpha)
        {
            if (!m_has_mvn)
            {
                m_X = sample_mvn(n_student);
            }

            m_alpha = thresh_(m_X); // N x K
            m_has_mvn = true;
            m_has_alpha = true;
        }
        else
        {
            if (m_alpha.n_rows != n_student)
                Rcpp::stop("alpha rows != n_student.");
        }

        m_simulation_output.set_size(n_student, m_nitem);
        m_simulation_output.zeros();

        // Precompute "necessary attributes" per item depending on rule
        arma::ivec need(m_nitem);
        if (m_rule == Rule::DINA)
        {
            need = arma::sum(m_Q, 1); // rowSums(Q)
        }
        else
        {
            // Rcpp::Rcout << "DINO rule (rcdm)\n";
            need.fill(1);
        }

        for (size_t j = 0; j < m_nitem; ++j)
        {
            double p_master = 1.0 - m_slip(j);
            double p_non = m_guess(j);
            int req = need(j);

            for (size_t i = 0; i < n_student; ++i)
            {
                int cnt = 0;
                for (size_t k = 0; k < m_nskill; ++k)
                {
                    cnt += (m_alpha(i, k) && m_Q(j, k));
                }

                double p = (cnt >= req) ? p_master : p_non;
                m_simulation_output(i, j) = (int)R::rbinom(1.0, p);
            }
        }
    }
};

} // namespace cdm
