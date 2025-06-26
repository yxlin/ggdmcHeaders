#pragma once
#include <RcppArmadillo.h>

///////////////////////////////////////////////////
/* ------------- Design helpers  ------------- */
///////////////////////////////////////////////////

std::shared_ptr<design::design_class> new_design_light(const Rcpp::S4 &dmi)
{
    // 1. The user enterd elements
    Rcpp::S4 model_r = dmi.slot("model");
    Rcpp::List parameter_map_r = model_r.slot("parameter_map");
    Rcpp::CharacterVector accumulators_r = model_r.slot("accumulators");
    Rcpp::List factors_r = model_r.slot("factors");
    Rcpp::List match_map_r = model_r.slot("match_map");
    Rcpp::List constants_r = model_r.slot("constants");
    std::string model_str = model_r.slot("type");

    // Convert to C++ types
    std::map<std::string, std::vector<std::string>> parameter_map =
        list_to_map<std::string>(parameter_map_r);

    std::vector<std::string> accumulators =
        Rcpp::as<std::vector<std::string>>(accumulators_r);

    std::map<std::string, std::vector<std::string>> factors =
        list_to_map<std::string>(factors_r);

    std::map<std::string, std::map<std::string, std::string>> match_map =
        nested_list_to_map(match_map_r);

    std::map<std::string, double> constants = constants_to_map(constants_r);

    // Light additional arguments
    auto model_boolean = R_ucube_to_std_ucube(model_r.slot("model_boolean"));

    auto cell_names =
        Rcpp::as<std::vector<std::string>>(model_r.slot("cell_names"));
    auto parameter_x_condition_names = Rcpp::as<std::vector<std::string>>(
        model_r.slot("parameter_x_condition_names"));

    auto out = std::make_shared<design::design_class>(
        parameter_map, accumulators, cell_names, parameter_x_condition_names,
        constants, model_boolean, model_str);
    // if (model_str == "lba" && dmi.hasSlot("node_1_index"))
    // {
    //     // auto node_1_index = r_mat_to_std_mat<Rcpp::IntegerMatrix, unsigned
    //     // int>(
    //     //     dmi.slot("node_1_index"));

    //     out = std::make_shared<design::design_class>(
    //         parameter_map, accumulators, cell_names,
    //         parameter_x_condition_names, constants, model_boolean,
    //         model_str);
    //     // node_1_index);
    //     // Rcpp::Rcout << "lba new_design_light type = " << model_str <<
    //     "\n";
    // }
    // else if (model_str == "fastdm")
    // {
    //     out = std::make_shared<design::design_class>(
    //         parameter_map, accumulators, cell_names,
    //         parameter_x_condition_names, constants, model_boolean,
    //         model_str);
    //     // Rcpp::Rcout << "fastdm new_design_light type = " << model_str <<
    //     // "\n";
    // }
    // else if (model_str == "hyper")
    // {

    //     out = std::make_shared<design::design_class>(
    //         parameter_map, accumulators, cell_names,
    //         parameter_x_condition_names, constants, model_boolean,
    //         model_str);
    //     // Rcpp::Rcout << "hyper new_design_light Model type = " << model_str
    //     //             << "\n";
    // }
    // else
    // {
    //     throw std::runtime_error("Undefined model type = ");
    // }

    return out;
}

///////////////////////////////////////////////////
/* ------------- Likelihood  ------------- */
///////////////////////////////////////////////////
std::vector<std::shared_ptr<likelihood::likelihood_class>>
new_likelihoods(const Rcpp::List &dmis)
{
    unsigned int n_subject = dmis.size();
    std::vector<std::shared_ptr<likelihood::likelihood_class>> likelihood_ptr(
        n_subject);

    for (size_t subject_idx = 0; subject_idx < n_subject; ++subject_idx)
    {
        Rcpp::S4 dmi = dmis[subject_idx];
        Rcpp::S4 model_r = dmi.slot("model");
        Rcpp::List data_r = dmi.slot("data");

        auto d_ptr = new_design_light(dmi);
        std::string model_str = model_r.slot("type");

        // Rcpp::List data_r = data_list_r[subject_idx];
        std::vector<std::string> cell_names = data_r.names();
        size_t ncell = cell_names.size();

        std::vector<std::vector<double>> data_cpp(ncell);

        for (size_t cell_idx = 0; cell_idx < ncell; ++cell_idx)
        {
            // Must assign in two steps?
            std::vector<double> rt = data_r[cell_idx];
            data_cpp[cell_idx] = rt;
        }

        std::shared_ptr<likelihood::likelihood_class> out;
        // if (model_str == "lba" && dmi.hasSlot("is_positive_drift"))
        // {
        auto is_positive_drift =
            Rcpp::as<std::vector<bool>>(dmi.slot("is_positive_drift"));

        likelihood_ptr[subject_idx] =
            std::make_shared<likelihood::likelihood_class>(
                d_ptr, data_cpp, cell_names, model_str, is_positive_drift);
        // Rcpp::Rcout << "lba new_likelihood type = " << model_str << "\n";
        // }
        // else if (model_str == "fastdm")
        // {
        //     likelihood_ptr[subject_idx] =
        //         std::make_shared<likelihood::likelihood_class>(
        //             d_ptr, data_cpp, cell_names, model_str);
        //     Rcpp::Rcout << "ddm new_likelihoods type = " << model_str <<
        //     "\n";
        // }
        // else
        // {
        //     Rcpp::Rcout << "Undefined new_likelihoods type = " << model_str
        //                 << "\n";
        // }

        // likelihood_ptr[subject_idx] =
        //     std::make_shared<likelihood::likelihood_class>(
        //         d_ptr, data_cpp, cell_names, model_str, is_positive_drift);
    }
    return likelihood_ptr;
}

/**
 * Creates a likelihood object from DMI (Data Model Interface)
 *
 * @param dmi An S4 object containing model and data specifications
 * @param p_prior Optional shared pointer to prior object (default nullptr)
 * @return Shared pointer to likelihood_class object
 */
std::shared_ptr<likelihood::likelihood_class>
new_likelihood(const Rcpp::S4 &dmi,
               const std::shared_ptr<prior::prior_class> &p_prior = nullptr)
{
    if (!dmi.hasSlot("model"))
    {
        Rcpp::stop("DMI must have the slot: model");
    }

    Rcpp::S4 model_r = dmi.slot("model");
    std::string model_str = model_r.slot("type");
    auto d_ptr = new_design_light(dmi);

    if (p_prior)
    {
        // Hyper-parameter version (second original function)
        arma::mat theta_data = Rcpp::as<arma::mat>(dmi.slot("data"));
        return std::make_shared<likelihood::likelihood_class>(
            d_ptr, p_prior, theta_data, model_str);
    }
    else
    {
        // Rcpp::Rcout << "new_likelihood type = " << model_str << "\n";

        if (!dmi.hasSlot("is_positive_drift"))
        {
            Rcpp::stop("DMI missing required slot: is_positive_drift");
        }

        Rcpp::List data_r = dmi.slot("data");
        std::vector<std::string> cell_names = data_r.names();
        size_t ncell = cell_names.size();
        std::vector<std::vector<double>> data_cpp(ncell);

        for (size_t cell_idx = 0; cell_idx < ncell; ++cell_idx)
        {
            // Must assign in two steps?
            std::vector<double> rt = data_r[cell_idx];
            data_cpp[cell_idx] = rt;
        }

        std::shared_ptr<likelihood::likelihood_class> out;
        auto is_positive_drift =
            Rcpp::as<std::vector<bool>>(dmi.slot("is_positive_drift"));
        return std::make_shared<likelihood::likelihood_class>(
            d_ptr, data_cpp, cell_names, model_str, is_positive_drift);
    }
}
