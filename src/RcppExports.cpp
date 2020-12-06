// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// soft_thresh_cpp
double soft_thresh_cpp(double z, double alpha, double tau);
RcppExport SEXP _regDIF_soft_thresh_cpp(SEXP zSEXP, SEXP alphaSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(soft_thresh_cpp(z, alpha, tau));
    return rcpp_result_gen;
END_RCPP
}
// firm_thresh_cpp
double firm_thresh_cpp(double z, double alpha, double tau, double gamma);
RcppExport SEXP _regDIF_firm_thresh_cpp(SEXP zSEXP, SEXP alphaSEXP, SEXP tauSEXP, SEXP gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(firm_thresh_cpp(z, alpha, tau, gamma));
    return rcpp_result_gen;
END_RCPP
}
// dnormCpp
arma::vec dnormCpp(arma::vec theta, double alpha_impact, double phi_impact, int num_quadpts);
RcppExport SEXP _regDIF_dnormCpp(SEXP thetaSEXP, SEXP alpha_impactSEXP, SEXP phi_impactSEXP, SEXP num_quadptsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_impact(alpha_impactSEXP);
    Rcpp::traits::input_parameter< double >::type phi_impact(phi_impactSEXP);
    Rcpp::traits::input_parameter< int >::type num_quadpts(num_quadptsSEXP);
    rcpp_result_gen = Rcpp::wrap(dnormCpp(theta, alpha_impact, phi_impact, num_quadpts));
    return rcpp_result_gen;
END_RCPP
}
// bernoulli_traceline_deriv
arma::mat bernoulli_traceline_deriv(arma::vec p_item, arma::mat theta, arma::mat predictors, int samp_size, int num_quadpts);
RcppExport SEXP _regDIF_bernoulli_traceline_deriv(SEXP p_itemSEXP, SEXP thetaSEXP, SEXP predictorsSEXP, SEXP samp_sizeSEXP, SEXP num_quadptsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type p_item(p_itemSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type predictors(predictorsSEXP);
    Rcpp::traits::input_parameter< int >::type samp_size(samp_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type num_quadpts(num_quadptsSEXP);
    rcpp_result_gen = Rcpp::wrap(bernoulli_traceline_deriv(p_item, theta, predictors, samp_size, num_quadpts));
    return rcpp_result_gen;
END_RCPP
}
// bernoulli_traceline_cpp
List bernoulli_traceline_cpp(arma::vec p_item, arma::mat theta, arma::mat predictors, int samp_size, int num_quadpts);
RcppExport SEXP _regDIF_bernoulli_traceline_cpp(SEXP p_itemSEXP, SEXP thetaSEXP, SEXP predictorsSEXP, SEXP samp_sizeSEXP, SEXP num_quadptsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type p_item(p_itemSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type predictors(predictorsSEXP);
    Rcpp::traits::input_parameter< int >::type samp_size(samp_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type num_quadpts(num_quadptsSEXP);
    rcpp_result_gen = Rcpp::wrap(bernoulli_traceline_cpp(p_item, theta, predictors, samp_size, num_quadpts));
    return rcpp_result_gen;
END_RCPP
}
// categorical_traceline_cpp
List categorical_traceline_cpp(arma::vec p_item, arma::mat theta, arma::mat predictors, int samp_size, int num_responses_item, int num_quadpts);
RcppExport SEXP _regDIF_categorical_traceline_cpp(SEXP p_itemSEXP, SEXP thetaSEXP, SEXP predictorsSEXP, SEXP samp_sizeSEXP, SEXP num_responses_itemSEXP, SEXP num_quadptsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type p_item(p_itemSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type predictors(predictorsSEXP);
    Rcpp::traits::input_parameter< int >::type samp_size(samp_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type num_responses_item(num_responses_itemSEXP);
    Rcpp::traits::input_parameter< int >::type num_quadpts(num_quadptsSEXP);
    rcpp_result_gen = Rcpp::wrap(categorical_traceline_cpp(p_item, theta, predictors, samp_size, num_responses_item, num_quadpts));
    return rcpp_result_gen;
END_RCPP
}
// cumulative_traceline_cpp
List cumulative_traceline_cpp(arma::vec p_item, arma::mat theta, arma::mat predictors, int samp_size, int num_responses_item, int num_quadpts);
RcppExport SEXP _regDIF_cumulative_traceline_cpp(SEXP p_itemSEXP, SEXP thetaSEXP, SEXP predictorsSEXP, SEXP samp_sizeSEXP, SEXP num_responses_itemSEXP, SEXP num_quadptsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type p_item(p_itemSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type predictors(predictorsSEXP);
    Rcpp::traits::input_parameter< int >::type samp_size(samp_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type num_responses_item(num_responses_itemSEXP);
    Rcpp::traits::input_parameter< int >::type num_quadpts(num_quadptsSEXP);
    rcpp_result_gen = Rcpp::wrap(cumulative_traceline_cpp(p_item, theta, predictors, samp_size, num_responses_item, num_quadpts));
    return rcpp_result_gen;
END_RCPP
}
// d_alpha_cpp
List d_alpha_cpp(arma::vec p_alpha, arma::vec p_phi, arma::mat etable_all, arma::mat theta, arma::mat mean_predictors, arma::mat var_predictors, int cov, int samp_size, int num_items, int num_quadpts);
RcppExport SEXP _regDIF_d_alpha_cpp(SEXP p_alphaSEXP, SEXP p_phiSEXP, SEXP etable_allSEXP, SEXP thetaSEXP, SEXP mean_predictorsSEXP, SEXP var_predictorsSEXP, SEXP covSEXP, SEXP samp_sizeSEXP, SEXP num_itemsSEXP, SEXP num_quadptsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type p_alpha(p_alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type p_phi(p_phiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type etable_all(etable_allSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mean_predictors(mean_predictorsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type var_predictors(var_predictorsSEXP);
    Rcpp::traits::input_parameter< int >::type cov(covSEXP);
    Rcpp::traits::input_parameter< int >::type samp_size(samp_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type num_items(num_itemsSEXP);
    Rcpp::traits::input_parameter< int >::type num_quadpts(num_quadptsSEXP);
    rcpp_result_gen = Rcpp::wrap(d_alpha_cpp(p_alpha, p_phi, etable_all, theta, mean_predictors, var_predictors, cov, samp_size, num_items, num_quadpts));
    return rcpp_result_gen;
END_RCPP
}
// d_phi_cpp
List d_phi_cpp(arma::vec p_alpha, arma::vec p_phi, arma::mat etable_all, arma::mat theta, arma::mat mean_predictors, arma::mat var_predictors, int cov, int samp_size, int num_items, int num_quadpts);
RcppExport SEXP _regDIF_d_phi_cpp(SEXP p_alphaSEXP, SEXP p_phiSEXP, SEXP etable_allSEXP, SEXP thetaSEXP, SEXP mean_predictorsSEXP, SEXP var_predictorsSEXP, SEXP covSEXP, SEXP samp_sizeSEXP, SEXP num_itemsSEXP, SEXP num_quadptsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type p_alpha(p_alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type p_phi(p_phiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type etable_all(etable_allSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mean_predictors(mean_predictorsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type var_predictors(var_predictorsSEXP);
    Rcpp::traits::input_parameter< int >::type cov(covSEXP);
    Rcpp::traits::input_parameter< int >::type samp_size(samp_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type num_items(num_itemsSEXP);
    Rcpp::traits::input_parameter< int >::type num_quadpts(num_quadptsSEXP);
    rcpp_result_gen = Rcpp::wrap(d_phi_cpp(p_alpha, p_phi, etable_all, theta, mean_predictors, var_predictors, cov, samp_size, num_items, num_quadpts));
    return rcpp_result_gen;
END_RCPP
}
// d_bernoulli_cpp
List d_bernoulli_cpp(std::string parm, arma::vec p_item, arma::mat etable1, arma::mat etable2, arma::mat theta, arma::mat predictors, int cov, int samp_size, int num_items, int num_quadpts);
RcppExport SEXP _regDIF_d_bernoulli_cpp(SEXP parmSEXP, SEXP p_itemSEXP, SEXP etable1SEXP, SEXP etable2SEXP, SEXP thetaSEXP, SEXP predictorsSEXP, SEXP covSEXP, SEXP samp_sizeSEXP, SEXP num_itemsSEXP, SEXP num_quadptsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type parm(parmSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type p_item(p_itemSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type etable1(etable1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type etable2(etable2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type predictors(predictorsSEXP);
    Rcpp::traits::input_parameter< int >::type cov(covSEXP);
    Rcpp::traits::input_parameter< int >::type samp_size(samp_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type num_items(num_itemsSEXP);
    Rcpp::traits::input_parameter< int >::type num_quadpts(num_quadptsSEXP);
    rcpp_result_gen = Rcpp::wrap(d_bernoulli_cpp(parm, p_item, etable1, etable2, theta, predictors, cov, samp_size, num_items, num_quadpts));
    return rcpp_result_gen;
END_RCPP
}
// d_categorical_cpp
List d_categorical_cpp(std::string parm, arma::vec p_item, Rcpp::List etable, arma::mat theta, arma::mat predictors, int thr, int cov, int samp_size, int num_responses_item, int num_items, int num_quadpts);
RcppExport SEXP _regDIF_d_categorical_cpp(SEXP parmSEXP, SEXP p_itemSEXP, SEXP etableSEXP, SEXP thetaSEXP, SEXP predictorsSEXP, SEXP thrSEXP, SEXP covSEXP, SEXP samp_sizeSEXP, SEXP num_responses_itemSEXP, SEXP num_itemsSEXP, SEXP num_quadptsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type parm(parmSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type p_item(p_itemSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type etable(etableSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type predictors(predictorsSEXP);
    Rcpp::traits::input_parameter< int >::type thr(thrSEXP);
    Rcpp::traits::input_parameter< int >::type cov(covSEXP);
    Rcpp::traits::input_parameter< int >::type samp_size(samp_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type num_responses_item(num_responses_itemSEXP);
    Rcpp::traits::input_parameter< int >::type num_items(num_itemsSEXP);
    Rcpp::traits::input_parameter< int >::type num_quadpts(num_quadptsSEXP);
    rcpp_result_gen = Rcpp::wrap(d_categorical_cpp(parm, p_item, etable, theta, predictors, thr, cov, samp_size, num_responses_item, num_items, num_quadpts));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_regDIF_soft_thresh_cpp", (DL_FUNC) &_regDIF_soft_thresh_cpp, 3},
    {"_regDIF_firm_thresh_cpp", (DL_FUNC) &_regDIF_firm_thresh_cpp, 4},
    {"_regDIF_dnormCpp", (DL_FUNC) &_regDIF_dnormCpp, 4},
    {"_regDIF_bernoulli_traceline_deriv", (DL_FUNC) &_regDIF_bernoulli_traceline_deriv, 5},
    {"_regDIF_bernoulli_traceline_cpp", (DL_FUNC) &_regDIF_bernoulli_traceline_cpp, 5},
    {"_regDIF_categorical_traceline_cpp", (DL_FUNC) &_regDIF_categorical_traceline_cpp, 6},
    {"_regDIF_cumulative_traceline_cpp", (DL_FUNC) &_regDIF_cumulative_traceline_cpp, 6},
    {"_regDIF_d_alpha_cpp", (DL_FUNC) &_regDIF_d_alpha_cpp, 10},
    {"_regDIF_d_phi_cpp", (DL_FUNC) &_regDIF_d_phi_cpp, 10},
    {"_regDIF_d_bernoulli_cpp", (DL_FUNC) &_regDIF_d_bernoulli_cpp, 10},
    {"_regDIF_d_categorical_cpp", (DL_FUNC) &_regDIF_d_categorical_cpp, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_regDIF(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
