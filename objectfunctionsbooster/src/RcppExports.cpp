// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// sim_vsveipdr_Gs_booster
arma::mat sim_vsveipdr_Gs_booster(int t_max, arma::vec y, arma::vec alpha, arma::vec beta, arma::vec gamma_d, arma::vec gamma_r, arma::vec thetas, int M, arma::vec dG1s, arma::vec dG2s, arma::vec dG3s, double mu_1, double mu_2, double mu_3, double kappa, double rho, double omega, int r_beta);
RcppExport SEXP _objectfunctionsbooster_sim_vsveipdr_Gs_booster(SEXP t_maxSEXP, SEXP ySEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP gamma_dSEXP, SEXP gamma_rSEXP, SEXP thetasSEXP, SEXP MSEXP, SEXP dG1sSEXP, SEXP dG2sSEXP, SEXP dG3sSEXP, SEXP mu_1SEXP, SEXP mu_2SEXP, SEXP mu_3SEXP, SEXP kappaSEXP, SEXP rhoSEXP, SEXP omegaSEXP, SEXP r_betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type t_max(t_maxSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma_d(gamma_dSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma_r(gamma_rSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type thetas(thetasSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type dG1s(dG1sSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type dG2s(dG2sSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type dG3s(dG3sSEXP);
    Rcpp::traits::input_parameter< double >::type mu_1(mu_1SEXP);
    Rcpp::traits::input_parameter< double >::type mu_2(mu_2SEXP);
    Rcpp::traits::input_parameter< double >::type mu_3(mu_3SEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< int >::type r_beta(r_betaSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_vsveipdr_Gs_booster(t_max, y, alpha, beta, gamma_d, gamma_r, thetas, M, dG1s, dG2s, dG3s, mu_1, mu_2, mu_3, kappa, rho, omega, r_beta));
    return rcpp_result_gen;
END_RCPP
}
// object_kappa_tra_booster
arma::rowvec object_kappa_tra_booster(arma::vec x, arma::mat dat_sim, int t_V, int iter, arma::vec thetas, int M, double mu_1, double mu_2, double mu_3, int r_beta, arma::mat bspline2);
RcppExport SEXP _objectfunctionsbooster_object_kappa_tra_booster(SEXP xSEXP, SEXP dat_simSEXP, SEXP t_VSEXP, SEXP iterSEXP, SEXP thetasSEXP, SEXP MSEXP, SEXP mu_1SEXP, SEXP mu_2SEXP, SEXP mu_3SEXP, SEXP r_betaSEXP, SEXP bspline2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type dat_sim(dat_simSEXP);
    Rcpp::traits::input_parameter< int >::type t_V(t_VSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type thetas(thetasSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< double >::type mu_1(mu_1SEXP);
    Rcpp::traits::input_parameter< double >::type mu_2(mu_2SEXP);
    Rcpp::traits::input_parameter< double >::type mu_3(mu_3SEXP);
    Rcpp::traits::input_parameter< int >::type r_beta(r_betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type bspline2(bspline2SEXP);
    rcpp_result_gen = Rcpp::wrap(object_kappa_tra_booster(x, dat_sim, t_V, iter, thetas, M, mu_1, mu_2, mu_3, r_beta, bspline2));
    return rcpp_result_gen;
END_RCPP
}
// object_kappa_booster
arma::rowvec object_kappa_booster(arma::vec x, arma::mat dat_sim, int t_V, int iter, arma::vec thetas, int M, double mu_1, double mu_2, double mu_3, int r_beta, arma::mat bspline2, int inter_b_min, int inter_b_max, int inter_l);
RcppExport SEXP _objectfunctionsbooster_object_kappa_booster(SEXP xSEXP, SEXP dat_simSEXP, SEXP t_VSEXP, SEXP iterSEXP, SEXP thetasSEXP, SEXP MSEXP, SEXP mu_1SEXP, SEXP mu_2SEXP, SEXP mu_3SEXP, SEXP r_betaSEXP, SEXP bspline2SEXP, SEXP inter_b_minSEXP, SEXP inter_b_maxSEXP, SEXP inter_lSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type dat_sim(dat_simSEXP);
    Rcpp::traits::input_parameter< int >::type t_V(t_VSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type thetas(thetasSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< double >::type mu_1(mu_1SEXP);
    Rcpp::traits::input_parameter< double >::type mu_2(mu_2SEXP);
    Rcpp::traits::input_parameter< double >::type mu_3(mu_3SEXP);
    Rcpp::traits::input_parameter< int >::type r_beta(r_betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type bspline2(bspline2SEXP);
    Rcpp::traits::input_parameter< int >::type inter_b_min(inter_b_minSEXP);
    Rcpp::traits::input_parameter< int >::type inter_b_max(inter_b_maxSEXP);
    Rcpp::traits::input_parameter< int >::type inter_l(inter_lSEXP);
    rcpp_result_gen = Rcpp::wrap(object_kappa_booster(x, dat_sim, t_V, iter, thetas, M, mu_1, mu_2, mu_3, r_beta, bspline2, inter_b_min, inter_b_max, inter_l));
    return rcpp_result_gen;
END_RCPP
}
// sim_vsveipdr_Gs_booster_reinfect
arma::mat sim_vsveipdr_Gs_booster_reinfect(int t_max, arma::vec y, arma::vec alpha, arma::vec beta, arma::vec gamma_d, arma::vec gamma_r, arma::vec thetas, int M, arma::vec dG1s, arma::vec dG2s, arma::vec dG3s, double mu_1, double mu_2, double mu_3, double mu_r, double kappa, double rho, double omega, int r_beta);
RcppExport SEXP _objectfunctionsbooster_sim_vsveipdr_Gs_booster_reinfect(SEXP t_maxSEXP, SEXP ySEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP gamma_dSEXP, SEXP gamma_rSEXP, SEXP thetasSEXP, SEXP MSEXP, SEXP dG1sSEXP, SEXP dG2sSEXP, SEXP dG3sSEXP, SEXP mu_1SEXP, SEXP mu_2SEXP, SEXP mu_3SEXP, SEXP mu_rSEXP, SEXP kappaSEXP, SEXP rhoSEXP, SEXP omegaSEXP, SEXP r_betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type t_max(t_maxSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma_d(gamma_dSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma_r(gamma_rSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type thetas(thetasSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type dG1s(dG1sSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type dG2s(dG2sSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type dG3s(dG3sSEXP);
    Rcpp::traits::input_parameter< double >::type mu_1(mu_1SEXP);
    Rcpp::traits::input_parameter< double >::type mu_2(mu_2SEXP);
    Rcpp::traits::input_parameter< double >::type mu_3(mu_3SEXP);
    Rcpp::traits::input_parameter< double >::type mu_r(mu_rSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< int >::type r_beta(r_betaSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_vsveipdr_Gs_booster_reinfect(t_max, y, alpha, beta, gamma_d, gamma_r, thetas, M, dG1s, dG2s, dG3s, mu_1, mu_2, mu_3, mu_r, kappa, rho, omega, r_beta));
    return rcpp_result_gen;
END_RCPP
}
// object_kappa_tra_booster_reinfect
arma::rowvec object_kappa_tra_booster_reinfect(arma::vec x, arma::mat dat_sim, int t_V, int iter, arma::vec thetas, int M, double mu_1, double mu_2, double mu_3, double mu_r, int r_beta, arma::mat bspline2);
RcppExport SEXP _objectfunctionsbooster_object_kappa_tra_booster_reinfect(SEXP xSEXP, SEXP dat_simSEXP, SEXP t_VSEXP, SEXP iterSEXP, SEXP thetasSEXP, SEXP MSEXP, SEXP mu_1SEXP, SEXP mu_2SEXP, SEXP mu_3SEXP, SEXP mu_rSEXP, SEXP r_betaSEXP, SEXP bspline2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type dat_sim(dat_simSEXP);
    Rcpp::traits::input_parameter< int >::type t_V(t_VSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type thetas(thetasSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< double >::type mu_1(mu_1SEXP);
    Rcpp::traits::input_parameter< double >::type mu_2(mu_2SEXP);
    Rcpp::traits::input_parameter< double >::type mu_3(mu_3SEXP);
    Rcpp::traits::input_parameter< double >::type mu_r(mu_rSEXP);
    Rcpp::traits::input_parameter< int >::type r_beta(r_betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type bspline2(bspline2SEXP);
    rcpp_result_gen = Rcpp::wrap(object_kappa_tra_booster_reinfect(x, dat_sim, t_V, iter, thetas, M, mu_1, mu_2, mu_3, mu_r, r_beta, bspline2));
    return rcpp_result_gen;
END_RCPP
}
// object_kappa_booster_reinfect
arma::rowvec object_kappa_booster_reinfect(arma::vec x, arma::mat dat_sim, int t_V, int iter, arma::vec thetas, int M, double mu_1, double mu_2, double mu_3, double mu_r, int r_beta, arma::mat bspline2, int inter_b_min, int inter_b_max, int inter_l);
RcppExport SEXP _objectfunctionsbooster_object_kappa_booster_reinfect(SEXP xSEXP, SEXP dat_simSEXP, SEXP t_VSEXP, SEXP iterSEXP, SEXP thetasSEXP, SEXP MSEXP, SEXP mu_1SEXP, SEXP mu_2SEXP, SEXP mu_3SEXP, SEXP mu_rSEXP, SEXP r_betaSEXP, SEXP bspline2SEXP, SEXP inter_b_minSEXP, SEXP inter_b_maxSEXP, SEXP inter_lSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type dat_sim(dat_simSEXP);
    Rcpp::traits::input_parameter< int >::type t_V(t_VSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type thetas(thetasSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< double >::type mu_1(mu_1SEXP);
    Rcpp::traits::input_parameter< double >::type mu_2(mu_2SEXP);
    Rcpp::traits::input_parameter< double >::type mu_3(mu_3SEXP);
    Rcpp::traits::input_parameter< double >::type mu_r(mu_rSEXP);
    Rcpp::traits::input_parameter< int >::type r_beta(r_betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type bspline2(bspline2SEXP);
    Rcpp::traits::input_parameter< int >::type inter_b_min(inter_b_minSEXP);
    Rcpp::traits::input_parameter< int >::type inter_b_max(inter_b_maxSEXP);
    Rcpp::traits::input_parameter< int >::type inter_l(inter_lSEXP);
    rcpp_result_gen = Rcpp::wrap(object_kappa_booster_reinfect(x, dat_sim, t_V, iter, thetas, M, mu_1, mu_2, mu_3, mu_r, r_beta, bspline2, inter_b_min, inter_b_max, inter_l));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_objectfunctionsbooster_sim_vsveipdr_Gs_booster", (DL_FUNC) &_objectfunctionsbooster_sim_vsveipdr_Gs_booster, 18},
    {"_objectfunctionsbooster_object_kappa_tra_booster", (DL_FUNC) &_objectfunctionsbooster_object_kappa_tra_booster, 11},
    {"_objectfunctionsbooster_object_kappa_booster", (DL_FUNC) &_objectfunctionsbooster_object_kappa_booster, 14},
    {"_objectfunctionsbooster_sim_vsveipdr_Gs_booster_reinfect", (DL_FUNC) &_objectfunctionsbooster_sim_vsveipdr_Gs_booster_reinfect, 19},
    {"_objectfunctionsbooster_object_kappa_tra_booster_reinfect", (DL_FUNC) &_objectfunctionsbooster_object_kappa_tra_booster_reinfect, 12},
    {"_objectfunctionsbooster_object_kappa_booster_reinfect", (DL_FUNC) &_objectfunctionsbooster_object_kappa_booster_reinfect, 15},
    {NULL, NULL, 0}
};

RcppExport void R_init_objectfunctionsbooster(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
