// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::export]]
arma::mat sim_vsveipdr_Gs_booster(int t_max, arma::vec y, arma::vec alpha, arma::vec beta, arma::vec gamma_d, 
                                  arma::vec gamma_r, arma::vec thetas, int M, arma::vec dG1s, arma::vec dG2s, arma::vec dG3s,
                                  double mu_1, double mu_2, double mu_3, double kappa, double rho, double omega, int r_beta) {
  arma::mat dat_gen(t_max+1, 16, arma::fill::zeros);
  //data.frame(t = t, S = S, V1 = V1, V2 = V2 , Ea = Ea, Ep = Ep, I = I, R.a = R.a, 
  //           R.d = R.d, R.r = R.r, G1 = G1, G2 = G2, R = R.d + R.r, N = I + R.d + R.r, V3 = V3, G3 = G3)
  dat_gen(0, 0) = 1;
  dat_gen(0, 2) = y(0);
  dat_gen(0, 3) = y(1);
  dat_gen(0, 4) = y(2);
  dat_gen(0, 5) = y(3);
  dat_gen(0, 6) = y(4);
  dat_gen(0, 8) = y(5);
  dat_gen(0, 9) = y(6);
  dat_gen(0, 7) = y(7);
  dat_gen(0, 14) = y(8);
  dat_gen(0, 1) = M - arma::sum(y);
  
  double dG1, dG2, dG3, dV1, dV2, dV3, dI_S, dI_V1, dI_V2, dI_V3, dEa, dEp, dR_r, dR_d, theta;
  
  for (int i = 0; i < t_max; i = i + 1) {
    theta = thetas(i);
    dG1 = dG1s(i);
    dG2 = dG2s(i);
    dG3 = dG3s(i);
    dV1 = abs(mu_1 * dat_gen(i, 2));
    dV2 = abs(mu_2 * dat_gen(i, 3));
    dV3 = abs(mu_3 * dat_gen(i, 14));
    dI_S = abs(dat_gen(i, 1) / M * (beta(i) * dat_gen(i, 5) + 
      beta(i) / r_beta * dat_gen(i, 6) + beta(i) / r_beta * dat_gen(i, 4)) );
    dI_V1 = abs(kappa * rho * dat_gen(i, 2) / M * (beta(i) * dat_gen(i, 5) + 
      beta(i) / r_beta * dat_gen(i, 6) + beta(i) / r_beta * dat_gen(i, 4)));
    dI_V2 = abs(kappa * dat_gen(i, 3) / M * (beta(i) * dat_gen(i, 5) + 
      beta(i) / r_beta * dat_gen(i, 6) + beta(i) / r_beta * dat_gen(i, 4)));
    dI_V3 = abs(kappa * omega * dat_gen(i, 14) / M * (beta(i) * dat_gen(i, 5) + 
      beta(i) / r_beta * dat_gen(i, 6) + beta(i) / r_beta * dat_gen(i, 4)));
    dEa = abs(gamma_r(i) * dat_gen(i, 4));
    dEp = abs(alpha(i) * dat_gen(i, 5));
    dR_r = abs(gamma_r(i) * dat_gen(i, 6));
    dR_d = abs(gamma_d(i) * dat_gen(i, 6));
    
    dV1 = std::max(std::max(std::min(R::rpois(dV1), dat_gen(i, 2) + dG1 - dG2), 0.0), dG1 - dat_gen(i, 1));
    dV2 = std::max(std::max(std::min(R::rpois(dV2), dat_gen(i, 3) + dG2 - dG3), 0.0), dG1 - dat_gen(i, 1) - dV1);
    dV3 = std::max(std::max(std::min(R::rpois(dV3), dat_gen(i, 14) + dG3), 0.0), dG1 - dat_gen(i, 1) - dV1 - dV2);
    dI_S = std::max(std::min(R::rpois(dI_S), dat_gen(i, 1) - dG1 + dV1 + dV2 + dV3), 0.0);
    dI_V1 = std::max(std::min(R::rpois(dI_V1), dat_gen(i, 2) + dG1 - dV1 - dG2), 0.0);
    dI_V2 = std::max(std::min(R::rpois(dI_V2), dat_gen(i, 3) + dG2 - dV2 - dG3), 0.0);
    dI_V3 = std::max(std::min(R::rpois(dI_V3), dat_gen(i, 14) + dG3 - dV3), 0.0);
    dEa = std::max(std::min(R::rpois(dEa), dat_gen(i, 4) + (dI_S + dI_V1 + dI_V2 + dI_V3) * (1 - theta)), 0.0);
    dEp = std::max(std::min(R::rpois(dEp), dat_gen(i, 5) + (dI_S + dI_V1 + dI_V2 + dI_V3) * theta), 0.0);
    dR_d = std::max(std::min(R::rpois(dR_d), dat_gen(i, 6) + dEp), 0.0);
    dR_r = std::max(std::min(R::rpois(dR_r), dat_gen(i, 6) + dEp - dR_d), 0.0);
    
    dat_gen(i+1, 0) = i+2;
    dat_gen(i+1, 10) = std::max(dat_gen(i, 10) + dG1, 0.0);
    dat_gen(i+1, 11) = std::max(dat_gen(i, 11) + dG2, 0.0);
    dat_gen(i+1, 15) = std::max(dat_gen(i, 15) + dG3, 0.0);
    dat_gen(i+1, 1) = std::max(dat_gen(i, 1) - dI_S - dG1 + dV1 + dV2 + dV3, 0.0);
    dat_gen(i+1, 2) = std::max(dat_gen(i, 2) + dG1 - dV1 - dI_V1 - dG2, 0.0);
    dat_gen(i+1, 3) = std::max(dat_gen(i, 3) + dG2 - dV2 - dI_V2 - dG3, 0.0);
    dat_gen(i+1, 14) = std::max(dat_gen(i, 14) + dG3 - dV3 - dI_V3, 0.0);
    dat_gen(i+1, 4) = std::max(dat_gen(i, 4) + (dI_S + dI_V1 + dI_V2 + dI_V3) * (1 - theta) - dEa, 0.0);
    dat_gen(i+1, 5) = std::max(dat_gen(i, 5) + (dI_S + dI_V1 + dI_V2 + dI_V3) * theta - dEp, 0.0);
    dat_gen(i+1, 6) = std::max(dat_gen(i, 6) + dEp - dR_d - dR_r, 0.0) ;
    dat_gen(i+1, 8) = std::max(dat_gen(i, 8) + dR_d, 0.0);
    dat_gen(i+1, 9) = std::max(dat_gen(i, 9) + dR_r, 0.0);
    dat_gen(i+1, 7) = std::max(dat_gen(i, 7) + dEa, 0.0);
  }
  
  dat_gen.col(12) = dat_gen.col(8) + dat_gen.col(9);
  dat_gen.col(13) = dat_gen.col(6) + dat_gen.col(12);
  return dat_gen;
}

// [[Rcpp::export]]
arma::rowvec object_kappa_tra_booster(arma::vec x, arma::mat dat_sim, int t_V, int iter, 
                                      arma::vec thetas, int M, double mu_1, double mu_2, double mu_3, int r_beta, arma::mat bspline2){
  
  arma::vec lb(t_V+1, arma::fill::value(0.001));
  arma::vec ub(t_V+1, arma::fill::value(1.0));
  arma::vec beta_try = arma::min(arma::max(bspline2 * x.subvec(3, x.size()-1), lb), ub);
  arma::mat traject(t_V+1, 16);
  arma::mat traject_s(t_V+1, 3, arma::fill::zeros);
  arma::uvec ind_init1 = {0};
  arma::uvec ind_init2 = {14, 15, 11, 10, 5, 2, 3, 9, 17};
  arma::vec y = dat_sim(ind_init1, ind_init2).t();
  arma::uvec ind_traject = {5, 6, 13};
  arma::vec alpha_try(t_V+1, arma::fill::value(dat_sim(0, 8)));
  
  for (int jj = 0; jj < iter; jj = jj + 1){
    traject = sim_vsveipdr_Gs_booster(t_V, y, alpha_try, beta_try, dat_sim.col(12), 
                                      dat_sim.col(13), thetas, M, diff(dat_sim.col(6)), diff(dat_sim.col(7)), diff(dat_sim.col(16)),
                                      mu_1, mu_2, mu_3, x(1), x(0), x(2), r_beta);
    traject_s = traject_s + traject.cols(ind_traject);
  }
  
  arma::uvec ind_crit = {10, 5, 1};
  arma::mat truth = dat_sim.cols(ind_crit);
  arma::mat object_k = traject_s/iter/truth-1;
  object_k.elem( arma::find(truth == 0) ).zeros();
  
  return arma::pow(object_k.rows(1, t_V).as_col(), 2).t(); 
}

// [[Rcpp::export]]
arma::rowvec object_kappa_booster(arma::vec x, arma::mat dat_sim, int t_V, int iter, 
                                  arma::vec thetas, int M, double mu_1, double mu_2, double mu_3, int r_beta, arma::mat bspline2,
                                  int inter_b_min, int inter_b_max, int inter_l){
  arma::vec lb(t_V+1, arma::fill::value(0.001));
  arma::vec ub(t_V+1, arma::fill::value(1.0));
  arma::vec beta_try = arma::min(arma::max(bspline2 * x.subvec(3, x.size()-1), lb), ub);
  arma::mat traject(t_V+1, 16);
  arma::mat traject_s(t_V+1, 3, arma::fill::zeros);
  arma::uvec ind_init1 = {0};
  arma::uvec ind_init2 = {14, 15, 11, 10, 5, 2, 3, 9, 17};
  arma::vec y = dat_sim(ind_init1, ind_init2).t();
  arma::uvec ind_traject = {5, 6, 13};
  arma::vec alpha_try(t_V+1, arma::fill::value(dat_sim(0, 8)));
  
  for (int jj = 0; jj < iter; jj = jj + 1){
    traject = sim_vsveipdr_Gs_booster(t_V, y, alpha_try, beta_try, dat_sim.col(12), 
                                      dat_sim.col(13), thetas, M, diff(dat_sim.col(6)), diff(dat_sim.col(7)), diff(dat_sim.col(16)),
                                      mu_1, mu_2, mu_3, x(1), x(0), x(2), r_beta);
    traject_s = traject_s + traject.cols(ind_traject);
  }
  
  arma::uvec ind_crit = {10, 5, 1};
  arma::mat truth = dat_sim.cols(ind_crit);
  arma::mat object_k = traject_s/iter/truth-1;
  object_k.elem( arma::find(truth == 0) ).zeros();
  
  arma::ivec ind_inter_bmin = Rcpp::seq_len( inter_l+1 ) + inter_b_min -1;
  arma::ivec exten1(inter_b_max - inter_b_min + 1, arma::fill::ones);
  arma::ivec lag_b = Rcpp::seq_len(inter_b_max - inter_b_min + 1) -1;
  arma::ivec exten2(inter_l + 1, arma::fill::ones);
  arma::imat ind_criter = ind_inter_bmin * exten1.t() + exten2 * lag_b.t();
  arma::uvec ind_Ep_po = arma::conv_to< arma::uvec >::from(ind_criter.as_col());
  arma::uvec ind_Ep = {0};
  arma::mat ind_all_cri = arma::pow( object_k(ind_Ep_po, ind_Ep), 2 );
  ind_all_cri.reshape( inter_l + 1, inter_b_max - inter_b_min + 1);
  
  return arma::sqrt(arma::mean(ind_all_cri, 0));
}

// [[Rcpp::export]]
arma::mat sim_vsveipdr_Gs_booster_reinfect(int t_max, arma::vec y, arma::vec alpha, arma::vec beta, arma::vec gamma_d, 
                                           arma::vec gamma_r, arma::vec thetas, int M, arma::vec dG1s, arma::vec dG2s, arma::vec dG3s,
                                           double mu_1, double mu_2, double mu_3, double mu_r, double kappa, double rho, double omega, int r_beta) {
  arma::mat dat_gen(t_max+1, 16, arma::fill::zeros);
  //data.frame(t = t, S = S, V1 = V1, V2 = V2 , Ea = Ea, Ep = Ep, I = I, R.a = R.a, 
  //           R.d = R.d, R.r = R.r, G1 = G1, G2 = G2, R = R.d + R.r, N = N, V3 = V3, G3 = G3)
  dat_gen(0, 0) = 1;
  dat_gen(0, 2) = y(0);
  dat_gen(0, 3) = y(1);
  dat_gen(0, 4) = y(2);
  dat_gen(0, 5) = y(3);
  dat_gen(0, 6) = y(4);
  dat_gen(0, 8) = y(5);
  dat_gen(0, 9) = y(6);
  dat_gen(0, 7) = y(7);
  dat_gen(0, 14) = y(8);
  dat_gen(0, 1) = M - arma::sum(y.subvec(0, 8));
  dat_gen(0, 13) = y(9);
  
  double dG1, dG2, dG3, dV1, dV2, dV3, dI_S, dI_V1, dI_V2, dI_V3, dEa, dEp, dR_r, dR_d, theta, RrS, RaS;
  
  for (int i = 0; i < t_max; i = i + 1) {
    theta = thetas(i);
    dG1 = dG1s(i);
    dG2 = dG2s(i);
    dG3 = dG3s(i);
    dV1 = abs(mu_1 * dat_gen(i, 2));
    dV2 = abs(mu_2 * dat_gen(i, 3));
    dV3 = abs(mu_3 * dat_gen(i, 14));
    dI_S = abs(dat_gen(i, 1) / M * (beta(i) * dat_gen(i, 5) + 
      beta(i) / r_beta * dat_gen(i, 6) + beta(i) / r_beta * dat_gen(i, 4)) );
    dI_V1 = abs(kappa * rho * dat_gen(i, 2) / M * (beta(i) * dat_gen(i, 5) + 
      beta(i) / r_beta * dat_gen(i, 6) + beta(i) / r_beta * dat_gen(i, 4)));
    dI_V2 = abs(kappa * dat_gen(i, 3) / M * (beta(i) * dat_gen(i, 5) + 
      beta(i) / r_beta * dat_gen(i, 6) + beta(i) / r_beta * dat_gen(i, 4)));
    dI_V3 = abs(kappa * omega * dat_gen(i, 14) / M * (beta(i) * dat_gen(i, 5) + 
      beta(i) / r_beta * dat_gen(i, 6) + beta(i) / r_beta * dat_gen(i, 4)));
    dEa = abs(gamma_r(i) * dat_gen(i, 4));
    dEp = abs(alpha(i) * dat_gen(i, 5));
    dR_r = abs(gamma_r(i) * dat_gen(i, 6));
    dR_d = abs(gamma_d(i) * dat_gen(i, 6));
    RrS = abs(mu_r * dat_gen(i, 9));
    RaS = abs(mu_r * dat_gen(i, 7));
    
    RrS = std::max(std::max(std::min(R::rpois(RrS), dat_gen(i, 9)), 0.0), dG1 - dat_gen(i, 1));
    RaS = std::max(std::max(std::min(R::rpois(RaS), dat_gen(i, 7)), 0.0), dG1 - dat_gen(i, 1) - RrS);
    dV1 = std::max(std::max(std::min(R::rpois(dV1), dat_gen(i, 2) + dG1 - dG2), 0.0), dG1 - dat_gen(i, 1) - RrS - RaS);
    dV2 = std::max(std::max(std::min(R::rpois(dV2), dat_gen(i, 3) + dG2 - dG3), 0.0), dG1 - dat_gen(i, 1) - RrS - RaS - dV1);
    dV3 = std::max(std::max(std::min(R::rpois(dV3), dat_gen(i, 14) + dG3), 0.0), dG1 - dat_gen(i, 1) - RrS - RaS - dV1 - dV2);
    dI_S = std::max(std::min(R::rpois(dI_S), dat_gen(i, 1) - dG1 + dV1 + dV2 + dV3 + RrS + RaS), 0.0);
    dI_V1 = std::max(std::min(R::rpois(dI_V1), dat_gen(i, 2) + dG1 - dV1 - dG2), 0.0);
    dI_V2 = std::max(std::min(R::rpois(dI_V2), dat_gen(i, 3) + dG2 - dV2 - dG3), 0.0);
    dI_V3 = std::max(std::min(R::rpois(dI_V3), dat_gen(i, 14) + dG3 - dV3), 0.0);
    dEa = std::max(std::min(R::rpois(dEa), dat_gen(i, 4) + (dI_S + dI_V1 + dI_V2 + dI_V3) * (1 - theta)), 0.0);
    dEp = std::max(std::min(R::rpois(dEp), dat_gen(i, 5) + (dI_S + dI_V1 + dI_V2 + dI_V3) * theta), 0.0);
    dR_d = std::max(std::min(R::rpois(dR_d), dat_gen(i, 6) + dEp), 0.0);
    dR_r = std::max(std::min(R::rpois(dR_r), dat_gen(i, 6) + dEp - dR_d), 0.0);
    
    dat_gen(i+1, 0) = i+2;
    dat_gen(i+1, 10) = std::max(dat_gen(i, 10) + dG1, 0.0);
    dat_gen(i+1, 11) = std::max(dat_gen(i, 11) + dG2, 0.0);
    dat_gen(i+1, 15) = std::max(dat_gen(i, 15) + dG3, 0.0);
    dat_gen(i+1, 1) = std::max(dat_gen(i, 1) - dI_S - dG1 + dV1 + dV2 + dV3 + RrS + RaS, 0.0);
    dat_gen(i+1, 2) = std::max(dat_gen(i, 2) + dG1 - dV1 - dI_V1 - dG2, 0.0);
    dat_gen(i+1, 3) = std::max(dat_gen(i, 3) + dG2 - dV2 - dI_V2 - dG3, 0.0);
    dat_gen(i+1, 14) = std::max(dat_gen(i, 14) + dG3 - dV3 - dI_V3, 0.0);
    dat_gen(i+1, 4) = std::max(dat_gen(i, 4) + (dI_S + dI_V1 + dI_V2 + dI_V3) * (1 - theta) - dEa, 0.0);
    dat_gen(i+1, 5) = std::max(dat_gen(i, 5) + (dI_S + dI_V1 + dI_V2 + dI_V3) * theta - dEp, 0.0);
    dat_gen(i+1, 6) = std::max(dat_gen(i, 6) + dEp - dR_d - dR_r, 0.0) ;
    dat_gen(i+1, 8) = std::max(dat_gen(i, 8) + dR_d, 0.0);
    dat_gen(i+1, 9) = std::max(dat_gen(i, 9) + dR_r - RrS, 0.0);
    dat_gen(i+1, 7) = std::max(dat_gen(i, 7) + dEa - RaS, 0.0);
    dat_gen(i+1, 13) = std::max(dat_gen(i, 13) + dEp, 0.0);
  }
  
  dat_gen.col(12) = dat_gen.col(8) + dat_gen.col(9);
  
  return dat_gen;
}

// [[Rcpp::export]]
arma::rowvec object_kappa_tra_booster_reinfect(arma::vec x, arma::mat dat_sim, int t_V, int iter, 
                                               arma::vec thetas, int M, double mu_1, double mu_2, double mu_3, double mu_r, int r_beta, arma::mat bspline2){
  
  arma::vec lb(t_V+1, arma::fill::value(0.001));
  arma::vec ub(t_V+1, arma::fill::value(1.0));
  arma::vec beta_try = arma::min(arma::max(bspline2 * x.subvec(3, x.size()-1), lb), ub);
  arma::mat traject(t_V+1, 16);
  arma::mat traject_s(t_V+1, 3, arma::fill::zeros);
  arma::uvec ind_init1 = {0};
  arma::uvec ind_init2 = {14, 15, 11, 10, 5, 2, 3, 9, 17, 1};
  arma::vec y = dat_sim(ind_init1, ind_init2).t();
  arma::uvec ind_traject = {5, 6, 13};
  arma::vec alpha_try(t_V+1, arma::fill::value(dat_sim(0, 8)));
  
  for (int jj = 0; jj < iter; jj = jj + 1){
    traject = sim_vsveipdr_Gs_booster_reinfect(t_V, y, alpha_try, beta_try, dat_sim.col(12), 
                                               dat_sim.col(13), thetas, M, diff(dat_sim.col(6)), diff(dat_sim.col(7)), diff(dat_sim.col(16)),
                                               mu_1, mu_2, mu_3, mu_r, x(1), x(0), x(2), r_beta);
    traject_s = traject_s + traject.cols(ind_traject);
  }
  
  arma::uvec ind_crit = {10, 5, 1};
  arma::mat truth = dat_sim.cols(ind_crit);
  arma::mat object_k = traject_s/iter/truth-1;
  object_k.elem( arma::find(truth == 0) ).zeros();
  
  return arma::pow(object_k.rows(1, t_V).as_col(), 2).t(); 
}

// [[Rcpp::export]]
arma::rowvec object_kappa_booster_reinfect(arma::vec x, arma::mat dat_sim, int t_V, int iter, 
                                           arma::vec thetas, int M, double mu_1, double mu_2, double mu_3, double mu_r, int r_beta, arma::mat bspline2,
                                           int inter_b_min, int inter_b_max, int inter_l){
  arma::vec lb(t_V+1, arma::fill::value(0.001));
  arma::vec ub(t_V+1, arma::fill::value(1.0));
  arma::vec beta_try = arma::min(arma::max(bspline2 * x.subvec(3, x.size()-1), lb), ub);
  arma::mat traject(t_V+1, 16);
  arma::mat traject_s(t_V+1, 3, arma::fill::zeros);
  arma::uvec ind_init1 = {0};
  arma::uvec ind_init2 = {14, 15, 11, 10, 5, 2, 3, 9, 17, 1};
  arma::vec y = dat_sim(ind_init1, ind_init2).t();
  arma::uvec ind_traject = {5, 6, 13};
  arma::vec alpha_try(t_V+1, arma::fill::value(dat_sim(0, 8)));
  
  for (int jj = 0; jj < iter; jj = jj + 1){
    traject = sim_vsveipdr_Gs_booster_reinfect(t_V, y, alpha_try, beta_try, dat_sim.col(12), 
                                               dat_sim.col(13), thetas, M, diff(dat_sim.col(6)), diff(dat_sim.col(7)), diff(dat_sim.col(16)),
                                               mu_1, mu_2, mu_3, mu_r, x(1), x(0), x(2), r_beta);
    traject_s = traject_s + traject.cols(ind_traject);
  }
  
  arma::uvec ind_crit = {10, 5, 1};
  arma::mat truth = dat_sim.cols(ind_crit);
  arma::mat object_k = traject_s/iter/truth-1;
  object_k.elem( arma::find(truth == 0) ).zeros();
  
  arma::ivec ind_inter_bmin = Rcpp::seq_len( inter_l+1 ) + inter_b_min -1;
  arma::ivec exten1(inter_b_max - inter_b_min + 1, arma::fill::ones);
  arma::ivec lag_b = Rcpp::seq_len(inter_b_max - inter_b_min + 1) -1;
  arma::ivec exten2(inter_l + 1, arma::fill::ones);
  arma::imat ind_criter = ind_inter_bmin * exten1.t() + exten2 * lag_b.t();
  arma::uvec ind_Ep_po = arma::conv_to< arma::uvec >::from(ind_criter.as_col());
  arma::uvec ind_Ep = {0};
  arma::mat ind_all_cri = arma::pow( object_k(ind_Ep_po, ind_Ep), 2 );
  ind_all_cri.reshape( inter_l + 1, inter_b_max - inter_b_min + 1);
  
  return arma::sqrt(arma::mean(ind_all_cri, 0));
}
