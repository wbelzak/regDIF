#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using std::exp;

//[[Rcpp::export]]
List bernoulli_traceline(
    double p_c0,
    double p_a0,
    arma::vec p_c1,
    arma::vec p_a1,
    arma::vec theta,
    arma::mat predictors,
    double samp_size,
    double num_quadpts
){

  // There are two loops: 1) responses; 2) quadrature points


  List traceline(2);
  // 1. Loop through number of responses
  for(int r = 0; r < 2; ++r){

    // create matrix to hold tracelines
    arma::mat temp_mat = arma::zeros(samp_size,num_quadpts);

    // 2. Loop through quadrature points
    for(int q = 0; q < num_quadpts; ++q){


      if(r == 0){ // response one
        temp_mat.col(q) = 1 - 1/(1 + exp(-(p_c0 + predictors * p_c1 + (p_a0 + predictors * p_a1) * theta[q])));
      } else if(r == 1){ // last response
        temp_mat.col(q) = 1/(1 + exp(-(p_c0 + predictors * p_c1 + (p_a0 + predictors * p_a1) * theta[q])));
      }

    }

    traceline[r] = temp_mat;

  }

  return traceline;
}

//[[Rcpp::export]]
List categorical_traceline(
    double p_c0,
    double p_a0,
    arma::vec p_thr,
    arma::vec p_c1,
    arma::vec p_a1,
    arma::vec theta,
    arma::mat predictors,
    double samp_size,
    double num_quadpts,
    int num_responses_item
){

  // There are two loops: 1) responses; 2) quadrature points


  List traceline(num_responses_item);
  // 1. Loop through number of responses
  for(int r = 0; r < num_responses_item; ++r){

    // create matrix to hold tracelines
    arma::mat temp_mat = arma::zeros(samp_size,num_quadpts);

    // 2. Loop through quadrature points
    for(int q = 0; q < num_quadpts; ++q){


      if(r == 0){ // response one
        temp_mat.col(q) = 1 - 1/(1 + exp(-(p_c0 + predictors * p_c1 + (p_a0 + predictors * p_a1) * theta[q])));
      } else if((r == 1) && (r < (num_responses_item - 1))){ // response two (if this is ordinal item)
        temp_mat.col(q) = 1/(1 + exp(-(p_c0 + predictors * p_c1 + (p_a0 + predictors * p_a1) * theta[q]))) -
                          1/(1 + exp(-(p_c0 - p_thr[0] + predictors * p_c1 + (p_a0 + predictors * p_a1) * theta[q])));
      } else if((r > 1) && (r < (num_responses_item - 1))){ // responses 3 to (num_responses_item - 1)
        temp_mat.col(q) = 1/(1 + exp(-(p_c0 - p_thr[r-2] + predictors * p_c1 + (p_a0 + predictors * p_a1) * theta[q]))) -
                          1/(1 + exp(-(p_c0 - p_thr[r-1] + predictors * p_c1 + (p_a0 + predictors * p_a1) * theta[q])));
      } else if(r == (num_responses_item - 1)){ // last response
        temp_mat.col(q) = 1/(1 + exp(-(p_c0 - p_thr[(num_responses_item - 3)] + predictors * p_c1 + (p_a0 + predictors * p_a1) * theta[q])));
      }

    }

    traceline[r] = temp_mat;

  }


  return traceline;
}

//[[Rcpp::export]]
List cumulative_traceline(
    double p_c0,
    double p_a0,
    arma::vec p_thr,
    arma::vec p_c1,
    arma::vec p_a1,
    arma::vec theta,
    arma::mat predictors,
    double samp_size,
    double num_quadpts,
    int num_responses_item
){

  // There are two loops: 1) responses; 2) quadrature points


  List traceline(num_responses_item - 1);
  // 1. Loop through number of responses
  for(int r = 0; r < (num_responses_item - 1); ++r){

    // create matrix to hold tracelines
    arma::mat temp_mat = arma::zeros(samp_size,num_quadpts);

    // 2. Loop through quadrature points
    for(int q = 0; q < num_quadpts; ++q){

      if(r == 0){
        temp_mat.col(q) = 1/(1 + exp(-(p_c0 + predictors * p_c1 + (p_a0 + predictors * p_a1) * theta[q])));
      } else if((r > 0)){
        temp_mat.col(q) = 1/(1 + exp(-(p_c0 - p_thr[r-1] + predictors * p_c1 + (p_a0 + predictors * p_a1) * theta[q])));
      }


    }

    traceline[r] = temp_mat;

  }

  return traceline;
}

/*
//[[Rcpp::export]]
arma::mat continuous_traceline(
    double p_c0,
    double p_a0,
    arma::vec p_c1,
    arma::vec p_a1,
    double p_s0,
    arma::vec p_s1,
    arma::vec theta,
    arma::mat predictors,
    double samp_size,
    double num_quadpts
){
    arma::mat traceline = arma::zeros(samp_size,num_quadpts);
    arma::mat mu = arma::zeros(samp_size,num_quadpts);
    arma::mat sigma = arma::zeros(samp_size);

    sigma = p_s0 * exp(predictors * p_s1);

    for(int q = 0; q < num_quadpts; ++q){
      mu.col(q) = p_c0 + predictors * p_c1 + (p_a0 + predictors * p_a1) * theta[q];
    }

    for(int i = 0; i < samp_size; ++i){
      traceline.col(i) = 1/sigma;
    }


  return traceline;
}
*/
