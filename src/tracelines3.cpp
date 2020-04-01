#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using std::exp;

//[[Rcpp::export]]
List bernoulli_traceline3(
    arma::vec p_item,
    arma::vec theta,
    arma::mat predictors,
    double samp_size,
    double num_quadpts
){

  // There are two loops: 1) responses; 2) quadrature points
  double p_c0 = p_item(0,0);
  double p_a0 = p_item(1,0);
  arma::vec p_c1 = p_item.subvec(2,1+predictors.n_cols);
  arma::vec p_a1 = p_item.subvec(2+predictors.n_cols,2*predictors.n_cols+1);

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
