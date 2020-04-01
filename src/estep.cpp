#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using std::exp;
using std::sqrt;

//[[Rcpp::export]]
List bernoulli_traceline4(
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


//[[Rcpp::export]]
List estep(
    List p,
    arma::vec theta,
    arma::mat responses,
    arma::mat predictors,
    double samp_size,
    double num_items,
    arma::vec num_responses,
    double num_quadpts
){

  List traceline(num_items);
  List etable(num_items);
  arma::mat etable_all = arma::zeros(samp_size,num_quadpts);

  arma::vec p_g = p[num_items];
  arma::vec p_b = p[num_items+1];
  arma::vec alpha = predictors * p_g;
  arma::vec phi = exp(predictors * p_b);

  List temp_ls(num_items);
  for(int j = 0; j < num_items; j++){
    traceline[j] = bernoulli_traceline4(p[j],theta,predictors,samp_size,num_quadpts);
  }

  arma::vec posterior(num_quadpts);

  for(int i = 0; i < samp_size; i++){
    for(int q = 0; q < num_quadpts; q++){
     posterior[q] = R::dnorm(theta[q], alpha[i], sqrt(phi[i]), FALSE);
    }


    for(int j = 0; j < num_items; j++){
    // double x = responses(i,j);
     temp_ls[j] = traceline[j];
    // arma::mat temp_mat = temp_ls[x];
    // for(int q = 0; q < num_quadpts; q++){
    //   posterior[q] = posterior[q] * temp_mat(i,q);
    // }

    }
  }


  return List::create(alpha,phi,traceline,temp_ls);
}
