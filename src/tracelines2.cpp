#include <Rcpp.h>
using namespace Rcpp;
using std::exp;

//[[Rcpp::export]]
List bernoulli_traceline2(
    NumericVector& p_item,
    NumericVector& theta,
    NumericMatrix& predictors,
    double samp_size,
    double num_quadpts,
    double num_items
){

  // There are two loops: 1) responses; 2) quadrature points

  List traceline(2);
  double p_c0 = p_item[0];
  double p_a0 = p_item[1];
  IntegerVector idx = IntegerVector::create(2,3,4);
  IntegerVector idx2 = IntegerVector::create(5,6,7);
  NumericVector p_c1 = p_item[idx];
  NumericVector p_a1 = p_item[idx2];

  // 1. Loop through number of responses
  for(int r = 0; r < 2; ++r){

    // create matrix to hold tracelines
    NumericMatrix temp_mat(samp_size,num_quadpts);


    // 2. Loop through quadrature points
    for(int q = 0; q < num_quadpts; ++q){

      for(int i = 0; i < 3; ++i){
        if(r == 0){ // response one
          temp_mat(_,q) = 1 - 1/(1 + exp(-(p_c0 + predictors(_,i) * p_c1[i] + (p_a0 + predictors(_,i) * p_a1[i]) * theta[q])));
        } else if(r == 1){ // last response
          temp_mat(_,q) = temp_mat(_,q) + 1/(1 + exp(-(p_c0 + predictors(_,i) * p_c1[i] + (p_a0 + predictors(_,i) * p_a1[i]) * theta[q])));
        }
      }
    }

    traceline[r] = temp_mat;

  }

  return traceline;
}
