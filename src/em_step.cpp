#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using std::exp;
using std::sqrt;
using std::pow;


////////////////////////////////////////////
/////////////// TRESHOLDING ////////////////
////////////////////////////////////////////

template <typename T> double sgn(T val) {
  return (T(0) < val) - (val < T(0));
}
template <class T> const T& max (const T& a, const T& b) {
  return (a<b)?b:a;     // or: return comp(a,b)?b:a; for version (2)
}

//[[Rcpp::export]]
double soft_thresh_est(
    double z,
    double alpha,
    double tau
){

  double t = (fabs(z) - tau*alpha)/(1+tau*(1-alpha));
  double p_new = sgn(z)*max(t,0.0);
  return p_new;

}

//[[Rcpp::export]]
double firm_thresh_est(
    double z,
    double alpha,
    double tau,
    double gamma
){

  double p_new;
  if(fabs(z)/(1+tau*(1-alpha)) <= gamma*tau){
    p_new = (gamma/(gamma-1))*soft_thresh_est(z,alpha,tau);
  }else{
    p_new = z/(1+tau*(1-alpha));
  }

  return(p_new);
}

//[[Rcpp::export]]
arma::vec dnormCpp(
    arma::vec theta,
    double alpha_impact,
    double phi_impact,
    int num_quadpts
){
  arma::vec posterior(num_quadpts);
  static const double pi = 3.14159265;

  for(int q = 0; q < num_quadpts; q++){
    posterior[q] = 1 / sqrt(2*pi*phi_impact)*exp(-pow((theta[q] -
      alpha_impact),2) / (2*phi_impact));
  }

  return posterior;
}

///////////////////////////////////////////
/////////////// TRACELINES ////////////////
///////////////////////////////////////////

//[[Rcpp::export]]
arma::mat bernoulli_traceline_est(
    arma::vec p_item,
    arma::mat theta,
    arma::mat predictors,
    int samp_size,
    int num_quadpts
){

  // split item parameter vector
  double p_c0 = p_item(0,0);
  double p_a0 = p_item(1,0);
  arma::vec p_c1 = p_item.subvec(2,1+predictors.n_cols);
  arma::vec p_a1 = p_item.subvec(2+predictors.n_cols,2*predictors.n_cols+1);

  // create matrix to hold traceline
  arma::mat traceline = arma::zeros(samp_size,num_quadpts);

  // Loop through quadrature points
  for(int q = 0; q < num_quadpts; ++q){
    traceline.col(q) = 1 / (1 + exp(-(p_c0 + predictors * p_c1 +
      (p_a0 + predictors * p_a1) % theta.col(q))));
  }

  return traceline;
}


//[[Rcpp::export]]
List bernoulli_traceline_cpp(
    arma::vec p_item,
    arma::mat theta,
    arma::mat predictors,
    int samp_size,
    int num_quadpts
){

  // split item parameter vector
  double p_c0 = p_item(0,0);
  double p_a0 = p_item(1,0);
  arma::vec p_c1 = p_item.subvec(2,1+predictors.n_cols);
  arma::vec p_a1 = p_item.subvec(2+predictors.n_cols,2*predictors.n_cols+1);

  // create matrix to hold traceline
  arma::mat traceline0 = arma::zeros(samp_size,num_quadpts);
  arma::mat traceline1 = arma::zeros(samp_size,num_quadpts);

  // Loop through quadrature points
  for(int q = 0; q < num_quadpts; ++q){
    traceline1.col(q) = 1/(1 + exp(-(p_c0 + predictors * p_c1 +
      (p_a0 + predictors * p_a1) % theta.col(q))));
  }
  traceline0 = 1 - traceline1;

  List traceline_list(2);
  traceline_list[0] = traceline0;
  traceline_list[1] = traceline1;
  return traceline_list;
}


////////////////////////////////////////////
/////////////// DERIVATIVES ////////////////
////////////////////////////////////////////


  /////////////// ALPHA - IMPACT ////////////////

//[[Rcpp::export]]
List d_alpha_est(
    arma::vec p_alpha,
    arma::vec p_phi,
    arma::mat etable_all,
    arma::mat theta,
    arma::mat mean_predictors,
    arma::mat var_predictors,
    int cov,
    int samp_size,
    int num_items,
    int num_quadpts
){

  // obtain latent mean and variance vectors
  arma::vec alpha_impact = mean_predictors * p_alpha;
  arma::vec phi_impact = exp(var_predictors * p_phi);

  arma::mat eta_d = arma::zeros(samp_size,num_quadpts);
  for(int q = 0; q < num_quadpts; q++){
    eta_d.col(q) = mean_predictors.col(cov);
  }

  arma::mat d1_trace = arma::zeros(samp_size,num_quadpts);
  arma::mat d2_trace = arma::zeros(samp_size,num_quadpts);
  for(int i = 0; i < samp_size; i++){
    d1_trace.row(i) = eta_d.row(i)/phi_impact[i] % (theta.row(i) -
      alpha_impact[i]);
    d2_trace.row(i) = -pow(eta_d.row(i),2)/phi_impact[i];
  }

  List dlist(2);
  dlist[0] = sum(sum(etable_all % d1_trace));
  dlist[1] = sum(sum(etable_all % d2_trace));

  return(dlist);

}


/////////////// PHI - IMPACT ////////////////

//[[Rcpp::export]]
List d_phi_est(
    arma::vec p_alpha,
    arma::vec p_phi,
    arma::mat etable_all,
    arma::mat theta,
    arma::mat mean_predictors,
    arma::mat var_predictors,
    int cov,
    int samp_size,
    int num_items,
    int num_quadpts
){

  // obtain latent mean and variance vectors
  arma::vec alpha_impact = mean_predictors * p_alpha;
  arma::vec phi_impact = exp(var_predictors * p_phi);

  arma::vec eta_d1 = sqrt(phi_impact) / 2 % var_predictors.col(cov);
  arma::vec eta_d2 = sqrt(phi_impact) / 4 % pow(var_predictors.col(cov),2);

  arma::mat d1_trace = arma::zeros(samp_size,num_quadpts);
  arma::mat d2_trace = arma::zeros(samp_size,num_quadpts);
  for(int i = 0; i < samp_size; i++){
    d1_trace.row(i) = eta_d1[i] * (pow((theta.row(i) - alpha_impact[i]),2) /
      pow(phi_impact[i],(3.0/2)) - 1 / sqrt(phi_impact[i]));
    d2_trace.row(i) = -2 * eta_d2[i] * (pow(phi_impact[i],(-3.0/2)) *
      pow((theta.row(i) - alpha_impact[i]),2));
  }

  List dlist(2);
  dlist[0] = sum(sum(etable_all % d1_trace));
  dlist[1] = sum(sum(etable_all % d2_trace));

  return(dlist);

}


/////////////// BERNOULLI - DIF ////////////////

//[[Rcpp::export]]
List d_bernoulli_est(
    std::string parm,
    arma::vec p_item,
    arma::mat etable1,
    arma::mat etable2,
    arma::mat theta,
    arma::mat predictors,
    int cov,
    int samp_size,
    int num_items,
    int num_quadpts
){

  arma::mat eta_d = arma::zeros(samp_size,num_quadpts);
  if(parm == "c0"){
    eta_d.ones();
  } else if(parm == "a0"){
    for(int q = 0; q < num_quadpts; q++){
      eta_d.col(q) = arma::ones(samp_size) % theta.col(q);
    }
  } else if(parm == "c1"){
    for(int q = 0; q < num_quadpts; q++){
      eta_d.col(q) = predictors.col(cov);
    }
  } else if(parm == "a1"){
    for(int q = 0; q < num_quadpts; q++){
      eta_d.col(q) = predictors.col(cov) % theta.col(q);
    }
  }

  arma::mat traceline = bernoulli_traceline_est(p_item,
                                                theta,
                                                predictors,
                                                samp_size,
                                                num_quadpts);

  List dlist(2);
  dlist[0] = sum(sum((1-traceline) % eta_d % etable2)) +
    sum(sum(-traceline % eta_d % etable1));
  dlist[1] = sum(sum(-traceline % (1-traceline) % pow(eta_d,2) % etable1)) +
    sum(sum(-traceline % (1-traceline) % pow(eta_d,2) % etable2));

  return(dlist);

}
