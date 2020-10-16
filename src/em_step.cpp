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
    double lambda
){

  double t = (fabs(z) - lambda*alpha)/(1+lambda*(1-alpha));
  double p_new = sgn(z)*max(t,0.0);
  return p_new;

}

//[[Rcpp::export]]
double firm_thresh_est(
    double z,
    double alpha,
    double lambda,
    double gamma
){

  double p_new;
  if(fabs(z)/(1+lambda*(1-alpha)) <= gamma*lambda){
    p_new = (gamma/(gamma-1))*soft_thresh_est(z,alpha,lambda);
  }else{
    p_new = z/(1+lambda*(1-alpha));
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
    posterior[q] = 1/sqrt(2*pi*phi_impact)*exp(-pow((theta[q] - alpha_impact),2)/(2*phi_impact));
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
    traceline.col(q) = 1/(1 + exp(-(p_c0 + predictors * p_c1 + (p_a0 + predictors * p_a1) % theta.col(q))));
  }

  return traceline;
}


//[[Rcpp::export]]
List bernoulli_traceline_cpp(
    arma::vec p_item,
    arma::vec theta,
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
    traceline1.col(q) = 1/(1 + exp(-(p_c0 + predictors * p_c1 + (p_a0 + predictors * p_a1) * theta[q])));
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
    d1_trace.row(i) = eta_d.row(i)/phi_impact[i] % (theta.row(i) - alpha_impact[i]);
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
    d1_trace.row(i) = eta_d1[i] * (pow((theta.row(i) - alpha_impact[i]),2) / pow(phi_impact[i],(3.0/2)) - 1/sqrt(phi_impact[i]));
    d2_trace.row(i) = -2 * eta_d2[i] * (pow(phi_impact[i],(-3.0/2)) * pow((theta.row(i) - alpha_impact[i]),2));
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

  arma::mat traceline = bernoulli_traceline_est(p_item,theta,predictors,samp_size,num_quadpts);

  List dlist(2);
  dlist[0] = sum(sum((1-traceline) % eta_d % etable2)) + sum(sum(-traceline % eta_d % etable1));
  dlist[1] = sum(sum(-traceline % (1-traceline) % pow(eta_d,2) % etable1)) + sum(sum(-traceline % (1-traceline) % pow(eta_d,2) % etable2));

  return(dlist);

}


////////////////////////////////////////
/////////////// EM-STEP ////////////////
////////////////////////////////////////


//[[Rcpp::export]]
arma::mat em_step2(
    List p,
    const arma::vec& theta,
    const arma::mat& responses,
    const arma::mat& predictors,
    const arma::mat& mean_predictors,
    const arma::mat& var_predictors,
    int samp_size,
    int num_items,
    const arma::vec& num_responses,
    int num_quadpts
){

  ///////////////////////////////////////
  /////////////// E-STEP ////////////////
  ///////////////////////////////////////

  List p_working(num_items+2);
  for(int k = 0; k < (num_items+2); k++){
    arma::vec p_vec = p[k];
    p_working[k] = p_vec;
  }

  arma::vec p_alpha = p_working[num_items];
  arma::vec p_phi = p_working[num_items+1];
  arma::vec alpha_impact = mean_predictors * p_alpha;
  arma::vec phi_impact = exp(var_predictors * p_phi);

  List traceline(num_items);
  List etable(num_items);
  for(int j = 0; j < num_items; j++){
    traceline[j] = bernoulli_traceline_est(p_working[j],theta,predictors,samp_size,num_quadpts);
    etable[j] = arma::zeros(samp_size,num_quadpts);
  }

  arma::mat etable_all = arma::zeros(samp_size,num_quadpts);
  arma::vec posterior = arma::zeros(num_quadpts);
  for(int i = 0; i < samp_size; i++){

    posterior = dnormCpp(theta, alpha_impact[i], phi_impact[i], num_quadpts);

    for(int j = 0; j < num_items; j++){
      double resp = responses(i,j);
      if(resp == 1){
        posterior = posterior % (1 - as<arma::mat>(traceline[j]).row(i).t());
      } else if(resp == 2){
        posterior = posterior % as<arma::mat>(traceline[j]).row(i).t();
      } else{
        posterior = posterior;
      }
    }

    double marginal = sum(posterior);
    posterior = posterior/marginal;
    etable_all.row(i) = posterior.t();

    for(int j = 0; j < num_items; j++){
      as<arma::mat>(etable[j]).row(i) = as<arma::mat>(etable[j]).row(i) + posterior.t();
    }

  }

  return etable_all;
}





//[[Rcpp::export]]
List em_step(
    List p,
    arma::vec theta,
    arma::mat responses,
    arma::mat predictors,
    arma::mat mean_predictors,
    arma::mat var_predictors,
    StringVector itemtypes,
    StringVector penalty,
    arma::vec lambda,
    int pen,
    double alpha,
    double gamma,
    arma::vec anchor,
    bool rasch,
    int samp_size,
    int num_items,
    arma::vec num_responses,
    int num_quadpts,
    int num_predictors
){

  ///////////////////////////////////////
  /////////////// E-STEP ////////////////
  ///////////////////////////////////////

  List traceline(num_items);
  List etable(num_items);
  arma::mat etable_all = arma::zeros(samp_size,num_quadpts);
  List p_working(num_items+2);
  for(int k = 0; k < (num_items+2); k++){
    arma::vec p_vec = p[k];
    p_working[k] = p_vec;
  }

  arma::vec p_alpha = p_working[num_items];
  arma::vec p_phi = p_working[num_items+1];
  arma::vec alpha_impact = mean_predictors * p_alpha;
  arma::vec phi_impact = exp(var_predictors * p_phi);

  for(int j = 0; j < num_items; j++){
    traceline[j] = bernoulli_traceline_est(p_working[j],theta,predictors,samp_size,num_quadpts);
    etable[j] = arma::zeros(samp_size,num_quadpts);
  }

  arma::vec posterior(num_quadpts);

  for(int i = 0; i < samp_size; i++){
    for(int q = 0; q < num_quadpts; q++){
     posterior[q] = R::dnorm(theta[q], alpha_impact[i], sqrt(phi_impact[i]), FALSE);
    }


    for(int j = 0; j < num_items; j++){
      double resp = responses(i,j);
      arma::mat item_traceline = traceline[j];
      if(resp == 1){
        posterior = posterior % (1 - item_traceline.row(i).t());
      } else if(resp == 2){
        posterior = posterior % item_traceline.row(i).t();
      } else{
        posterior = posterior;
      }
    }

    double marginal = sum(posterior);
    posterior = posterior/marginal;
    etable_all.row(i) = posterior.t();

    for(int j = 0; j < num_items; j++){
      arma::mat item_etable = etable[j];
      item_etable.row(i) = item_etable.row(i) + posterior.t();
      etable[j] = item_etable;
    }

  }

  ///////////////////////////////////////
  /////////////// M-STEP ////////////////
  ///////////////////////////////////////


  /////////////// IMPACT ////////////////

  List anl_deriv_alpha(1);
  for(int c = 0; c < num_predictors; c++){
    anl_deriv_alpha = d_alpha_est(p_alpha,p_phi,etable_all,theta,mean_predictors,var_predictors,c,samp_size,num_items,num_quadpts);
    double deriv1 = anl_deriv_alpha[0];
    double deriv2 = anl_deriv_alpha[1];
    double p_new = p_alpha[c] - deriv1/deriv2;
    p_alpha[c] = p_new;
  }

  List anl_deriv_phi(1);
  for(int c = 0; c < num_predictors; c++){
    anl_deriv_phi = d_phi_est(p_alpha,p_phi,etable_all,theta,mean_predictors,var_predictors,c,samp_size,num_items,num_quadpts);
    double deriv1 = anl_deriv_phi[0];
    double deriv2 = anl_deriv_phi[1];
    double p_new = p_phi[c] - deriv1/deriv2;
    p_phi[c] = p_new;
  }

  p_working[num_items] = p_alpha;
  p_working[num_items+1] = p_phi;


  //////////////// DIF /////////////////

  List etable_items(num_items);

  arma::mat etable_item1 = arma::zeros(samp_size,num_quadpts);
  arma::mat etable_item2 = arma::zeros(samp_size,num_quadpts);
  arma::mat vec_zeros = arma::zeros(num_quadpts).t();

  // loop through items
  for(int j = 0; j < num_items; j++){

    arma::vec responses_item = responses.col(j);
    arma::mat etable_item = etable[j];
    etable_item1 = etable_item;
    etable_item2 = etable_item;

    for(int i = 0; i < samp_size; i++){
      if(responses_item[i] == 1){
        etable_item2.row(i) = vec_zeros;
      } else if(responses_item[i] == 2)
        etable_item1.row(i) = vec_zeros;
    }

    List anl_deriv(1);
    arma::vec p_item = p_working[j];

    //intercept base updates
      anl_deriv = d_bernoulli_est("c0",p_item,etable_item1,etable_item2,theta,predictors,0,samp_size,num_items,num_quadpts);
      double deriv1 = anl_deriv[0];
      double deriv2 = anl_deriv[1];
      double p_new = p_item[0] - deriv1/deriv2;
      p_item[0] = p_new;

    //slope base updates
    if(rasch == FALSE){
      anl_deriv = d_bernoulli_est("a0",p_item,etable_item1,etable_item2,theta,predictors,0,samp_size,num_items,num_quadpts);
      double deriv1 = anl_deriv[0];
      double deriv2 = anl_deriv[1];
      double p_new = p_item[1] - deriv1/deriv2;
      p_item[1] = p_new;
    }

    int item_j = j+1;
    if(!any(item_j == anchor)){

    //intercept dif updates
    for(int c = 0; c < num_predictors; c++){
      anl_deriv = d_bernoulli_est("c1",p_item,etable_item1,etable_item2,theta,predictors,c,samp_size,num_items,num_quadpts);
      double deriv1 = anl_deriv[0];
      double deriv2 = anl_deriv[1];
      double z = p_item[2+c] - deriv1/deriv2;
      double p_new;
      if(penalty[0] == "lasso"){
        p_new = soft_thresh_est(z,alpha,lambda[pen-1]);
      } else{
        p_new = firm_thresh_est(z,alpha,lambda[pen-1],gamma);
      }
      p_item[2+c] = p_new;
    }

    //slope dif updates
    for(int c = 0; c < num_predictors; c++){
      anl_deriv = d_bernoulli_est("a1",p_item,etable_item1,etable_item2,theta,predictors,c,samp_size,num_items,num_quadpts);
      double deriv1 = anl_deriv[0];
      double deriv2 = anl_deriv[1];
      double z = p_item[2+num_predictors+c] - deriv1/deriv2;
      double p_new;
      if(penalty[0] == "lasso"){
        p_new = soft_thresh_est(z,alpha,lambda[pen-1]);
      } else{
        p_new = firm_thresh_est(z,alpha,lambda[pen-1],gamma);
      }
      p_item[2+num_predictors+c] = p_new;
    }

    }


    p_working[j] = p_item;
  }


  return(p_working);
}

//update <- em_step(p,theta,responses,predictors,itemtypes,penalty,lambda,pen,alpha,gamma,anchor,rasch,samp_size,num_items,num_responses,num_quadpts,num_predictors)
