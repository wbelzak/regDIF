// #include <RcppArmadillo.h>
// //[[Rcpp::depends(RcppArmadillo)]]
// using namespace Rcpp;
// using std::exp;
// using std::sqrt;
// using std::pow;



////////////////////////////////////////////
/////////////// TRESHOLDING ////////////////
////////////////////////////////////////////

// template <typename T> double sgn(T val) {
//   return (T(0) < val) - (val < T(0));
// }
// template <class T> const T& max (const T& a, const T& b) {
//   return (a<b)?b:a;     // or: return comp(a,b)?b:a; for version (2)
// }
//
// //[[Rcpp::export]]
// double soft_thresh_cpp(
//     double z,
//     double alpha,
//     double tau
// ){
//
//   double t = (fabs(z) - tau*alpha)/(1+tau*(1-alpha));
//   double p_new = sgn(z)*max(t,0.0);
//   return p_new;
//
// }
//
// //[[Rcpp::export]]
// double firm_thresh_cpp(
//     double z,
//     double alpha,
//     double tau,
//     double gamma
// ){
//
//   double p_new;
//   if(fabs(z)/(1+tau*(1-alpha)) <= gamma*tau){
//     p_new = (gamma/(gamma-1))*soft_thresh_cpp(z,alpha,tau);
//   }else{
//     p_new = z/(1+tau*(1-alpha));
//   }
//
//   return(p_new);
// }
//
// //[[Rcpp::export]]
// arma::vec dnormCpp(
//     arma::vec theta,
//     double alpha_impact,
//     double phi_impact,
//     int num_quadpts
// ){
//   arma::vec posterior(num_quadpts);
//   static const double pi = 3.14159265;
//
//   for(int q = 0; q < num_quadpts; q++){
//     posterior[q] = 1 / sqrt(2*pi*phi_impact)*exp(-pow((theta[q] -
//       alpha_impact),2) / (2*phi_impact));
//   }
//
//   return posterior;
// }
//
// ///////////////////////////////////////////
// /////////////// TRACELINES ////////////////
// ///////////////////////////////////////////
//
// //[[Rcpp::export]]
// arma::mat bernoulli_traceline_cpp(
//     arma::vec p_item,
//     arma::vec theta,
//     arma::mat predictors,
//     int samp_size,
//     int num_quadpts
// ){
//
//   // split item parameter vector
//   double p_c0 = p_item(0,0);
//   double p_a0 = p_item(1,0);
//   arma::vec p_c1 = p_item.subvec(2,1+predictors.n_cols);
//   arma::vec p_a1 = p_item.subvec(2+predictors.n_cols,2*predictors.n_cols+1);
//
//   // create matrix to hold traceline
//   arma::mat traceline = arma::zeros(samp_size,num_quadpts);
//
//   // Loop through quadrature points
//   for(int q = 0; q < num_quadpts; ++q){
//     traceline.col(q) = 1/(1 + exp(-(p_c0 + predictors * p_c1 +
//       (p_a0 + predictors * p_a1) * theta[q])));
//   }
//
//   return traceline;
// }
//
//
// //[[Rcpp::export]]
// List categorical_traceline_cpp(
//     arma::vec p_item,
//     arma::vec theta,
//     arma::mat predictors,
//     int samp_size,
//     int num_responses_item,
//     int num_quadpts
// ){
//
//   // split item parameter vector
//   double p_c0_int = p_item(0,0);
//   arma::vec p_c0_thr = p_item.subvec(1,num_responses_item - 2);
//   double p_a0 = p_item(num_responses_item - 1,0);
//   arma::vec p_c1 = p_item.subvec(num_responses_item,
//                                  num_responses_item + predictors.n_cols - 1);
//   arma::vec p_a1 = p_item.subvec(num_responses_item + predictors.n_cols,
//                                  num_responses_item + predictors.n_cols * 2 -
//                                    1);
//
//   // Create matrices to hold tracelines
//   List traceline_list(num_responses_item);
//   for(int i = 0; i < num_responses_item; ++i){
//     traceline_list[i] = arma::zeros(samp_size,num_quadpts);
//   }
//
//   // Item Response 1
//   arma::mat first_response = traceline_list[0];
//   for(int q = 0; q < num_quadpts; ++q){
//     first_response.col(q) =
//       1 - 1 / (1 + exp(-(p_c0_int + predictors * p_c1 +
//       (p_a0 + predictors * p_a1) * theta[q])));
//   }
//   traceline_list[0] = first_response;
//
//   // Item Response J
//   arma::mat last_response = traceline_list[num_responses_item - 1];
//   for(int q = 0; q < num_quadpts; ++q){
//     last_response.col(q) =
//       1 / (1 + exp(-(p_c0_int - p_c0_thr(num_responses_item - 3,0) +
//       predictors * p_c1 + (p_a0 + predictors * p_a1) * theta[q])));
//   }
//   traceline_list[num_responses_item - 1] = last_response;
//
//   if(num_responses_item > 2){
//
//     // Item Response 2 (if more than 2 responses)
//     arma::mat second_response = traceline_list[1];
//     for(int q = 0; q < num_quadpts; ++q){
//       second_response.col(q) =
//         1 / (1 + exp(-(p_c0_int + predictors * p_c1 +
//         (p_a0 + predictors * p_a1) * theta[q]))) -
//         1 / (1 + exp(-(p_c0_int - p_c0_thr(0,0) + predictors * p_c1 +
//         (p_a0 + predictors * p_a1) * theta[q])));
//     }
//     traceline_list[1] = second_response;
//
//     if(num_responses_item > 3)
//
//       // Item Response 3 to J-1 Response (if more than 3 responses)
//       for(int r = 2; r < num_responses_item - 1; r++){
//
//         arma::mat next_response = traceline_list[r];
//
//         for(int q = 0; q < num_quadpts; q++){
//           next_response.col(q) =
//             1 / (1 + exp(-(p_c0_int - p_c0_thr(r-2,0) + predictors * p_c1 +
//             (p_a0 + predictors * p_a1) * theta[q]))) -
//             1 / (1 + exp(-(p_c0_int - p_c0_thr(r-1,0) + predictors * p_c1 +
//             (p_a0 + predictors * p_a1) * theta[q])));
//         }
//
//         traceline_list[r] = next_response;
//       }
//   }
//   return traceline_list;
// }
//
//
// //[[Rcpp::export]]
// List cumulative_traceline_cpp(
//     arma::vec p_item,
//     arma::vec theta,
//     arma::mat predictors,
//     int samp_size,
//     int num_responses_item,
//     int num_quadpts
// ){
//
//   // split item parameter vector
//   double p_c0_int = p_item(0,0);
//   arma::vec p_c0_thr = p_item.subvec(1,num_responses_item - 2);
//   double p_a0 = p_item(num_responses_item - 1,0);
//   arma::vec p_c1 = p_item.subvec(num_responses_item,
//                                  num_responses_item + predictors.n_cols - 1);
//   arma::vec p_a1 = p_item.subvec(num_responses_item + predictors.n_cols,
//                                  num_responses_item + predictors.n_cols * 2 -
//                                    1);
//
//   // Create matrices to hold tracelines
//   List traceline_list(num_responses_item-1);
//   for(int i = 0; i < num_responses_item-1; ++i){
//     traceline_list[i] = arma::zeros(samp_size,num_quadpts);
//   }
//
//   // Item Response 1
//   arma::mat first_response = traceline_list[0];
//   for(int q = 0; q < num_quadpts; ++q){
//     first_response.col(q) =
//       1 / (1 + exp(-(p_c0_int + predictors * p_c1 +
//       (p_a0 + predictors * p_a1) * theta[q])));
//   }
//   traceline_list[0] = first_response;
//
//   if(num_responses_item > 2){
//
//     // Item Response 2 to J-1 Response
//     for(int r = 1; r < num_responses_item - 1; r++){
//
//       arma::mat next_response = traceline_list[r];
//
//       for(int q = 0; q < num_quadpts; q++){
//         next_response.col(q) =
//           1 / (1 + exp(-(p_c0_int - p_c0_thr(r-1,0) + predictors * p_c1 +
//           (p_a0 + predictors * p_a1) * theta[q])));
//       }
//
//       traceline_list[r] = next_response;
//     }
//   }
//   return traceline_list;
// }
//
// ////////////////////////////////////////////
// /////////////// DERIVATIVES ////////////////
// ////////////////////////////////////////////
//
//
// /////////////// ALPHA - IMPACT ////////////////
//
// //[[Rcpp::export]]
// List d_alpha_cpp(
//     arma::vec p_alpha,
//     arma::vec p_phi,
//     arma::mat etable,
//     arma::vec theta,
//     arma::mat mean_predictors,
//     arma::mat var_predictors,
//     int cov,
//     int samp_size,
//     int num_items,
//     int num_quadpts
// ){
//
//   // obtain latent mean and variance vectors
//   arma::vec alpha_impact = mean_predictors * p_alpha;
//   arma::vec phi_impact = exp(var_predictors * p_phi);
//
//   arma::mat eta_d = arma::zeros(samp_size,num_quadpts);
//   for(int q = 0; q < num_quadpts; q++){
//     eta_d.col(q) = mean_predictors.col(cov);
//   }
//
//
//   arma::mat d1_trace = arma::zeros(samp_size,num_quadpts);
//   arma::mat d2_trace = arma::zeros(samp_size,num_quadpts);
//   for(int i = 0; i < samp_size; i++){
//     d1_trace.row(i) = eta_d.row(i)/phi_impact[i] % (theta.t() -
//       alpha_impact[i]);
//     d2_trace.row(i) = -pow(eta_d.row(i),2)/phi_impact[i];
//   }
//
//
//
//   List dlist(2);
//   dlist[0] = sum(sum(etable % d1_trace));
//   dlist[1] = sum(sum(etable % d2_trace));
//
//   return(dlist);
//
// }
//
//
// /////////////// PHI - IMPACT ////////////////
//
// //[[Rcpp::export]]
// List d_phi_cpp(
//     arma::vec p_alpha,
//     arma::vec p_phi,
//     arma::mat etable,
//     arma::vec theta,
//     arma::mat mean_predictors,
//     arma::mat var_predictors,
//     int cov,
//     int samp_size,
//     int num_items,
//     int num_quadpts
// ){
//
//   // obtain latent mean and variance vectors
//   arma::vec alpha_impact = mean_predictors * p_alpha;
//   arma::vec phi_impact = exp(var_predictors * p_phi);
//
//   arma::vec eta_d1 = sqrt(phi_impact) / 2 % var_predictors.col(cov);
//   arma::vec eta_d2 = sqrt(phi_impact) / 4 % pow(var_predictors.col(cov),2);
//
//   arma::mat d1_trace = arma::zeros(samp_size,num_quadpts);
//   arma::mat d2_trace = arma::zeros(samp_size,num_quadpts);
//   for(int i = 0; i < samp_size; i++){
//     d1_trace.row(i) = eta_d1[i] * (pow((theta.t() - alpha_impact[i]),2) /
//       pow(phi_impact[i],(3.0/2)) - 1 / sqrt(phi_impact[i]));
//     d2_trace.row(i) = -2 * eta_d2[i] * (pow(phi_impact[i],(-3.0/2)) *
//       pow((theta.t() - alpha_impact[i]),2));
//   }
//
//   List dlist(2);
//   dlist[0] = sum(sum(etable % d1_trace));
//   dlist[1] = sum(sum(etable % d2_trace));
//
//   return(dlist);
//
// }
//
//
// /////////////// BERNOULLI - DIF ////////////////
//
// //[[Rcpp::export]]
// List d_bernoulli_cpp(
//     std::string parm,
//     arma::vec p_item,
//     arma::mat etable1,
//     arma::mat etable2,
//     arma::vec theta,
//     arma::mat predictors,
//     int cov,
//     int samp_size,
//     int num_items,
//     int num_quadpts
// ){
//
//   arma::mat eta_d = arma::zeros(samp_size,num_quadpts);
//   if(parm == "c0"){
//     eta_d.ones();
//   } else if(parm == "a0"){
//     for(int q = 0; q < num_quadpts; q++){
//       eta_d.col(q) = eta_d.col(q) + theta[q];
//     }
//   } else if(parm == "c1"){
//     for(int q = 0; q < num_quadpts; q++){
//       eta_d.col(q) = predictors.col(cov);
//     }
//   } else if(parm == "a1"){
//     for(int q = 0; q < num_quadpts; q++){
//       eta_d.col(q) = predictors.col(cov) * theta[q];
//     }
//   }
//
//
//   arma::mat traceline = bernoulli_traceline_cpp(p_item,
//                                                 theta,
//                                                 predictors,
//                                                 samp_size,
//                                                 num_quadpts);
//
//   List dlist(2);
//   dlist[0] = sum(sum((1-traceline) % eta_d % etable2)) +
//     sum(sum(-traceline % eta_d % etable1));
//   dlist[1] = sum(sum(-traceline % (1-traceline) % pow(eta_d,2) % etable1)) +
//     sum(sum(-traceline % (1-traceline) % pow(eta_d,2) % etable2));
//
//   return(dlist);
//
// }
//
// /////////////// CATEGORICAL - DIF ////////////////
//
// //[[Rcpp::export]]
// List d_categorical_cpp(
//     std::string parm,
//     arma::vec p_item,
//     Rcpp::List etable,
//     arma::vec theta,
//     arma::mat predictors,
//     int thr,
//     int cov,
//     int samp_size,
//     int num_responses_item,
//     int num_items,
//     int num_quadpts
// ){
//
//   arma::mat eta_d = arma::zeros(samp_size,num_quadpts);
//   if(parm == "c0"){
//     eta_d.ones();
//   } else if(parm == "a0"){
//     for(int q = 0; q < num_quadpts; q++){
//       eta_d.col(q) = eta_d.col(q) + theta[q];
//     }
//   } else if(parm == "c1"){
//     for(int q = 0; q < num_quadpts; q++){
//       eta_d.col(q) = predictors.col(cov-1);
//     }
//   } else if(parm == "a1"){
//     for(int q = 0; q < num_quadpts; q++){
//       eta_d.col(q) = predictors.col(cov-1) * theta[q];
//     }
//   }
//
//   List cum_traceline = cumulative_traceline_cpp(p_item,
//                                                 theta,
//                                                 predictors,
//                                                 samp_size,
//                                                 num_responses_item,
//                                                 num_quadpts);
//
//   arma::mat d1 = arma::zeros(samp_size,num_quadpts);
//   arma::mat d2 = arma::zeros(samp_size,num_quadpts);
//   double d1_sum;
//   double d2_sum;
//
//   if(thr < 0) {
//     arma::mat etable_first = etable[0];
//     arma::mat etable_last = etable[num_responses_item - 1];
//     arma::mat cum_traceline_first = cum_traceline[0];
//     arma::mat cum_traceline_last = cum_traceline[num_responses_item - 2];
//
//     arma::mat d1 = eta_d % (-1 * etable_first % cum_traceline_first +
//       etable_last % (1 - cum_traceline_last));
//     arma::mat d2 = pow(eta_d,2) % (-1 * etable_first %
//       (cum_traceline_first % (1 - cum_traceline_first)) +
//       etable_last % (-1 * cum_traceline_last % (1 - cum_traceline_last)));
//
//     for(int i = 1; i < num_responses_item - 1; i++) {
//       arma::mat etable_next = etable[i];
//       arma::mat cum_traceline_next = cum_traceline[i];
//       arma::mat cum_traceline_next_minus1 = cum_traceline[i-1];
//
//       d1 = d1 + eta_d %
//         etable_next %
//         ((1 - cum_traceline_next) -
//         cum_traceline_next_minus1);
//       d2 = d2 + pow(eta_d,2) %
//         etable_next %
//         (pow(cum_traceline_next_minus1,2) +
//         pow(cum_traceline_next,2) -
//         cum_traceline_next_minus1 -
//         cum_traceline_next);
//     }
//
//     d1_sum = arma::accu(d1);
//     d2_sum = arma::accu(d2);
//
//   } else {
//     List cat_traceline = categorical_traceline_cpp(p_item,
//                                                    theta,
//                                                    predictors,
//                                                    samp_size,
//                                                    num_responses_item,
//                                                    num_quadpts);
//
//     arma::mat etable_current_thr = etable[thr-1];
//     arma::mat etable_next_thr = etable[thr];
//     arma::mat cum_traceline_current_thr = cum_traceline[thr-1];
//     arma::mat cat_traceline_current_thr = cat_traceline[thr-1];
//     arma::mat cat_traceline_next_thr = cat_traceline[thr];
//
//     arma::mat d1 = (-1 * etable_current_thr % cum_traceline_current_thr %
//       (1 - cum_traceline_current_thr) / cat_traceline_current_thr) +
//       (etable_next_thr % cum_traceline_current_thr %
//       (1 - cum_traceline_current_thr) / cat_traceline_next_thr);
//     arma::mat d2 = etable_current_thr / cat_traceline_current_thr %
//       (cum_traceline_current_thr % pow(1 - cum_traceline_current_thr, 2) -
//       pow(cum_traceline_current_thr,2) % (1 - cum_traceline_current_thr) +
//       pow(cum_traceline_current_thr,2) % pow(1 - cum_traceline_current_thr,2) /
//         cat_traceline_current_thr) -
//           etable_next_thr / cat_traceline_next_thr %
//           (cum_traceline_current_thr % pow(1 - cum_traceline_current_thr,2) -
//           pow(cum_traceline_current_thr,2) % (1 - cum_traceline_current_thr) -
//           pow(cum_traceline_current_thr,2) %
//           pow(1 - cum_traceline_current_thr,2) / cat_traceline_next_thr);
//
//     d1_sum = arma::accu(d1);
//     d2_sum = arma::accu(d2);
//   }
//
//   List dlist(2);
//   dlist[0] = d1_sum;
//   dlist[1] = d2_sum;
//
//   return(dlist);
//
// }
//
//
