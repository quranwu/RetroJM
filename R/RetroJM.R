##' Joint model on the retrospective time scale
##'
##' @description
##' \loadmathjax
##' RetroJM is used for fitting joint model on retrospective time scale starting from death time.
##' The outcome model is:
##' \mjdeqn{Y_i(t^*)=\beta_{\mu}(t^*)+A_i\beta_A(t^*)+X_i^T\psi_X+U_i^*(t^*)+\epsilon_i(t^*)}{}
##' where \mjeqn{t^*=D-t}{}, D is death time and t is observed time. 
##' The hazard for competing risks are: 
##' \mjdeqn{\lambda_i^{\delta}(t)=\exp{\Big(\alpha_0^{\delta}(t)+A_i\alpha_A^{\delta}(t)+\alpha_U^{\delta} U_i \Big)},\hspace{0.2cm}\text{for}\hspace{0.12cm} \delta=1,2,}{}
##'
##' @param data The dataset used in analysis. It could be either long format or wide format. 
##' It should contain longitudinal measurement, measurement time, a binary variable as 
##' parameter of interest, ID for subjects, event/censoring time, an indicator
##' for death status, an indicator for drop out status, and any covariates. 
##' @param long_format  A logical argument takes value TRUE or FALSE. If TRUE, the input dataset
##' is in long format. If FALSE, the input dataset should be wide format. Defaut is TRUE.  
##' @param outcome_var   Name of outcome variables. If the input data is long format, it should 
##' be a single character; If the input is wide format, it should be a character vector. 
##' @param time_var Name of observed time variables. If the input data is long format, it should 
##' be a single character; If the input is wide format, it should be a character vector. It should 
##' have same length as "outcome_var". 
##' @param Ai_var Name of parameter of interest. Have to be a binary variable coded as 0/1, 
##' such as control/treatment.
##' @param last_obs_time_var Name of variable for event/censoring time. 
##' @param ID_var Name of ID variable. 
##' @param surv_censor_var Name of indicator variable for death status. 
##' @param drop_censor_var Name of indicator variable for drop out status. 
##' @param cov_var Name of covariates. Default is NULL, meaning no covariates. 
##' @param longi_end_time The end time for longitudinal trajectory. It take numeric as input.
##' If NULL, then it equals to study end time. Default is NULL.  
##' @param surv_drop_knots The number of additional pieces for piecewise constant hazard model in survival
##' and dropout hazard. Default is 0, meaning no additional pieces and the hazard is constant. 
##' @param longi_inter_knots The number of additional pieces for piecewise linear model in 
##' longitudinal trajectory. Default is 0, meaning no additional pieces and the longitudinal
##' trajectory is linear. 
##' @param num_cores The number of cores used for parallel computing. 
##'
##' @return
##' A `list` of `3` datasets containing the parameter estiamtes for survival, dropout, and 
##' longitudinal measurement.
##' Each dataset has row representing each parameter, 5 columns for `Estimates`, `Standard Error`,
##' `Lower bound for 95\% Confidence Interval`, `Upper bound for 95\% Confidence Interval`,
##' `P value`.
##'
##' @examples 
##' 
##' \dontrun{
##' library(RetroJM)
##'
##' ## load the simulated example data Retro_data
##' data(Retro_data)
##' 
##' test_run<-RetroJM(data = Retro_data,
##' long_format = FALSE,
##' outcome_var = paste0("Y",1:13),
##' time_var = paste0("time",1:13),
##' Ai_var = "Ai",
##' last_obs_time_var = "observed.t",
##' cov_var = paste0("x",1:2),
##' ID_var = "ID",
##' surv_censor_var = "surv_ind",
##' drop_censor_var = "drop_ind",
##' longi_end_time = 25,
##' surv_drop_knots = 1,
##' longi_inter_knots = 1,
##' num_cores = 2)
##' 
##' }
##'
##'
##' @importFrom foreach foreach %dopar%
##' @importFrom parallel makeCluster clusterExport stopCluster clusterSetRNGStream detectCores
##' @importFrom doParallel registerDoParallel stopImplicitCluster
##' @importFrom NlcOptim solnl
##' @importFrom tidyr pivot_wider pivot_longer
##' @importFrom lme4 lmer fixef
##' @import Rcpp
##' @import stats
##' @import mathjaxr
##'
##'
##'
##' @export
##' @useDynLib RetroJM, .registration=TRUE

RetroJM<-function(data,
                  long_format = TRUE,
                  outcome_var,
                  time_var,
                  Ai_var,
                  last_obs_time_var,
                  cov_var,
                  ID_var,
                  surv_censor_var,
                  drop_censor_var,
                  longi_end_time = NULL,
                  surv_drop_knots = 0,
                  longi_inter_knots = 0,
                  num_cores = detectCores() - 2
                  ) {
  if (long_format) {
    meta_list <- metaData_long(
      long_dat = data,
      longi_outcome_var = outcome_var,
      time_var = time_var,
      Ai_var = Ai_var,
      last_obs_time_var = last_obs_time_var,
      cov_var = cov_var,
      ID_var = ID_var,
      death_ind_var = surv_censor_var,
      drop_ind_var = drop_censor_var,
      longi_end_time = longi_end_time
    )
  } else {
    meta_list <- metaData_wide(
      wide_dat = data,
      longi_outcome_var = outcome_var,
      time_var = time_var,
      Ai_var = Ai_var,
      last_obs_time_var = last_obs_time_var,
      cov_var = cov_var,
      ID_var = ID_var,
      death_ind_var = surv_censor_var,
      drop_ind_var = drop_censor_var,
      longi_end_time = longi_end_time
    )
  }

  opt_res <- run_func(
    meta_list = meta_list,
    surv_drop_knots = surv_drop_knots,
    longi_inter_knots = longi_inter_knots,
    Ai_var = Ai_var,
    last_obs_time_var = last_obs_time_var,
    cov_var = cov_var,
    ID_var = ID_var,
    num_cores = num_cores
  )
  
  output_res<-output_func(
    meta_list = meta_list,
    opt_res_list = opt_res,
    surv_drop_knots = surv_drop_knots,
    longi_inter_knots = longi_inter_knots,
    cov_var = cov_var
  )
  
  output_res$opt_res<-opt_res
  
  return(output_res)
  
}
