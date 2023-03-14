output_func<-function(meta_list,
                      opt_res_list,
                      surv_drop_knots,
                      longi_inter_knots,
                      cov_var
                      
) {
  sd_cov_be_scale<-meta_list$sd_cov_be_scale
  mean_cov_be_scale<-meta_list$mean_cov_be_scale
  mean_outcome<-meta_list$mean_outcome
  sd_outcome<-meta_list$sd_outcome
  
  par_est<-opt_res_list$par
  hessian_est<-opt_res_list$hessian
  valid_par_base<-seq_len(7)
  
  
  longi_knots_num <- longi_inter_knots+1
  valid_par_beta0<-(8:(8+longi_knots_num))
  valid_par_beta1<-((max(valid_par_beta0)+1):(max(valid_par_beta0)+longi_knots_num+1))
  valid_par_betaU<-((max(valid_par_beta1)+1):(max(valid_par_beta1)+longi_knots_num+1))
  if (surv_drop_knots == 0){
    valid_par_back<-NULL
    add_par_back<-max(valid_par_betaU)+4
  } else {
    valid_par_back<-((max(valid_par_betaU)+1):(max(valid_par_betaU)+4*surv_drop_knots))
    add_par_back<-max(valid_par_back)
  }
  
  if (is.null(cov_var)) {
    cov_length <- 1
    valid_par_cov <- NULL
    add_par_cov <- add_par_back + 3
  } else {
    cov_length <- length(cov_var)
    valid_par_cov <- ((add_par_back+1):(add_par_back+3*cov_length))
    add_par_cov <- add_par_back+3*cov_length
  }
  
  se_hess<-sqrt(diag(solve(hessian_est[c(valid_par_base,
                               valid_par_beta0,
                               valid_par_beta1,
                               valid_par_betaU,
                               valid_par_back,
                               valid_par_cov),
                             c(valid_par_base,
                               valid_par_beta0,
                               valid_par_beta1,
                               valid_par_betaU,
                               valid_par_back,
                               valid_par_cov)])))
  par_use<-par_est[c(valid_par_base,
                     valid_par_beta0,
                     valid_par_beta1,
                     valid_par_betaU,
                     valid_par_back,
                     valid_par_cov)]
  
  if (surv_drop_knots==0) {
    surv_par<-par_use[c(1:3)]
    surv_hess<-se_hess[c(1:3)]
    
    drop_par<-par_use[c(5:7)]
    drop_hess<-se_hess[c(5:7)]
    
    surv_name<-c("eita10","eita1A","alpha1U")
    drop_name<-c("gamma10","gamma1A","alpha2U")
  } else {
    surv_par<-par_use[c(1,(max(valid_par_betaU)+1):(max(valid_par_betaU)+surv_drop_knots),
                        2,(max(valid_par_betaU)+surv_drop_knots+1):(max(valid_par_betaU)+2*surv_drop_knots),
                        3)]
    surv_se<-se_hess[c(1,(max(valid_par_betaU)+1):(max(valid_par_betaU)+surv_drop_knots),
                         2,(max(valid_par_betaU)+surv_drop_knots+1):(max(valid_par_betaU)+2*surv_drop_knots),
                         3)]
    
    drop_par<-par_use[c(5,(max(valid_par_betaU)+2*surv_drop_knots+1):(max(valid_par_betaU)+3*surv_drop_knots),
                        6,(max(valid_par_betaU)+3*surv_drop_knots+1):(max(valid_par_betaU)+4*surv_drop_knots),
                        7)]
    drop_se<-se_hess[c(5,(max(valid_par_betaU)+2*surv_drop_knots+1):(max(valid_par_betaU)+3*surv_drop_knots),
                         6,(max(valid_par_betaU)+3*surv_drop_knots+1):(max(valid_par_betaU)+4*surv_drop_knots),
                         7)]
    
    surv_name<-c(paste0("eita",1:(surv_drop_knots+1),"0"),
                 paste0("eita",1:(surv_drop_knots+1),"A"),
                 "alpha1U")
    drop_name<-c(paste0("gamma",1:(surv_drop_knots+1),"0"),
                 paste0("gamma",1:(surv_drop_knots+1),"A"),
                 "alpha2U")
  }
  surv_out_dat<-data.frame(Estimate=surv_par,
                           SE=surv_se)
  rownames(surv_out_dat)<-surv_name
  
  drop_out_dat<-data.frame(Estimate=drop_par,
                           SE=drop_se)
  rownames(drop_out_dat)<-drop_name
  
  tau_par<-par_use[4] * sd_outcome
  tau_se<-se_hess[4] * sd_outcome
  tau_name<-"tau"
  
  betamu_par<-par_use[valid_par_beta0] * sd_outcome
  betamu_par[1] <- betamu_par[1] + mean_outcome
  betamu_se<-se_hess[valid_par_beta0] * sd_outcome
  betamu_name<-c("betamu_base", paste0("betamu_slope",1:(longi_inter_knots+1)))
  
  betaA_par<-par_use[valid_par_beta1] * sd_outcome
  betaA_se<-se_hess[valid_par_beta1] * sd_outcome
  betaA_name<-c("betaA_base", paste0("betaA_slope",1:(longi_inter_knots+1)))
  
  longiout_dat<-data.frame(Estimate=c(betamu_par,betaA_par,tau_par),
                           SE=c(betamu_se,betaA_se,tau_se))
  rownames(longiout_dat)<-c(betamu_name,betaA_name,tau_name)
  
  if (!is.null(cov_var)) {
    surv_X_var<-par_use[(length(par_use)-cov_length*3+1):(length(par_use)-cov_length*2)]
    surv_X_var<-surv_X_var/sd_cov_be_scale
    
    surv_X_se<-se_hess[(length(par_use)-cov_length*3+1):(length(par_use)-cov_length*2)]
    surv_X_se<-surv_X_se/sd_cov_be_scale
    
    surv_X_name<-cov_var
    surv_X_dat<-data.frame(Estimate=surv_X_var,
                           SE=surv_X_se)
    rownames(surv_X_dat)<-surv_X_name
    
    surv_out_dat$Estimate[1:(surv_drop_knots+1)]<-
      surv_out_dat$Estimate[1:(surv_drop_knots+1)] - sum(surv_X_var*mean_cov_be_scale)
    surv_out_dat<-rbind(surv_out_dat,surv_X_dat)
    
    drop_X_var<-par_use[(length(par_use)-cov_length*2+1):(length(par_use)-cov_length*1)]
    drop_X_var<-drop_X_var/sd_cov_be_scale
    
    drop_X_se<-se_hess[(length(par_use)-cov_length*2+1):(length(par_use)-cov_length*1)]
    drop_X_se<-drop_X_se/sd_cov_be_scale
    
    drop_X_name<-cov_var
    drop_X_dat<-data.frame(Estimate=drop_X_var,
                           SE=drop_X_se)
    rownames(drop_X_dat)<-drop_X_name
    
    drop_out_dat$Estimate[1:(surv_drop_knots+1)]<-
      drop_out_dat$Estimate[1:(surv_drop_knots+1)] - sum(drop_X_var*mean_cov_be_scale)
    drop_out_dat<-rbind(drop_out_dat,drop_X_dat)
    
    beta_X_var<-par_use[(length(par_use)-cov_length+1):(length(par_use))] * sd_outcome/sd_cov_be_scale
    beta_X_se<-se_hess[(length(par_use)-cov_length+1):(length(par_use))] * sd_outcome/sd_cov_be_scale

    beta_X_name<-cov_var
    beta_X_dat<-data.frame(Estimate=beta_X_var,
                           SE=beta_X_se)
    rownames(beta_X_dat)<-beta_X_name
    longiout_dat$Estimate[1] <- longiout_dat$Estimate[1] - sum(beta_X_var*mean_cov_be_scale)
    longiout_dat<-rbind(longiout_dat,beta_X_dat)
    
  }
  
  surv_out_dat$CI.low<-surv_out_dat$Estimate-1.96*surv_out_dat$SE
  surv_out_dat$CI.up<-surv_out_dat$Estimate+1.96*surv_out_dat$SE
  surv_out_dat$p.value<-(1-pnorm(abs(surv_out_dat$Estimate/surv_out_dat$SE)))*2
  
  drop_out_dat$CI.low<-drop_out_dat$Estimate-1.96*drop_out_dat$SE
  drop_out_dat$CI.up<-drop_out_dat$Estimate+1.96*drop_out_dat$SE
  drop_out_dat$p.value<-(1-pnorm(abs(drop_out_dat$Estimate/drop_out_dat$SE)))*2
  
  longiout_dat$CI.low<-longiout_dat$Estimate-1.96*longiout_dat$SE
  longiout_dat$CI.up<-longiout_dat$Estimate+1.96*longiout_dat$SE
  longiout_dat$p.value<-(1-pnorm(abs(longiout_dat$Estimate/longiout_dat$SE)))*2
  
  return(list(surv_out_dat=surv_out_dat,
              drop_out_dat=drop_out_dat,
              longiout_dat=longiout_dat))
}



