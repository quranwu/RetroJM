ini_value<-function(meta_list,
                    surv_drop_knots,
                    longi_inter_knots,
                    Ai_var,
                    last_obs_time_var,
                    cov_var,
                    ID_var
                    ) {
  longi_dat_cov_comp <- meta_list$longi_dat_cov_comp
  tstari1s <- meta_list$tstari1s
  longi_end_time <- meta_list$longi_end_time
  sty_end_time <- meta_list$sty_end_time
  
  if (longi_inter_knots==0) {
    longi_knots_use<-c(0,longi_end_time)
  } else {
    quan_len<-seq(from=0,to=1,length.out = longi_inter_knots+2)
    longi_knots<-as.numeric(quantile(unlist(tstari1s),quan_len))
    longi_knots_use<-c(0,longi_knots[2:(longi_inter_knots+1)],longi_end_time)
  }
  
  
  ini_base<-c(0,0,0,sd(longi_dat_cov_comp$outcome),0,0,0)
  lb_base<-c(rep(-5,3),0.1,rep(-5,3))
  ub_base<-c(rep(5,3),30,rep(5,3))
  
  if (surv_drop_knots==0) {
    ini_base_back<-c(rep(0,4))
    lb_base_back<-c(rep(0,4))
    ub_base_back<-c(rep(0,4))
  } else {
    ini_base_back<-c(rep(0,surv_drop_knots*4))
    lb_base_back<-rep(-5,surv_drop_knots*4)
    ub_base_back<-rep(5,surv_drop_knots*4)
  }
  
  lb_ini_mu<-c(rep(-100,longi_inter_knots+2))
  lb_ini_A<-c(rep(-100,longi_inter_knots+2))
  lb_ini_mu_Ui<-c(rep(-100,longi_inter_knots+2))
  
  ub_ini_mu<-c(rep(100,longi_inter_knots+2))
  ub_ini_A<-c(rep(100,longi_inter_knots+2))
  ub_ini_mu_Ui<-c(rep(100,longi_inter_knots+2))
  
  t_star_dat_mod<-piecewise_time(longi_dat_cov_comp$longi_tstar,longi_knots_use)
  
  if (is.null(cov_var)) {
    ini_res<-lmer(longi_dat_cov_comp$outcome ~ as.matrix(t_star_dat_mod) + 
                    (1|ID_new) + 
                    longi_dat_cov_comp[[Ai_var]] + 
                    as.matrix(longi_dat_cov_comp[[Ai_var]] * t_star_dat_mod),
                  data=longi_dat_cov_comp)
    
    ini_mu<-as.vector(fixef(ini_res)[1:(longi_inter_knots+2)])
    ini_A<-as.vector(fixef(ini_res)[(longi_inter_knots+3):(2*longi_inter_knots+4)])
    ini_mu_Ui<-rep(0,longi_inter_knots+2)
    ini_alphaX<-0
    ini_gammaX<-0
    ini_PsiX<-0
    
    lb_alphaX<-0
    lb_gammaX<-0
    lb_PsiX<-0
    
    ub_alphaX<-0
    ub_gammaX<-0
    ub_PsiX<-0
  } else {
    cov_length<-length(cov_var)
    ini_res<-lmer(longi_dat_cov_comp$outcome ~ as.matrix(t_star_dat_mod) + 
                    (1|ID_new) + 
                    longi_dat_cov_comp[[Ai_var]] + 
                    as.matrix(longi_dat_cov_comp[[Ai_var]] * t_star_dat_mod) + 
                    as.matrix(longi_dat_cov_comp[,cov_var,drop=FALSE]),
                  data=longi_dat_cov_comp)
    
    ini_mu <- as.vector(fixef(ini_res)[1:(longi_inter_knots + 2)])
    ini_A <-
      as.vector(fixef(ini_res)[(longi_inter_knots + 3):(2 * longi_inter_knots +
                                                          4)])
    ini_mu_Ui <- rep(0, longi_inter_knots + 2)
    ini_alphaX <- rep(0, cov_length)
    ini_gammaX <- rep(0, cov_length)
    ini_PsiX <-
      as.vector(fixef(ini_res)[(2 * longi_inter_knots + 5):(2 * longi_inter_knots +
                                                              5 + cov_length - 1)])
    
    lb_alphaX<- rep(-5, cov_length)
    lb_gammaX<- rep(-5, cov_length)
    lb_PsiX<- rep(-100, cov_length)
    
    ub_alphaX<- rep(5, cov_length)
    ub_gammaX<- rep(5, cov_length)
    ub_PsiX<- rep(100, cov_length)
    
  }
  
  ini.theta <-
    c(ini_base,
      ini_mu,
      ini_A,
      ini_mu_Ui,
      ini_base_back,
      ini_alphaX,
      ini_gammaX,
      ini_PsiX)
  
  lb<-c(lb_base,lb_ini_mu,lb_ini_A,lb_ini_mu_Ui,lb_base_back,lb_alphaX,lb_gammaX,lb_PsiX)
  ub<-c(ub_base,ub_ini_mu,ub_ini_A,ub_ini_mu_Ui,ub_base_back,ub_alphaX,ub_gammaX,ub_PsiX)
  
  return(list(ini.theta=ini.theta,
              lb=lb,
              ub=ub,
              longi_knots_use = longi_knots_use))
  
}
