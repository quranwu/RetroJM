run_func <- function(meta_list,
                     surv_drop_knots,
                     longi_inter_knots,
                     Ai_var,
                     last_obs_time_var,
                     cov_var,
                     ID_var,
                     num_cores) {
  ini_value_list <- ini_value(
    meta_list = meta_list,
    surv_drop_knots = surv_drop_knots,
    longi_inter_knots = longi_inter_knots,
    Ai_var = Ai_var,
    last_obs_time_var = last_obs_time_var,
    cov_var = cov_var,
    ID_var = ID_var
  )
  ini_value_use <- ini_value_list$ini.theta
  lb_use <- ini_value_list$lb
  ub_use <- ini_value_list$ub
  longi_knots_use <- ini_value_list$longi_knots_use
  
  Aeq<-matrix(NA,ncol = length(ini_value_use),nrow=longi_inter_knots+1)
  Aeq_vec<-rep(0,length(ini_value_use))
  for (i in 1:(longi_inter_knots+1)) {
    Aeq_vec[c(8+2*longi_inter_knots+4)]<-1
    Aeq_vec[c(8+2*longi_inter_knots+4+i)]<-longi_knots_use[i+1]
    Aeq[i,]<- -Aeq_vec
  }
  
  Beq<-rep(0,longi_inter_knots+1)
  
  ti1s <- meta_list$ti1s
  yi1s <- meta_list$yi1s
  Ai1s <- meta_list$Ai1s
  Xi1s <- meta_list$Xi1s
  tstari1s <- meta_list$tstari1s
  ti2s <- meta_list$ti2s
  yi2s <- meta_list$yi2s
  Ai2s <- meta_list$Ai2s
  Xi2s <- meta_list$Xi2s
  tstari2s <- meta_list$tstari2s
  longi.t2s <- meta_list$longi.t2s
  ti3s <- meta_list$ti3s
  yi3s <- meta_list$yi3s
  Ai3s <- meta_list$Ai3s
  Xi3s <- meta_list$Xi3s
  tstari3s <- meta_list$tstari3s
  longi.t3s <- meta_list$longi.t3s
  ti4s <- meta_list$ti4s
  Ai4s <- meta_list$Ai4s
  Xi4s <- meta_list$Xi4s
  ti5s <- meta_list$ti5s
  Ai5s <- meta_list$Ai5s
  Xi5s <- meta_list$Xi5s
  ti6s <- meta_list$ti6s
  Ai6s <- meta_list$Ai6s
  Xi6s <- meta_list$Xi6s
  longi_end_time <- meta_list$longi_end_time
  sty_end_time <- meta_list$sty_end_time
  
  sub_group<-c(length(ti1s),length(ti2s),length(ti3s),length(ti4s),length(ti5s),length(ti6s))
  NonEmptiGroups<-which(sub_group!=0)
  
  
  surv_drop_quantile <- seq(from=0,to=1,length.out=surv_drop_knots+2)
  
  if (surv_drop_knots == 0) {
    surv_knots <- sty_end_time
    drop_knots <- sty_end_time
  } else {
    surv_knots <- quantile(ti1s,surv_drop_quantile)[-c(1,length(surv_drop_quantile))]
    drop_knots <- quantile(ti2s,surv_drop_quantile)[-c(1,length(surv_drop_quantile))]
  }
  
  cov_length<-ifelse(is.null(cov_var),1,length(cov_var))
  
  registerDoParallel(num_cores)
  
  est1<-try(solnl(ini_value_use, objfun =  function(parameters){
    par<-numeric(9)
    par[1]<-1
    par[2]<-1
    par[3]<-parameters[1]
    par[4]<-parameters[2]
    par[5]<-parameters[3]
    par[6]<-parameters[4]
    par[7]<- parameters[5]
    par[8]<- parameters[6]
    par[9]<-parameters[7]
    # par[10]<-parameters[8]
    
    beta_mu_vec<-parameters[8:(8+longi_inter_knots+1)]
    beta_A_vec<-parameters[(8+longi_inter_knots+2):(8+2*longi_inter_knots+3)]
    beta_mu_Ui_vec<-parameters[(8+2*longi_inter_knots+4):(8+3*longi_inter_knots+5)]
    beta_A_Ui_vec<-rep(0,longi_inter_knots+2)
    
    surv_par_0<-parameters[(14+3*longi_inter_knots):
                             (13+3*longi_inter_knots+surv_drop_knots)]
    surv_par_1<-parameters[(14+3*longi_inter_knots+surv_drop_knots):
                             (13+3*longi_inter_knots+2*surv_drop_knots)]
    drop_par_0<-parameters[(14+3*longi_inter_knots+2*surv_drop_knots):
                             (13+3*longi_inter_knots+3*surv_drop_knots)]
    drop_par_1<-parameters[(14+3*longi_inter_knots+3*surv_drop_knots):
                             (13+3*longi_inter_knots+4*surv_drop_knots)]
    
    alphaX<-parameters[(14+3*longi_inter_knots+4*surv_drop_knots):
                         (13+3*longi_inter_knots+4*surv_drop_knots+cov_length)]
    
    gammaX<-parameters[(14+3*longi_inter_knots+4*surv_drop_knots+cov_length):
                         (13+3*longi_inter_knots+4*surv_drop_knots+2*cov_length)]
    
    PsiX<-parameters[(14+3*longi_inter_knots+4*surv_drop_knots+2*cov_length):
                       (13+3*longi_inter_knots+4*surv_drop_knots+3*cov_length)]
    
    
    hess_negloglike(par,alphaX,gammaX,PsiX, ti1s,ti2s,ti3s,ti4s,ti5s,ti6s,yi1s,yi2s,yi3s,Xi1s,Xi2s,Xi3s,Xi4s,Xi5s,Xi6s,Ai1s,Ai2s,Ai3s,Ai4s,Ai5s,Ai6s,
                    tstari1s,tstari2s,tstari3s,longi.t2s,longi.t3s,NonEmptiGroups,beta_mu_vec,beta_A_vec,beta_mu_Ui_vec,beta_A_Ui_vec,
                    longi_knots_use,longi_end_time,sty_end_time,surv_knots,drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1)},
    A = Aeq,
    B = Beq,
    lb = lb_use,
    ub = ub_use
  ),TRUE)
  stopImplicitCluster()
  
  est1$longi_knots_use <- longi_knots_use
  
  return(est1)
}
