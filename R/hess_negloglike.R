hess_negloglike=function(thetam,alphaX,gammaX,PsiX,ti1s,ti2s,ti3s,ti4s,ti5s,ti6s,yi1s,yi2s,yi3s,Xi1s,Xi2s,Xi3s,Xi4s,Xi5s,Xi6s,Ai1s,Ai2s,Ai3s,Ai4s,Ai5s,Ai6s,tstari1s,tstari2s,tstari3s,longi.t2s,longi.t3s,NonEmptiGroups,beta_mu_vec,beta_A_vec,beta_mu_Ui_vec,beta_A_Ui_vec,knots,t_end,sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1){
  # cat("thetam",c(thetam,beta_mu_vec,beta_A_vec,beta_mu_Ui_vec,surv_par_0, surv_par_1, drop_par_0, drop_par_1,alphaX,gammaX,PsiX),"\n")
  int.tol <- 10^(-5)
  like1s<-like2s<-like3s<-like4s<-like5s<-like6s<-c()
  lsum<-0
  i <- numeric(0)
  if(1%in%NonEmptiGroups==TRUE){
    like1s<-foreach (i = seq_len(length(ti1s)),
                     .combine = c) %dopar% {
      # cat("like1",i)
      like1_value<-Egudenintg1hesscpp(thetam,alphaX,gammaX,PsiX,ti=ti1s[i],yi=yi1s[[i]][!is.na(yi1s[[i]])],Xi=Xi1s[[i]],Ai=Ai1s[i],tstari=tstari1s[[i]],
                                    beta_mu_vec = beta_mu_vec,beta_A_vec = beta_A_vec, beta_mu_Ui_vec, beta_A_Ui_vec, knots = knots, t_end = t_end,
                                    sty_end = sty_end, surv_knots = surv_knots, drop_knots=drop_knots, surv_par_0 = surv_par_0, surv_par_1 = surv_par_1, drop_par_0 = drop_par_0, drop_par_1 = drop_par_1, 
                                    tol = int.tol,nodes = nodes,weights = weights)
      return(like1_value)

    }
    lsum=sum(like1s) 

  }

  if(2%in%NonEmptiGroups==TRUE){
    like2s<-foreach (i = seq_len(length(ti2s)),
                     .combine = c) %dopar% {
                       # cat("like1",i)
                       like2_value<-Egudenintg2hesscpp(thetam,alphaX,gammaX,PsiX,ti=ti2s[i],yi=yi2s[[i]][!is.na(yi2s[[i]])],Xi=Xi2s[[i]],Ai=Ai2s[i],tstari=tstari2s[[i]],longi_t=longi.t2s[[i]],
                                                       beta_mu_vec = beta_mu_vec,beta_A_vec = beta_A_vec, beta_mu_Ui_vec, beta_A_Ui_vec, knots = knots, t_end = t_end,
                                                       sty_end = sty_end, surv_knots = surv_knots, drop_knots=drop_knots, surv_par_0 = surv_par_0, surv_par_1 = surv_par_1, drop_par_0 = drop_par_0, drop_par_1 = drop_par_1, 
                                                       tol = int.tol,nodes = nodes,weights = weights)
                       return(like2_value)
                       
                     }
    
    lsum=lsum+sum(like2s)

  }

  if(3%in%NonEmptiGroups==TRUE){
    like3s<-foreach (i = seq_len(length(ti3s)),
                     .combine = c) %dopar% {
                       # cat("like1",i)
                       like3_value<-Egudenintg3hesscpp(thetam,alphaX,gammaX,PsiX,ti=ti3s[i],yi=yi3s[[i]][!is.na(yi3s[[i]])],Xi=Xi3s[[i]],Ai=Ai3s[i],tstari=tstari3s[[i]],longi_t=longi.t3s[[i]],
                                                       beta_mu_vec = beta_mu_vec,beta_A_vec = beta_A_vec, beta_mu_Ui_vec, beta_A_Ui_vec, knots = knots, t_end = t_end,
                                                       sty_end = sty_end, surv_knots = surv_knots, drop_knots=drop_knots, surv_par_0 = surv_par_0, surv_par_1 = surv_par_1, drop_par_0 = drop_par_0, drop_par_1 = drop_par_1, 
                                                       tol = int.tol,nodes = nodes,weights = weights)
                       return(like3_value)
                       
                     }
    
    
    lsum=lsum+sum(like3s)
    # print("like3s")
    # print(sum(like3s))
  }
  if(4%in%NonEmptiGroups==TRUE){
    for (i in 1:length(ti4s)) {
      #cat("like4",i)
      like4s[i] <- Egudenintg4hesscpp(thetam,alphaX,gammaX,PsiX,ti=ti4s[i],Xi=Xi4s[[i]],Ai=Ai4s[i],
                                      sty_end = sty_end, surv_knots = surv_knots, drop_knots=drop_knots, surv_par_0 = surv_par_0, surv_par_1 = surv_par_1, drop_par_0 = drop_par_0, drop_par_1 = drop_par_1, 
                                      tol = int.tol,nodes = nodes,weights = weights)
      #cat(like4s[i])
    }
    lsum=lsum+sum(like4s)
    # print("like4s")
    # print(sum(like4s))
  }
  if(5%in%NonEmptiGroups==TRUE){
    for (i in 1:length(ti5s)) {
      #cat("like5",i)
      like5s[i] <- Egudenintg5hesscpp(thetam,alphaX,gammaX,PsiX,ti=ti5s[i],Xi=Xi5s[[i]],Ai=Ai5s[i],
                                      sty_end = sty_end, surv_knots = surv_knots, drop_knots=drop_knots, surv_par_0 = surv_par_0, surv_par_1 = surv_par_1, drop_par_0 = drop_par_0, drop_par_1 = drop_par_1, 
                                      tol = int.tol,nodes = nodes,weights = weights)
      #      cat(like5s[i])
    }
    lsum=lsum+sum(like5s)
    # print("like5s")
    # print(sum(like5s))
  }
  if(6%in%NonEmptiGroups==TRUE){
    for (i in 1:length(ti6s)) {
      #cat("like6",i)
      like6s[i] <- Egudenintg6hesscpp(thetam,alphaX,gammaX,PsiX,ti=ti6s[i],Xi=Xi6s[[i]],Ai=Ai6s[i],
                                      sty_end = sty_end, surv_knots = surv_knots, drop_knots=drop_knots, surv_par_0 = surv_par_0, surv_par_1 = surv_par_1, drop_par_0 = drop_par_0, drop_par_1 = drop_par_1, 
                                      tol = int.tol,nodes = nodes,weights = weights)
      #cat(like6s[i])
    }
    lsum=lsum+sum(like6s)
    # print("like6s")
    # print(sum(like6s))
  }
  # cat("lsum",lsum,"\n")
  return(-lsum)
  
}
