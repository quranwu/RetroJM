piecewise_time<-function(time_use,knots) {
  time_mat<-matrix(nrow=length(time_use),ncol=length(knots)-1)
  for (i in 1:length(knots)-1) {
    time_mat[,i]<-ifelse(time_use<knots[i],0,ifelse(time_use<knots[i+1],time_use-knots[i],knots[i+1]-knots[i]))
  }
  #time_mat[,length(knots)]<-ifelse(obs_t<knots[length(knots)],0,knots[length(knots)])
  return(time_mat)
}
