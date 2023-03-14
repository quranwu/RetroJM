metaData_long <- function(long_dat,
                          longi_outcome_var,
                          time_var,
                          Ai_var,
                          last_obs_time_var,
                          cov_var,
                          ID_var,
                          death_ind_var,
                          drop_ind_var,
                          longi_end_time
) {
  
  if (any(length(longi_outcome_var)>1, length(time_var)>1)) {
    stop("The longitudinal variable should have length 1. Is it a wide data instead?")
  }
  long_dat<-long_dat[complete.cases(
    long_dat[,c(Ai_var,last_obs_time_var,
                cov_var,ID_var,death_ind_var,drop_ind_var)]),]
  
  
  long_dup_dat <- long_dat[,c(Ai_var,last_obs_time_var,
                              cov_var,ID_var,death_ind_var,drop_ind_var)]
  long_uni_dat <- unique(long_dup_dat)
  ID_uniq<-long_uni_dat[[ID_var]]
  if (any(duplicated(ID_uniq))) {
    dup_loc<-which(duplicated(ID_uniq))
    stop("Some base variables have different value at different time. Please double check. 
         It happens at observation ",dup_loc)
  }
  nSub <- nrow(long_uni_dat)
  
  ID_uniq<-long_uni_dat[[ID_var]]
  Ai <- long_uni_dat[[Ai_var]]
  last_obs_time <- long_uni_dat[[last_obs_time_var]]
  
  
  if (is.null(cov_var)) {
    Xi <- matrix(0, nrow = nSub, ncol = 1)
    Xi <- split(Xi, seq(nSub))
    sd_be_scale <- 1
    mean_be_scale <- 0
  } else {
    Xi <- long_uni_dat[, cov_var, drop = FALSE]
    sd_be_scale <- apply(Xi, 2, sd)
    mean_be_scale <- apply(Xi, 2, mean)
    Xi <- scale(Xi)
    Xi <- split(Xi, seq(nSub))
    
  }
  

  longi_outcome_time_long_dat <- long_dat[,c(ID_var,longi_outcome_var,time_var)]
  longi_outcome_time_long_dat_comp<-na.omit(longi_outcome_time_long_dat)
  longi_outcome_time_long_dat_cov_comp<-merge(
    longi_outcome_time_long_dat_comp,
    long_uni_dat,
    by=ID_var,
    all.x=TRUE)
  
  colnames(longi_outcome_time_long_dat_cov_comp)[
    which(colnames(longi_outcome_time_long_dat_cov_comp) == longi_outcome_var)] <- "outcome"
  colnames(longi_outcome_time_long_dat_cov_comp)[
    which(colnames(longi_outcome_time_long_dat_cov_comp) == time_var)] <- "time"
  longi_outcome_time_long_dat_cov_comp$longi_tstar<-
    longi_outcome_time_long_dat_cov_comp[[last_obs_time_var]]-
    longi_outcome_time_long_dat_cov_comp$time
  
  
  longi_outcome_list <- list()
  longi_time_list <- list()
  tstar_list <- list()
  ii <- 0
  ID_longi<-longi_outcome_time_long_dat_cov_comp[[ID_var]]
  
  longi_outcome_all<-longi_outcome_time_long_dat_cov_comp$outcome
  mean_outcome<-mean(longi_outcome_all)
  sd_outcome<-sd(longi_outcome_all)
  longi_outcome_all<-(longi_outcome_all-mean_outcome)/sd_outcome
  longi_outcome_time_long_dat_cov_comp$outcome<-longi_outcome_all
  
  longi_time_all<-longi_outcome_time_long_dat_cov_comp$time
  longi_tstar_all<-longi_outcome_time_long_dat_cov_comp$longi_tstar
  if (any(longi_tstar_all<0)) {
    warning("Some observed timepoints are later than death/censor time. They are set to be 
            death/censor time.")
    longi_time_all[longi_tstar_all<0]<-
      longi_outcome_time_long_dat_cov_comp[[last_obs_time_var]][longi_tstar_all<0]
    longi_tstar_all[longi_tstar_all<0]<-0
    longi_outcome_time_long_dat_cov_comp$longi_tstar<-longi_tstar_all
  }
  
  longi_outcome_time_long_dat_cov_comp$ID_new<-longi_outcome_time_long_dat_cov_comp[[ID_var]]
  for (i in ID_uniq) {
    ii <- ii + 1
    longi_outcome_list[[ii]] <- longi_outcome_all[ID_longi == i]
    longi_time_list[[ii]] <- longi_time_all[ID_longi == i]
    tstar_list[[ii]] <- longi_tstar_all[ID_longi == i]
  }
  
  death_ind <- long_uni_dat[[death_ind_var]]
  if (is.null(drop_ind_var)) {
    drop_ind <- rep(0, nSub)
  } else {
    drop_ind <- long_uni_dat[[drop_ind_var]]
  }
  
  sty_end_time <- max(last_obs_time)
  
  if (is.null(longi_end_time)) {
    longi_end_time <- sty_end_time
  }
  if (longi_end_time > sty_end_time) {
    longi_end_time <- sty_end_time
  }
  
  delta_value <- ifelse(death_ind == 1, 0,
                        ifelse(drop_ind == 1, 1, 2))
  
  
  ti1s <-
    last_obs_time[delta_value == 0 & sapply(longi_time_list, length) != 0]
  yi1s <-
    longi_outcome_list[delta_value == 0 &
                         sapply(longi_time_list, length) != 0]
  Ai1s <- Ai[delta_value == 0 & sapply(longi_time_list, length) != 0]
  Xi1s <- Xi[delta_value == 0 & sapply(longi_time_list, length) != 0]
  tstari1s <-
    tstar_list[delta_value == 0 & sapply(longi_time_list, length) != 0]
  
  ti2s <-
    last_obs_time[delta_value == 1 & sapply(longi_time_list, length) != 0]
  yi2s <-
    longi_outcome_list[delta_value == 1 &
                         sapply(longi_time_list, length) != 0]
  Ai2s <- Ai[delta_value == 1 & sapply(longi_time_list, length) != 0]
  Xi2s <- Xi[delta_value == 1 & sapply(longi_time_list, length) != 0]
  tstari2s <-
    tstar_list[delta_value == 1 & sapply(longi_time_list, length) != 0]
  longi.t2s <-
    longi_time_list[delta_value == 1 &
                      sapply(longi_time_list, length) != 0]
  
  ti3s <-
    last_obs_time[delta_value == 2 & sapply(longi_time_list, length) != 0]
  yi3s <-
    longi_outcome_list[delta_value == 2 &
                         sapply(longi_time_list, length) != 0]
  Ai3s <- Ai[delta_value == 2 & sapply(longi_time_list, length) != 0]
  Xi3s <- Xi[delta_value == 2 & sapply(longi_time_list, length) != 0]
  tstari3s <-
    tstar_list[delta_value == 2 & sapply(longi_time_list, length) != 0]
  longi.t3s <-
    longi_time_list[delta_value == 2 &
                      sapply(longi_time_list, length) != 0]
  
  ti4s <-
    last_obs_time[delta_value == 0 & sapply(longi_time_list, length) == 0]
  Ai4s <- Ai[delta_value == 0 & sapply(longi_time_list, length) == 0]
  Xi4s <- Xi[delta_value == 0 & sapply(longi_time_list, length) == 0]
  
  ti5s <-
    last_obs_time[delta_value == 1 & sapply(longi_time_list, length) == 0]
  Ai5s <- Ai[delta_value == 1 & sapply(longi_time_list, length) == 0]
  Xi5s <- Xi[delta_value == 1 & sapply(longi_time_list, length) == 0]
  
  ti6s <-
    last_obs_time[delta_value == 2 & sapply(longi_time_list, length) == 0]
  Ai6s <- Ai[delta_value == 2 & sapply(longi_time_list, length) == 0]
  Xi6s <- Xi[delta_value == 2 & sapply(longi_time_list, length) == 0]
  
  return(
    list(
      ti1s = ti1s,
      yi1s = yi1s,
      Ai1s = Ai1s,
      Xi1s = Xi1s,
      tstari1s = tstari1s,
      ti2s = ti2s,
      yi2s = yi2s,
      Ai2s = Ai2s,
      Xi2s = Xi2s,
      tstari2s = tstari2s,
      longi.t2s = longi.t2s,
      ti3s = ti3s,
      yi3s = yi3s,
      Ai3s = Ai3s,
      Xi3s = Xi3s,
      tstari3s = tstari3s,
      longi.t3s = longi.t3s,
      ti4s = ti4s,
      Ai4s = Ai4s,
      Xi4s = Xi4s,
      ti5s = ti5s,
      Ai5s = Ai5s,
      Xi5s = Xi5s,
      ti6s = ti6s,
      Ai6s = Ai6s,
      Xi6s = Xi6s,
      delta_value = delta_value,
      sd_cov_be_scale = sd_be_scale,
      mean_cov_be_scale = mean_be_scale,
      mean_outcome = mean_outcome,
      sd_outcome = sd_outcome,
      longi_dat_cov_comp = longi_outcome_time_long_dat_cov_comp,
      longi_end_time = longi_end_time,
      sty_end_time = sty_end_time
    )
  )
  
}
