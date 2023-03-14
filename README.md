# RetroJM

## Overview
This package provides a way to perform retrospective joint modeling method.


## Installation
```r
## From GitHub: 
devtools::install_github("quranwu/RetroJM")
```

## Usage
Detailed instructions can be found in the vignette file.

## Example
The example data is available in ./data folder.

```r
library(RetroJM)

## load the simulated example data Retro_data
data(Retro_data)

test_run<-RetroJM(data = Retro_data,
  long_format = FALSE,
  outcome_var = paste0("Y",1:13),
  time_var = paste0("time",1:13),
  Ai_var = "Ai",
  last_obs_time_var = "observed.t",
  cov_var = paste0("x",1:2),
  ID_var = "ID",
  surv_censor_var = "surv_ind",
  drop_censor_var = "drop_ind",
  longi_end_time = 25,
  surv_drop_knots = 1,
  longi_inter_knots = 1,
  num_cores = 2)
```
Once the analysis is done, `test_run` is a list with four elements, each for survival related parameters, drop out related parameters, longitudinal outcome related parameters, and list of raw results.


## Run simulated dataset in the paper: 

The dataset can be prepared as below. The raw dataset can be found in "./data/simulate_data". To prepare `data1.Rdata`, for example, 

``` r

data1$surv_ind<-ifelse(unlist(data1$delta)==1,1,0)
data1$drop_ind<-ifelse(unlist(data1$delta)==2,1,0)

yi<-data1$Yi
yi_dat<-do.call(rbind,yi)
colnames(yi_dat)<-paste0("Y",1:ncol(yi_dat))

longi_time<-data1$long.t
longi_time_full<-lapply(longi_time,function(x) {
  save_vec<-rep(NA,13)
  if (length(x)!=0) {
    save_vec[1:length(x)]<-x
  }
  return(save_vec)
})

Xi<-data1$Xi
Xi_dat<-do.call(rbind,Xi)
colnames(Xi_dat)<-paste0("x",1:2)

time_dat<-do.call(rbind,longi_time_full)

colnames(time_dat)<-paste0("time",1:ncol(time_dat))

data_use<-cbind(data1[,c("Ai","observed.t","surv_ind","drop_ind")],Xi_dat,yi_dat,time_dat)
data_use$ID<-rownames(data_use)
```

`data_use` is the prepared dataset. It can be fit into the model as: 

```r
test_run_0<-RetroJM(data = data_use,
                  long_format = FALSE,
                  outcome_var = paste0("Y",1:ncol(yi_dat)),
                  time_var = paste0("time",1:ncol(time_dat)),
                  Ai_var = "Ai",
                  last_obs_time_var = "observed.t",
                  cov_var = paste0("x",1:2),
                  ID_var = "ID",
                  surv_censor_var = "surv_ind",
                  drop_censor_var = "drop_ind",
                  longi_end_time = 25,
                  surv_drop_knots = 1,
                  longi_inter_knots = 0,
                  num_cores = 5)
```

The `longi_inter_knots` can be modified with different number of longitudinal knots, and then the model selection can be conducted. 

