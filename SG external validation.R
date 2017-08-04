# SG external validation
rm(list=ls())

# External validation of CU and UKPDS mortality model for diabetic patients
# Using HKU data
library(rms); library(ggplot2); library(survival); library(dplyr)

load("vfm/mortality_t2_5.Rdata")
load("mortality survey all_predictors.Rdata")
test_model <- approx$approx_model

# Data---------------
d$event <- as.numeric(d$event)
d$years <- d$time / 365.25
units(d$years) <- "year"
label(d$years) <- "Survival Time"

completed <- list()
 for (i in 1:imp$n.impute){
 # Extract the 1st imputed dataset
completed[[i]] <- d
imputed <- impute.transcan(imp, imputation=i, data=d, list.out=TRUE, pr=FALSE, check=FALSE) 
completed[[i]][names(imputed)] <- imputed
}
names(completed) <- paste0("i", 1:imp$n.impute)

data_manipulation <- function(dataset){
# Convert TRUE:FALSE to 1:0
      convert_TF <- which(sapply(dataset, function(x) levels(x)) %in% list(c('FALSE', 'TRUE')))
      for(i in convert_TF){
      	dataset[,i] <- as.numeric(as.logical(paste0(dataset[,i])))
      }
      
# Construct variables
      dataset <- dataset %>%
            mutate(insulin = as.numeric(meds == 'insulin'),
                   egfr = exp(log.egfr_chinese),
                   acr = exp(log.urine_acr),
                   non_hdl = tc - hdl,
                   smoke_binary = as.numeric(smoking == 'Yes'))
      return(dataset)
}
completed_ukcu <- lapply(completed, data_manipulation)

# Calculate the C-index
calc_c_index <- function (dataset){
pred <- predict(test_model, newdata=dataset)
score_result <- data.frame(pred, years=dataset$years,
                                       event=dataset$event)
c(survConcordance(Surv(years, event) ~ pred, score_result)$concordance)
}
c_index <- do.call(rbind, lapply(completed, calc_c_index))
c("mean"=mean(c_index), "min"=min(c_index), "max"=max(c_index))     

# Calibration using one (imputed) dataset
# cut data.frame into 10 group
library(pec)
library(reshape2)  

FUN.calib_mean_pred <- function(dataset) {
 dataset$pred <- 1 - predictSurvProb(test_model, dataset, times=5)
  # if we want to predict 10-year risk, change times=5 to times=10.
  # event after 5 years is set to censored status? If so, run following code 
  # completed$event[completed$event == 1 & completed$years > 5] <- 0
  # table(completed$event)
  # if above annotation code is executed, c-statistics for each group is slightly larger.
  # predict 10 group's risk probabilities
  
  # predict survival probability and observed probability for overall
  dataset <- dataset[order(dataset$pred),]
  dataset <- dataset[, c("event", "years", "pred"),]
  dataset$group <- cut2(dataset$pred, g = 10)
  dataset$event <- as.numeric(dataset$event)
  levels(dataset$group) <- c(1:10)
  observed <- rep(NA, 10)
  for (i in 1:10) {
    km <- survfit(Surv(years, event) ~ 1, data = dataset[dataset$group == i, ])
    survest <- stepfun(km$time, c(1, km$surv))
    observed[i] <- 1-survest(5) # risk at 5-years
  }
  mean_pred <- tapply(dataset$pred, dataset$group, mean)
  dat <- data.frame(observed, pred = mean_pred, Deciles=1:10)
  risk_group <- dataset$group
  melt(dat, id.vars = 'Deciles', variable.name = 'Group', value.name = 'Risk')
}

# calib_data <- lapply(completed, FUN.calib_mean_pred)
calib_data <- FUN.calib_mean_pred(completed[[1]])
colours=c("skyblue1", "navyblue")
plot_all_bars(calib_data)
