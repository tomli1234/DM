#----------------------------------------------
# Generate descriptive table 
# for risk factors stratified by 
# diagnosis cohort and observational period
# Tom, March 2017
#-----------------------------------------------

#########

rm(list=ls())

library(tidyverse)
library(data.table)

setwd('C:\\Users\\tomli\\Desktop\\Tom folder\\Chao\\Data\\risk factor')

# Read data---------------
load('haemoglobin_clean.RData')
# load('urine_acr_clean.RData')
# load('hba1c_clean.RData')
# load('cancer.RData')
# load('af.RData')
# load('patient.RData')
# load('wbc_clean.RData')

haemoglobin <- haemoglobin %>%
      mutate(Year = format(ref_date,"%Y"),
             Month = format(ref_date,"%m"),
             Obser_period = cut(as.numeric(Year), c(2006, 2009, 2012, 2015), right = FALSE),
             Obser_period = forcats::fct_recode(Obser_period, 
                                 "1" = "[2006,2009)",
                                 "2" = "[2009,2012)",
                                 "3" = "[2012,2015)"),
             Obser_period = as.numeric(paste0(Obser_period)),
             ref_date = as.numeric(ref_date),
             serial_no = as.numeric(serial_no)) %>%
      select(-Year, -Month) %>%
      data.table()

# haemoglobin <- haemoglobin %>%
# 
#       # Split by observation period
#       split(.$Obser_period) %>%
#       lapply(function(x) {
# 
#             # Split by serial_no      
#             x %>% split(.$serial_no) %>%
#                   lapply(function(y) {
# 
#                         # select the latest record
#                         y <- y %>% filter(ref_date == max(ref_date)) 
# 
#                         # Return the average value
#                         mean(y$haemoglobin)
#                   })
#       }) 

# The R code above consumes so much memory
# So try using rcpp
library(Rcpp)

# This c++ function calculate the mean of 
# the latest records of the patients in a given period
cppFunction('double mean_latest(NumericMatrix x) {
      int nrow = x.nrow();
      int index_no = 0;
      int n_patient = 0;
      int ref_date = 0;
      double cumcum = 0;
      double latest_i = 0;
      for(int i=0; i<nrow; i++) {
            if(x(i, 0) == index_no) {
                  if(x(i, 1) >= ref_date) {
                        latest_i = x(i, 2);
                  }
            }
            if(x(i, 0) != index_no) {
                  n_patient += 1;
                  cumcum += latest_i;
                  ref_date = 0;
                  if(x(i, 1) >= ref_date) {
                        latest_i = x(i, 2);
                  }
            }
            index_no = x(i, 0);
      }
      cumcum += latest_i;
      return cumcum/n_patient;
}')
mean_latest(as.matrix(haemoglobin[Obser_period == 1,]))
mean_latest(as.matrix(haemoglobin[Obser_period == 2,]))
mean_latest(as.matrix(haemoglobin[Obser_period == 3,]))








