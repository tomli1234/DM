rm(list=ls())

# External validation of CU and UKPDS mortality model for diabetic patients
# Using HKU data
# Tom Li, tomli123@hku.hk, 21 Mar 2017
library(rms); library(Hmisc); library(ggplot2); library(survival); library(dplyr)

# Data--------------------------------------------------------------------------------------------
data_manipulation <- function(d, imp) {
      d$event <- as.numeric(d$event)
      d$years <- d$time / 365.25
      units(d$years) <- "year"
      label(d$years) <- "Survival Time"
      
      # Extract the 1st imputed dataset
      completed <- d
      imputed <- impute.transcan(imp, imputation=1, data=d, list.out=TRUE, pr=FALSE, check=FALSE) 
      completed[names(imputed)] <- imputed
      
      # Convert TRUE:FALSE to 1:0
      convert_TF <- which(sapply(completed,function(x) levels(x)) %in% list(c('FALSE', 'TRUE')))
      for(i in convert_TF){
      	completed[,i] <- as.numeric(as.logical(paste0(completed[,i])))
      }
      
      # Construct variables
      completed <- completed %>%
            mutate(insulin = as.numeric(meds == 'insulin'),
                   egfr = exp(log.egfr_chinese),
                   acr = exp(log.urine_acr),
                   non_hdl = tc - hdl,
                   smoke_binary = as.numeric(smoking == 'Yes'))
      return(completed)
}


# CU and UKPDS models------------------------------------------------------------------------
# CU mortality
CU_mortality <- function(testdata){	 
      pred_score <- with(testdata, 0.0586*age - 0.7049*female + 0.1078*abs(bmi-26)+ 
                 0.1469*log(acr+1, 10)^2 - 0.0165*as.numeric(ifelse(egfr<60, egfr, 60)) +
                 -0.1913*haemoglobin + 0.5698*pad + 1.3384*cancer + 
                 0.7035*insulin)
      
      pred_surv <- 1-0.9567^exp(0.9768*(pred_score + 0.0415))
      
      return(list(pred_score = pred_score, pred_surv = pred_surv))
}

# CU stroke
CU_stroke <- function(testdata){
      pred_score <- with(testdata, 0.0634*age + 0.0897*hba1c + 0.5314*log(acr, 10) + 0.5636*chd)
      pred_surv <- 1-0.9707^exp(pred_score - 4.5674)

      return(list(pred_score = pred_score, pred_surv = pred_surv))
}

# CU CHD
# for non_hdl, use tc - hdl (total cholesterol - hdl cholesterol)
CU_chd <- function (d){
      pred_score <- with(d, 0.0267*age - 0.3536*female + 0.4373*smoke_binary + 0.0403*duration
              -0.4808*log10(egfr) + 0.1232 *log10(1+(acr)) + 0.2644*(non_hdl))
      
      pred_surv <- 1-0.9616^exp(0.9440*(pred_score-0.7082))
      
      return(list(pred_score = pred_score, pred_surv = pred_surv))
}

CU_models <- list(mortality=CU_mortality, stroke=CU_stroke, chd=CU_chd)

# UKPDS stroke
UKPDS_stroke <- function(d){
      pred_score <- with(d, 0.00186* 1.092^(age-55)* 0.7^(female)* 1.547^(smoke_binary)*
            8.554^(af)* 1.122^((sbp-135.5)/10)* 1.138^(lr-5.11))
      pred_surv <- 1-exp(-pred_score* 1.145^(d$duration)* (1-1.145^(5))/(1-1.145))
      
      return(list(pred_score = pred_score, pred_surv = pred_surv))
      
}

# UKPDS CHD
UKPDS_chd <- function(d){
      pred_score <- 0.0112* 1.059^(d$age-55)* 0.525^(d$female)* 1.35^(d$smoke_binary)* 
            1.183^(d$hba1c-6.72)* 1.088^((d$sbp-135.7)/10)* 3.845^(log(d$lr)-1.59)
      pred_surv <- 1-exp(-pred_score* 1.078^(d$duration)* (1-1.078^5)/(1-1.078))
      return(list(pred_score = pred_score, pred_surv = pred_surv))
      
}

UKPDS_models <- list(mortality=NULL, stroke=UKPDS_stroke, chd=UKPDS_chd)


# C-index and calibration---------------------------------------------------------------------
FUN.deciles_mean_pred <- function(data) {
      data <- data %>%
            arrange(pred_surv) %>%
            mutate(group = cut2(pred_surv, g=10),
                   years = time / 365.25,
                   event = as.numeric(event))
      levels(data$group) <- c(1:10)
      
      observed <- rep(NA, 10)
      for (i in 1:10) {
            km <- survfit(Surv(years, event)~1, data=data[data$group==i, ])
            survest <- stepfun(km$time, c(1, km$surv))
            observed[i] <- 1-survest(5) # risk at 5-years
      }
      mean_pred <- tapply(data$pred_surv, data$group, mean)
      data.frame(observed, pred = mean_pred, Deciles=1:10)
}

plot_bar <- function(dat.melt, colours = c('skyblue1','darkblue')){
      ggplot(dat.melt, aes(x = factor(Deciles), y = Risk*100, fill = Group)) + 
            geom_bar(stat = 'identity', position = 'dodge', width = 0.7) +  
            xlab("Tenth of predicted risk") + 
            ylab("5 year risk (%)") + 
            theme_classic()+
            theme(axis.title.x = element_text(size=14), 
                  axis.title.y = element_text(size=14), 
                  axis.text.y = element_text(size=12)) + 
            scale_x_discrete() + 
            scale_y_continuous(breaks=seq(0, 100, 10)) + 
            scale_fill_manual(values = c("observed" = colours[1],'pred' = colours[2]),
                              labels=c("observed" = 'Observed','pred' = 'Predicted')) +
            theme(legend.position="bottom", 
                  legend.title=element_blank(), 
                  legend.text = element_text(size = 12, face = "bold"),
                  panel.grid.major.y = element_line( size=.1, color="grey"))
      
}

plot_cal_bar <- function(dat.melt, kolor = colours){
  ggplot(dat.melt, aes(x = Deciles, y = Risk*100, fill = factor(Group, levels = c("pred", "observed")))) + 
    geom_bar(stat = 'identity', position = 'dodge', width = 0.7, color = "black") +  
    xlab("Tenth of predicted risk") + 
    ylab("5 year risk (%)") + 
    theme(axis.title.x = element_text(size=14), 
          axis.title.y = element_text(size=14), 
          axis.text.y = element_text(size=12)) + 
    scale_x_discrete() + 
    scale_y_continuous(breaks=seq(0, 100, 10)) + 
    scale_fill_manual(values = c('pred' = kolor[1], "observed" = kolor[2]),
                      labels=c("observed" = 'Observed','pred' = 'Predicted')) +
    theme_minimal() + 
    theme(legend.position="bottom", 
          legend.title=element_blank(), 
          legend.text = element_text(size = 12, face = "bold"),
          panel.grid.major.y = element_line( size=.1, color="grey"))
}

external_validation <- function(outcome, CU_UKPDS, sex) {
      
      # Select CU/UKPDS model
      if(CU_UKPDS == 'JADE') {
            model <- CU_models
      } else if(CU_UKPDS == 'UKPDS'){
            model <- UKPDS_models
      }
      
      # Select outcome
      if(outcome == 'mortality') {
            model_score <- model$mortality
            load(wd_mortality) 
      } else if(outcome == 'stroke') {
            model_score <- model$stroke
            load(wd_stroke) 
      } else if(outcome == 'chd') {
            model_score <- model$chd
            load(wd_chd) 
      }
      
      completed <- data_manipulation(d, imp)
      completed$pred_score <- model_score(completed)$pred_score
      completed$pred_surv <- model_score(completed)$pred_surv
      
      # Select sex
      if(sex == 'Female') {
            completed <- completed[completed$female == 1, ]
      } else if(sex == 'Male') {
            completed <- completed[completed$female == 0, ]
      } else if(sex == 'All') {
            
      }    
      
      # Discrimination-------------------------------------------------------
      score_result <- data.frame(pred_score=completed$pred_score, years=completed$years, event=completed$event)
      c_index <- c(survConcordance(Surv(years, event) ~ pred_score, score_result)$concordance)
      
      # Calibration----------------------------------------------------------
      # predicted survival at 5-years
      dat <-FUN.deciles_mean_pred(completed)
      library(reshape2)
      dat.melt <- melt(dat, id.vars = 'Deciles', variable.name = 'Group', value.name = 'Risk')
      p <- plot_cal_bar(dat.melt) + ggtitle(paste(CU_UKPDS, "model"))
      return(list(c_index = c_index, figure = p))
}

path <- "C:\\Users\\tomli\\Desktop\\Tom folder\\Chao\\Chao - 29Sep2016\\External validation\\CU"
setwd("E:/vfm")
#wd_mortality <- "mortality_20170320_5.Rdata" # Mortality data
wd_mortality <- "mortality_no-insulin_3" # Mortality data
#wd_mortality <- "mortality_t2_5.Rdata" # Mortality data
#wd_stroke <- "stroke_20170320_5.Rdata" # Stroke data
wd_stroke <- "stroke_no-insulin_3" # Stroke data
#wd_stroke <- "stroke_t2_3.Rdata" # Stroke data
#wd_chd <- "chd_20170320_5.Rdata" # CHD data
wd_chd <- "chd_no-insulin_3.Rdata" # CHD data
#wd_chd <- "chd_t2_5.Rdata" # CHD data

# outcome: 'chd', 'mortality', 'stroke'
# CU_UKPDS: 'CU', 'UKPDS'
# sex: 'Female', 'Male', 'All' 

colours <- c('skyblue1','darkblue')
m_jade <- external_validation(outcome = "mortality", CU_UKPDS = "JADE", sex = "All")
colours <- c('lightgreen','darkgreen')
stroke_jade <- external_validation(outcome = "stroke", CU_UKPDS = "JADE", sex = "All")
stroke_uk <- external_validation(outcome = "stroke", CU_UKPDS = "UKPDS", sex = "All")
colours <- c('pink','darkred')
chd_jade <- external_validation(outcome = "chd", CU_UKPDS = "JADE", sex = "All")
chd_uk <- external_validation(outcome = "chd", CU_UKPDS = "UKPDS", sex = "All")

m_jade$c_index
stroke_jade$c_index
stroke_uk$c_index
chd_uk$c_index
chd_jade$c_index

# Table of c index
c_index_all <- expand.grid(outcome = c('mortality','stroke','chd'), 
                           CU_UKPDS = c('JADE','UKPDS'), 
                           sex = c('Male','Female','All')) %>%
  
  filter(! (outcome == 'mortality' & CU_UKPDS == 'UKPDS')) %>%
  
  mutate(c_index = apply(., 1, function(x)
    external_validation(x[1],x[2],x[3])$c_index))

c_index_all %>%
  ggplot(aes(x = outcome, y = c_index, fill = CU_UKPDS)) +
  geom_bar(stat= 'identity', position = 'dodge', width = 0.5) +
  facet_grid(. ~ sex)

# multiple plot function
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL){
  library(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

load("E:/Mortality/all dm egfr sbp lr.Rdata")
colours <- c('skyblue1','darkblue')
mortality <- plot_cal_bar(boots_results$calibration) + ggtitle("HKU model")
load("E:/stroke/all dm creatinine map lr.Rdata")
colours <- c('lightgreen','darkgreen')
stroke <- plot_cal_bar(boots_results$calibration) + ggtitle("HKU model")
load("E:/CHD/all dm creatinine dbp lr.Rdata")
colours <- c('pink','darkred')
chd <- plot_cal_bar(boots_results$calibration) + ggtitle("HKU model")

multiplot(mortality, m_jade$figure)
multiplot(stroke, stroke_jade$figure, stroke_uk$figure)
multiplot(chd, chd_jade$figure, chd_uk$figure)

