# External validation of HKU mortality model for diabetic patients
# Using Singapore's data
# Tom Li, tomli123@hku.hk, 25 Oct 2016
library(rms); library(Hmisc); library(ggplot2); library(survival)
rm(list=ls())

setwd("C:\\Users\\tomli\\Desktop\\Tom folder\\Chao\\Chao - 29Sep2016\\External validation\\CU")

load("mortality_full_8.Rdata") 
load("SG_full_model.Rdata") # HKU model


imp$formula; imp$call
d$event <- as.numeric(d$event)
d$years <- d$time / 365.25
units(d$years) <- "year"
label(d$years) <- "Survival Time"

# Extract the 1st imputed dataset
completed <- d
imputed <- impute.transcan(imp, imputation=1, data=d, list.out=TRUE, pr=FALSE, check=FALSE) 
completed[names(imputed)] <- imputed

# Convert 1:0 to TRUE:FALSE
# convert_TF <- which(sapply(completed,function(x) levels(x)) %in% list(paste0(0:1)))
# for(i in convert_TF){
# 	completed[,i] <- completed[,i]==1
# }

# CU model risk score
CU_score <- function(testdata){	 
      with(testdata, 0.0586*age - 0.7049*as.numeric(female) + 0.1078*abs(bmi-26)+ 
                 0.1469*log(acr+1, 10)^2 - 0.0165*as.numeric(ifelse(exp(log.egfr)<60, exp(log.egfr), 60)) +
                 -0.1913*as.numeric(haemoglobin) + 0.5698*as.numeric(pad) + 1.3384*as.numeric(cancer) + 0.7035*insulin)
}

CU_score(completed)


# Discrimination-------------------------------------------------------
completed$pred <- CU_score(completed)
score_result <- data.frame(pred=completed$pred, years=completed$years, event=completed$event)
c_index <- c(survConcordance(Surv(years, event) ~ pred, score_result)$concordance)
#-------------------------------------------------------------------------


# Calibration----------------------------------------------------------
# predicted survival at 5-years
completed$pred <- 1-0.9567^exp(0.9768*(CU_score(completed) + 0.0415))

FUN.deciles_mean_pred <- function(data) {
  data <- data[order(data$pred),]
  data$group <- cut2(data$pred, g=10)
  data$years <- data$time / 365.25
  data$event <- as.numeric(data$event)
  levels(data$group) <- c(1:10)
  observed <- rep(NA, 10)
  for (i in 1:10) {
    km <- survfit(Surv(years, event)~1, data=data[data$group==i, ])
    survest <- stepfun(km$time, c(1, km$surv))
    observed[i] <- 1-survest(5) # risk at 5-years
  }
  mean_pred <- tapply(data$pred, data$group, mean)
  data.frame(observed, pred = mean_pred, Deciles=1:10)
}

dat <-FUN.deciles_mean_pred(completed)
library(reshape2)
dat.melt <- melt(dat, id.vars = 'Deciles', variable.name = 'Group', value.name = 'Risk')

p <- ggplot(dat.melt, aes(x = Deciles, y = Risk*100, color = Group)) + 
		geom_point(size = 4) +  xlab("Tenth of predicted risk") + 
		ylab("5 year risk (%)") + 
		theme(axis.title.x = element_text(size=14), 
			axis.title.y = element_text(size=14), 
			axis.text.y = element_text(size=12)) + 
		scale_x_continuous(breaks = NULL) + 
		scale_y_continuous(breaks=seq(0, 100, 10)) + 
		scale_color_manual(values = c("observed" = 'skyblue1','pred' = 'darkblue'), 
			labels=c("observed" = 'Observed','pred' = 'Predicted')) + 
		theme(legend.position="bottom", 
			legend.title=element_blank(), 
			legend.text = element_text(size = 12, face = "bold")) +
		annotate("text",3,100*max(dat.melt$Risk),size=8,
			label=paste0("C-statistic: ",round(c_index,3)))
#-------------------------------------------------------------------------

# Results output----------------------------------------------------------
ggsave(p, file="validation_results.pdf")
#-------------------------------------------------------------------------
