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

# UKPDS
FUN.stroke_ukpds <- function(period){
d <- period
d$smoking <- ifelse(d$smoking=="Yes", 1, 0) # No option for ex-smokers

q <-0.00186* 1.092^(d$age_entry-55)* 0.7^(d$female)* 1.547^(d$smoking)*
	8.554^(d$af)* 1.122^((d$sbp-135.5)/10)* 1.138^(d$lr-5.11)
(1-exp(-q* 1.145^(d$duration)* (1-1.145^(5))/(1-1.145)))*100
}


# stroke
CU_score <- function(testdata){
  with(testdata, 0.0634*age.entry + 0.0897*hba1c + 0.5314*log(acr, 10) + 0.5636*chd)
}



FUN.chd_ukpds <- function(period){
d <- period
d$smoking <- ifelse(d$smoking=="Yes", 1, 0) # No option for ex-smokers

q <-0.0112* 1.059^(d$age_entry-55)* 0.525^(d$female)* 1.35^(d$smoking)* 
	1.183^(d$hba1c-6.72)* 1.088^((d$sbp-135.7)/10)* 3.845^(log(d$lr)-1.59)
(1-exp(-q* 1.078^(d$duration)* (1-1.078^5)/(1-1.078)))*100
}


# for non_hdl, use tc - hdl (total cholesterol - hdl cholesterol)
# chd
FUN.chd_cuhk<-function (d, egfr, urine_acr, non_hdl){
names(d)[names(d)%in%egfr]<-"egfr"
names(d)[names(d)%in%urine_acr]<-"urine_acr"
names(d)[names(d)%in%non_hdl]<-"non_hdl"

q<-with(d, 0.0267*(age_entry) - 0.3536*(female) + 0.4373*(smoking) + 0.0403*(duration)
	-0.4808*log10(egfr) + 0.1232 *log10(1+(urine_acr)) + 0.2644*(non_hdl))
(1-0.9616^exp(0.9440*(q-0.7082)))*100
}


