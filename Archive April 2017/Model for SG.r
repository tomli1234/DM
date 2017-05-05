rm(list=ls())
library(rms); library(Hmisc)

setwd("")

load("mortality_full_8.Rdata")
# load("mortality_20160926_incident.Rdata")
# load("mortality_20160926_cancer.Rdata")
# load("mortality_20160926_nonincident.Rdata")

imp$formula; imp$call
d$event <- as.numeric(d$event)
d$years <- d$time / 365.25
units(d$years) <- "year"
label(d$years) <- "Survival Time"

# Dean's suggestion-----------------------------------
# Full model (with age:duration interaction)
nk <- 4
vars <- c("rcs(log.urine_acr, nk) + rcs(bmi, nk) + rcs(haemoglobin, nk) + rcs(wbc, nk) + rcs(pulse, nk) + rcs(sbp, nk) + rcs(lr, nk) + rcs(dbp, nk) + rcs(hba1c, nk) + rcs(log.creatinine, nk) + rcs(age,4)*rcs(duration,4) + smoking + af + cancer + pad + stroke + chd + ckd + meds + complications + rcs(map, nk) + rcs(ldl, nk) + rcs(hdl, nk) + rcs(tc, nk) + rcs(triglyceride, nk) + female")
# vars <- c("rcs(log.urine_acr, nk) + rcs(bmi, nk) + rcs(haemoglobin, nk) + rcs(wbc, nk) + rcs(pulse, nk) + rcs(sbp, nk) + rcs(lr, nk) + rcs(dbp, nk) + rcs(hba1c, nk) + rcs(log.creatinine, nk) + rcs(age, 4) + smoking + af + cancer + pad + stroke + chd + ckd + meds + complications + rcs(map, nk) + rcs(ldl, nk) + rcs(hdl, nk) + rcs(tc, nk) + rcs(triglyceride, nk) + female")
vars <- gsub("nk", nk, vars)
fm <- as.formula(paste("Surv(years, event) ~", vars))
# select variables
# remove sbp/dbp, use map
fm <- update(fm, ~. -rcs(sbp, 4) -rcs(dbp, 4))
# remove hdl/ldl/tc, use lr
fm <- update(fm, ~. -rcs(hdl, 4) -rcs(ldl, 4) -rcs(tc, 4))
# remove wbc
fm <- update(fm, ~. -rcs(wbc, 4))
# remove egfr, urine_acr, use creatinine
fm <- update(fm, ~. -rcs(egfr, 4) -rcs(log.urine_acr, 4))
# remove pre-existing medical conditions
fm <- update(fm, ~. - cancer - stroke - chd - meds - ckd - pad)
#-----------------------------------------------------

# SG available variables-----------------------------
nk <- 4
vars <- c("rcs(bmi, nk) + rcs(hba1c, nk) + rcs(age,4) + rcs(duration,4) + smoking + af + cancer + pad + stroke + chd + ckd + rcs(map, nk) + rcs(ldl, nk) + female")
# vars <- c("rcs(hba1c, nk) + rcs(age, 4) + af + cancer + pad + stroke + chd + ckd + rcs(map, nk) + rcs(sbp, nk) + rcs(dbp, nk) + rcs(ldl, nk) + female")
vars <- gsub("nk", nk, vars)
fm <- as.formula(paste("Surv(years, event) ~", vars))
#-----------------------------------------------------


# Stefan / SG available variables-----------------------------
nk <- 4
vars <- c("rcs(bmi, nk) + rcs(hba1c, nk) + rcs(age,4) + rcs(duration,4) + smoking + stroke + chd + rcs(map, nk) + rcs(lr, nk) + rcs(log.creatinine, nk) + female")
vars <- gsub("nk", nk, vars)
fm <- as.formula(paste("Surv(years, event) ~", vars))
#-----------------------------------------------------
completed <- d
imputed <- impute.transcan(imp, imputation=1, data=d, list.out=TRUE, pr=FALSE, check=FALSE) 
completed[names(imputed)] <- imputed
dd <<- datadist(completed); options(datadist ="dd")
fullmodel <- fit.mult.impute(fm, fitter = cph, xtrans = imp, data = d, x=TRUE, y=TRUE, surv=TRUE, time.inc=5) 
ggplot(Predict(fullmodel), sepdiscrete ='vertical', nlevels=4, vnames ='names')
plot(anova(fullmodel))
 
# Save the model---------------------------------------------
fullmodel$x <- NULL # remove raw data
fullmodel$y <- NULL # remove raw data
save(fullmodel, file="SG_full_model.Rdata")


# validate function (bootstrap) ---------------------------------------------
FunDiscrimination <- function (model, bootstrap) {
  set.seed (10)
  v <- validate (model, B=bootstrap)
  print(v)
  Dxy = v[rownames(v)=="Dxy", colnames(v)=="index.corrected"]
  c.stat <- Dxy/2+0.5
  print(paste("C-statistic =", round(c.stat, 4)))
  c.stat
}

c_index <- FunDiscrimination(fullmodel, 10)

# predicted scores
completed$pred <- (predict(fullmodel, newdata=completed))
# predicted survival at 5-years
library(pec)
completed$pred <- 1 - predictSurvProb(fullmodel, completed, times=5)

# Calibration - seperate for male and female?
# cal <- calibrate(fullmodel, u=5, B=10)
# plot(cal)
# cal.decile <- calibrate(model, u=5, B=10, cmethod="KM", m=(nrow(completed)/10))
# plot(cal.decile)

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