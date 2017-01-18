rm(list=ls())
library(rms); library(Hmisc)

setwd("C:\\Users\\tomli\\Desktop\\Tom folder\\Chao\\Chao - 29Sep2016\\")

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
# vars <- c("rcs(hba1c, nk) + rcs(age,4) + af + cancer + pad + stroke + chd + ckd + rcs(map, nk) + rcs(sbp, nk) + rcs(dbp, nk) + rcs(ldl, nk) + female")
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

# Save the model
fullmodel$x <- NULL # remove raw data
fullmodel$y <- NULL # remove raw data
save(fullmodel, file="External validation\\SG_full_model.Rdata")
