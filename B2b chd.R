#----------------------------------------------
# Model for chd
#-----------------------------------------------

library(rms)
library(Hmisc)
library(pec)
library(prodlim)
library(reshape2)

load("vfm/chd_20170320_5.Rdata")
imp$formula; imp$call
d$event <- as.numeric(d$event)
d$years <- d$time / 365.25
units(d$years) <- "year"
label(d$years) <- "Survival Time"

# imputing missing data
completed <- d
imputed <- impute.transcan(imp, imputation=1, data=d, list.out=TRUE, pr=FALSE, check=FALSE) 
completed[names(imputed)] <- imputed
# Since d$event and completed$event have two status, i.e., 
# 1 (censored indication) and 2 (event indication), 
# transform it to 0 and 1 status as usual.
completed$event <- completed$event-1
d$event <- d$event-1
dd <<- datadist(completed); options(datadist ="dd")

# (1) full model predictors for chd, no interaction terms 
vars <- c("rcs(age, nk) + rcs(duration, nk) + rcs(log.urine_acr, nk) + rcs(bmi, nk) + rcs(haemoglobin, nk) + rcs(wbc, nk) + rcs(pulse, nk) + rcs(lr, nk) + rcs(hba1c, nk) + rcs(log.creatinine, nk) + rcs(map, nk) + rcs(triglyceride, nk) + smoking + af + cancer + pad + ckd + meds + complications + female + stroke + rcs(sbp, nk) + rcs(dbp, nk) + rcs(ldl, nk) + rcs(hdl, nk) + rcs(tc, nk) + rcs(log.egfr_chinese, nk)")
nk <- 4
vars <- gsub("nk", nk, vars)
fm <- as.formula(paste("Surv(years, event) ~", vars))
fullmodel <- fit.mult.impute(fm, fitter = cph, xtrans = imp, data = d, x=TRUE, y=TRUE, surv=TRUE, time.inc=5)

# estimation curve for variable
ggplot(Predict(fullmodel), sepdiscrete ='vertical', nlevels=4, vnames ='names')
# estimation result
fullmodel

# figure 1 - importance of each variable
plot(anova(fullmodel, tol=1e-13))

# Model approximation
FUN.approx <- function(mod, completedata) {
  lp <- predict(mod)
  a <- ols(as.formula(paste0("lp ~ ", paste0(vars, collapse=""))), sigma=1, data=completedata)
  s <- fastbw(a, aic=10000000000)
  betas <- s$Coefficients
  X <- model.matrix(fm, data=completedata)
  ap <- X %*% t(betas)
  m <- ncol(ap) - 1
  fullchisq <- mod$stats[3]
  r2 <- frac <- numeric(m)
  for(i in 1:m) {
    lpa <- ap[,i]
    r2[i] <- cor(lpa, lp)^2
    fapprox <- cph(Surv(years, event) ~ lpa, data=completedata, x=TRUE, y=TRUE)
    frac[i] <- fapprox$stats[3]/fullchisq
  }
  plot(r2, frac, type="b")
  list(r2, frac)
}

# figure 2 - r^2 approximation
r <- FUN.approx(fullmodel, completed)
tab <- cbind(sort(r[[2]]), 1:length(r[[2]]))
tab
write.csv(tab, "table_r2.csv")

# 5-year risk prediction
# if we want to predict 10-year risk, change times=5 to times=10.
completed$pred <- 1 - predictSurvProb(fullmodel, completed, times=5)
table(completed$event)
# event after 5 years is setted to censored status? If so, run following code 
# completed$event[completed$event == 1 & completed$years > 5] <- 0
# table(completed$event)
# if above annotation code is executed, c-statistics for each group is slightly larger.

# C-statistic for overall
cstat <- function(completed){
  score_result <- cbind(completed$pred, completed$years, completed$event, completed$pred)
  cstats <- survConcordance(Surv(score_result[,2], score_result[,3])~ score_result[,1])$concordance
  cstats
}
# C-statistic for female group
cstat_female <- function(completed){
  score_result <- cbind(completed[completed$female==T, "pred"], completed[completed$female==T, "years"], completed[completed$female==T, "event"],completed[completed$female==T, "pred"])
  cstat_f <- survConcordance(Surv(score_result[,2], score_result[,3])~ score_result[,1])$concordance
  cstat_f
}
# C-statistics for male group
cstat_male <- function(completed){
  score_result <- cbind(completed[completed$female==F, "pred"],completed[completed$female==F, "years"], completed[completed$female==F, "event"],completed[completed$female==F, "pred"])
  cstat_m <- survConcordance(Surv(score_result[,2], score_result[,3])~score_result[,1])$concordance
  cstat_m
}
print(cat("c-stat=", cstat(completed), "\nfemale=", cstat_female(completed), "\nmale=", cstat_male(completed), "\n"))

source("B2_run_calibration.R")

# (2) age*duration interaction
fm <- update(fm, ~. -rcs(age, 4) -rcs(duration, 4) + rcs(age, 4)*rcs(duration, 4))
model_interaction <- fit.mult.impute(fm, fitter = cph, xtrans = imp, data = d, x=TRUE, y=TRUE, surv=TRUE, time.inc=5)

completed$pred <- 1 - predictSurvProb(model_interaction, completed, times=5)
print(cat("c-stat=", cstat(completed), "\nfemale=", cstat_female(completed), "\nmale=", cstat_male(completed), "\n"))
source("B2_run_calibration.R")

# (3) simplified model # 6 predictors (>90% r^2?) + female
vars <- c("rcs(duration, 4) + meds + rcs(age, 4) + smoking + rcs(log.urine_acr, 4) + rcs(hdl, 4) + female")
fm <- as.formula(paste("Surv(years, event) ~", vars))
model <- fit.mult.impute(fm, fitter = cph, xtrans = imp, data = d, x=TRUE, y=TRUE, surv=TRUE, time.inc=5)
completed$pred <- 1 - predictSurvProb(model, completed, times=5)
print(cat("c-stat=", cstat(completed), "\nfemale=", cstat_female(completed), "\nmale=", cstat_male(completed), "\n"))

ggplot(Predict(model), sepdiscrete ='vertical', nlevels=4, vnames ='names')
source("B2_run_calibration.R")

chd <- model
save(chd, file = "vfm/chd_model.Rdata")




