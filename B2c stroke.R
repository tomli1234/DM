#----------------------------------------------
# Model for stroke
#-----------------------------------------------

library(rms)
library(Hmisc)
library(pec)
library(prodlim)
library(reshape2)

load("vfm/stroke_t2_3.Rdata")
imp$formula; imp$call
d$event <- as.numeric(d$event)
d$years <- d$time / 365.25
units(d$years) <- "year"
label(d$years) <- "Survival Time"

# Model development-------------------------------------------------------------------------------------------------------
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

# (1) full model predictors for stroke, no interaction terms 
vars <- c("rcs(age, nk) + rcs(duration, nk) + rcs(log.urine_acr, nk) + rcs(bmi, nk) + rcs(haemoglobin, nk) + rcs(wbc, nk) + rcs(pulse, nk) + rcs(lr, nk) + rcs(hba1c, nk) + rcs(log.creatinine, nk) + rcs(map, nk) + rcs(triglyceride, nk) + smoking + af + cancer + pad + ckd + meds + complications + female + chd + rcs(sbp, nk) + rcs(dbp, nk) + rcs(ldl, nk) + rcs(hdl, nk) + rcs(tc, nk) + rcs(log.egfr_chinese, nk)")
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
  list(r2, frac, row.names(s[[1]]))
}

# figure 2 - r^2 approximation
r <- FUN.approx(fullmodel, completed)
# tab <- cbind(sort(r[[2]]), 1:length(r[[2]]))
# tab
# write.csv(tab, "table_r2.csv")

# 5-year risk prediction
# if we want to predict 10-year risk, change times=5 to times=10.
completed$pred <- 1 - predictSurvProb(fullmodel, completed, times=5)
table(completed$event)
# event after 5 years is setted to censored status? If so, run following code 
# completed$event[completed$event == 1 & completed$years > 5] <- 0
# table(completed$event)
# if above annotation code is executed, c-statistics for each group is slightly larger.


# Model Validation----------------------------------------------------------------------------------------
# For full model
FunDiscrimination <- function (model, bootstrap) {
  set.seed (10)
  v <- validate (model, B=bootstrap)
  print(v)
  Dxy = v[rownames(v)=="Dxy", colnames(v)=="index.corrected"]
  c.stat <- Dxy/2+0.5
  print(paste("C-statistic =", round(c.stat, 4)))
  c.stat
}
FunDiscrimination(fullmodel, 3)

# For the approximate model
Cindex_approx <- function(number_boostrap){
	c_index <- NULL
	for(i in 1:number_boostrap){
		# Bootstrap sampling
		boostrap_sample <- completed[sample(1:nrow(completed), nrow(completed), replace = TRUE), ]
		boostrap_fullmodel <- cph(fm, data=boostrap_sample, x=TRUE, y=TRUE, surv=TRUE, time.inc=5)
		r <- FUN.approx(boostrap_fullmodel, boostrap_sample)

		# Select variables that makes up 90% of R2
		cutoff <- max(which(r[[1]]>0.9))
		term_remove <- paste0(r[[3]][1:cutoff], collapse = '|')
		fm_terms <- terms(fm)
		term_remove_index <- grep(term_remove, attr(fm_terms, 'term.labels'))
		fm_approx <- drop.terms(fm_terms, term_remove_index, keep.response = TRUE)

		# Fit the approximate model
		approx_model <- cph(fm_approx, data=boostrap_sample, x=TRUE, y=TRUE, surv=TRUE, time.inc=5)

		# Calculate the C-index
		completed$pred <- predict(approx_model, newdata=completed)
		score_result <- data.frame(pred=completed$pred, years=completed$years, event=completed$event)
		c_index <- c(c_index, c(survConcordance(Surv(years, event) ~ pred, score_result)$concordance))
	}
	return(c_index)
}

Cindex_approx(3)







source("B2_run_calibration.R")

# (2) age*duration interaction
fm <- update(fm, ~. -rcs(age, 4) -rcs(duration, 4) + rcs(age, 4)*rcs(duration, 4))
model_interaction <- fit.mult.impute(fm, fitter = cph, xtrans = imp, data = d, x=TRUE, y=TRUE, surv=TRUE, time.inc=5)

completed$pred <- 1 - predictSurvProb(model_interaction, completed, times=5)
print(cat("c-stat=", cstat(completed), "\nfemale=", cstat_female(completed), "\nmale=", cstat_male(completed), "\n"))
source("B2_run_calibration.R")

# (3) simplified model # 7 predictors (>90% r^2?) + female
vars <- c("rcs(age, 4) + rcs(duration, 4) + af + meds + rcs(wbc, 4) + rcs(hba1c, 4) + rcs(hdl, 4) + female")
fm <- as.formula(paste("Surv(years, event) ~", vars))
model <- fit.mult.impute(fm, fitter = cph, xtrans = imp, data = d, x=TRUE, y=TRUE, surv=TRUE, time.inc=5)
completed$pred <- 1 - predictSurvProb(model, completed, times=5)
print(cat("c-stat=", cstat(completed), "\nfemale=", cstat_female(completed), "\nmale=", cstat_male(completed), "\n"))

ggplot(Predict(model), sepdiscrete ='vertical', nlevels=4, vnames ='names')
source("B2_run_calibration.R")

stroke <- model
save(fullmodel, stroke, file = "vfm/stroke_model.Rdata")

