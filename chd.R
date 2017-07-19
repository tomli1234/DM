# CHD
# load("E:/vfm/CHD_t2_5.Rdata"); folder <- "CHD/t2 dm"
load("E:/vfm/CHD_20170320_5.Rdata"); folder <- "CHD/all dm "

d$event <- as.numeric(d$event)-1
# Since d$event has two status, i.e. 1 (censored indication) and 2 (event indication), transform it to 0 and 1 status as usual.
d$years <- d$time / 365.25
units(d$years) <- "year"; label(d$years) <- "Survival Time"

fm <- as.formula(Surv(years, event) ~ rcs(age, 4) + rcs(duration, 4) + rcs(log.urine_acr, 4) + rcs(bmi, 4) + rcs(haemoglobin, 4) + rcs(wbc, 4) + rcs(pulse, 4) + rcs(lr, 4) + rcs(hba1c, 4) + rcs(log.creatinine, 4) + rcs(map, 4) + rcs(triglyceride, 4) + smoking + af + cancer + pad + ckd + meds + complications + female + stroke + rcs(sbp, 4) + rcs(dbp, 4) + rcs(ldl, 4) + rcs(hdl, 4) + rcs(tc, 4) + rcs(log.egfr_chinese, 4))
fm1 <- fm

filename <- "all predictors"
rmarkdown::render("C:/Users/Chao/Documents/GitHub/DM/input.Rmd", params=list(set_title=paste(folder, filename)), output_file=paste0(folder, filename, ".html"))

# change egfr_chinese, creatinine (+/-urine acr)
rm(list=setdiff(ls(), c("d", "imp", "fm1", "folder")))
fm <- update(fm1, ~. -rcs(log.egfr_chinese, 4))

rm(list=setdiff(ls(), c("d", "imp", "fm1", "folder")))
fm <- update(fm1, ~. -rcs(log.creatinine, 4))

rm(list=setdiff(ls(), c("d", "imp", "fm1", "folder")))
fm <- update(fm1, ~. -rcs(log.egfr_chinese, 4) -rcs(log.creatinine, 4))

rm(list=setdiff(ls(), c("d", "imp", "fm1", "folder")))
fm <- update(fm1, ~. -rcs(log.creatinine, 4) -rcs(log.urine_acr, 4) )

# change sbp, dbp, map
rm(list=setdiff(ls(), c("d", "imp", "fm1", "folder")))
fm <- update(fm1, ~. -rcs(log.egfr_chinese, 4) -rcs(dbp, 4) -rcs(map, 4))

rm(list=setdiff(ls(), c("d", "imp", "fm1", "folder")))
fm <- update(fm1, ~. -rcs(log.egfr_chinese, 4) -rcs(sbp, 4) -rcs(map, 4))

rm(list=setdiff(ls(), c("d", "imp", "fm1", "folder")))
fm <- update(fm1, ~. -rcs(log.egfr_chinese, 4) -rcs(sbp, 4) -rcs(dbp, 4))

# change lr, ldl, hdl, tc
rm(list=setdiff(ls(), c("d", "imp", "fm1", "folder")))
filename <- "creatinine dbp lr"
fm <- update(fm1, ~. -rcs(log.egfr_chinese, 4) -rcs(sbp, 4) -rcs(map, 4) -rcs(hdl, 4) -rcs(ldl, 4) -rcs(tc, 4))
rmarkdown::render("C:/Users/Chao/Documents/GitHub/DM/input.Rmd", params=list(set_title=paste(folder, filename)), output_file=paste0(folder, filename, ".html"))

rm(list=setdiff(ls(), c("d", "imp", "fm1", "folder")))
fm <- update(fm1, ~. -rcs(log.egfr_chinese, 4) -rcs(sbp, 4) -rcs(map, 4) -rcs(lr, 4) -rcs(ldl, 4) -rcs(tc, 4))

rm(list=setdiff(ls(), c("d", "imp", "fm1", "folder")))
fm <- update(fm1, ~. -rcs(log.egfr_chinese, 4) -rcs(sbp, 4) -rcs(map, 4) -rcs(lr, 4) -rcs(hdl, 4) -rcs(tc, 4))

rm(list=setdiff(ls(), c("d", "imp", "fm1", "folder")))
fm <- update(fm1, ~. -rcs(log.egfr_chinese, 4) -rcs(dbp, 4) -rcs(map, 4) -rcs(lr, 4) -rcs(hdl, 4) -rcs(ldl, 4))


# Interactions (full model)
fm <- as.formula(Surv(years, event) ~ rcs(age, 4)*rcs(duration, 4) + rcs(log.urine_acr, 4)*rcs(age, 4) + rcs(bmi, 4)*rcs(age, 4) + rcs(haemoglobin, 4)*rcs(age, 4) + rcs(wbc, 4)*rcs(age, 4) + rcs(pulse, 4)*rcs(age, 4) + rcs(lr, 4)*rcs(age, 4) + rcs(hba1c, 4)*rcs(age, 4) + rcs(log.creatinine, 4)*rcs(age, 4) + rcs(map, 4)*rcs(age, 4) + rcs(triglyceride, 4)*rcs(age, 4) + smoking + af + cancer + pad + ckd + meds + complications + female + stroke + rcs(sbp, 4)*rcs(age, 4) + rcs(dbp, 4)*rcs(age, 4) + rcs(ldl, 4)*rcs(age, 4) + rcs(hdl, 4)*rcs(age, 4) + rcs(tc, 4)*rcs(age, 4) + rcs(log.egfr_chinese, 4))

rmarkdown::render("C:/Users/Chao/Documents/GitHub/DM/input.Rmd", params=list(set_title=paste(folder, "interactions(full model)")), output_file=paste(folder, "interactions.html"))

# Interactions (simplified model)
completed <- d
imputed <- impute.transcan(imp, imputation=1, data=d, list.out=TRUE, pr=FALSE, check=FALSE) 
completed[names(imputed)] <- imputed
dd <<- datadist(completed); options(datadist ="dd")

fm <- as.formula(Surv(years, event) ~ rcs(age, 4)*rcs(duration, 4) + rcs(age, 4)*rcs(log.urine_acr, 4) + rcs(age, 4)*female)

FunDiscrimination <- function (model, bootstrap) {
  set.seed (10)
  v <- validate (model, B=bootstrap)
  print(v)
  Dxy = v[rownames(v)=="Dxy", colnames(v)=="index.corrected"]
  c.stat <- Dxy/2+0.5
  print(paste("C-statistic =", round(c.stat, 4)))
  c.stat
}
fullmodel <- fit.mult.impute(fm, fitter = cph, xtrans = imp, data = d, x=TRUE, y=TRUE, surv=TRUE, time.inc=5)
c.stat <- FunDiscrimination(fullmodel, 3)

