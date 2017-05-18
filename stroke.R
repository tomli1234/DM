# Stroke
load("E:/vfm/stroke_t2_3.Rdata"); folder <- "stroke/t2 dm"
# load("E:/vfm/stroke_h_t2_3.Rdata"); folder <- "stroke_h/t2 dm"
# load("E:/vfm/stroke_o_t2_3.Rdata"); folder <- "stroke_o/t2 dm"
# load("E:/vfm/stroke_20170320_5.Rdata"); folder <- "stroke/all dm"

setwd(paste0("C:/Users/Chao/Dropbox (Personal)/Diabetes/Results/", folder))
d$event <- as.numeric(d$event)-1
# Since d$event has two status, i.e. 1 (censored indication) and 2 (event indication), transform it to 0 and 1 status as usual.
d$years <- d$time / 365.25
units(d$years) <- "year"; label(d$years) <- "Survival Time"

fm <- as.formula(Surv(years, event) ~ rcs(age, 4) + rcs(duration, 4) + rcs(log.urine_acr, 4) + rcs(bmi, 4) + rcs(haemoglobin, 4) + rcs(wbc, 4) + rcs(pulse, 4) + rcs(lr, 4) + rcs(hba1c, 4) + rcs(log.creatinine, 4) + rcs(map, 4) + rcs(triglyceride, 4) + smoking + af + cancer + pad + ckd + meds + complications + female + chd + rcs(sbp, 4) + rcs(dbp, 4) + rcs(ldl, 4) + rcs(hdl, 4) + rcs(tc, 4) + rcs(log.egfr_chinese, 4))
fm1 <- fm

rmarkdown::render("C:/Users/Chao/Documents/GitHub/DM/input.Rmd", params=list(set_title=paste(folder, "all predictors")), output_file=paste(folder, "stroke.html"))

# change egfr_chinese, creatinine (+/-urine acr)
fm <- update(fm1, ~. -rcs(log.egfr_chinese, 4))
rmarkdown::render("C:/Users/Chao/Documents/GitHub/DM/input.Rmd", params=list(set_title=paste(folder, "creatinine")), output_file=paste(folder, "creatinine.html"))

fm <- update(fm1, ~. -rcs(creatinine, 4))
rmarkdown::render("C:/Users/Chao/Documents/GitHub/DM/input.Rmd", params=list(set_title=paste(folder, "egfr")), output_file=paste(folder, "egfr.html"))

fm <- update(fm1, ~. -rcs(log.egfr_chinese, 4) -rcs(creatinine, 4))
rmarkdown::render("C:/Users/Chao/Documents/GitHub/DM/input.Rmd", params=list(set_title=paste(folder, "uacr")), output_file=paste(folder, "uacr.html"))

fm <- update(fm1, ~. -rcs(creatinine, 4) -rcs(log.urine_acr, 4) )
rmarkdown::render("C:/Users/Chao/Documents/GitHub/DM/input.Rmd", params=list(set_title=paste(folder, "no uacr")), output_file=paste(folder, "no_uacr.html"))

# change sbp, dbp, map
fm <- update(fm, ~. -rcs(creatinine, 4) -rcs(dbp, 4) -rcs(map, 4))
rmarkdown::render("C:/Users/Chao/Documents/GitHub/DM/input.Rmd", params=list(set_title=paste(folder, "sbp")), output_file=paste(folder, "sbp.html"))

fm <- update(fm, ~. -rcs(creatinine, 4) -rcs(sbp, 4) -rcs(map, 4))
rmarkdown::render("C:/Users/Chao/Documents/GitHub/DM/input.Rmd", params=list(set_title=paste(folder, "dbp")), output_file=paste(folder, "dbp.html"))

fm <- update(fm, ~. -rcs(creatinine, 4) -rcs(sbp, 4) -rcs(dbp, 4))
rmarkdown::render("C:/Users/Chao/Documents/GitHub/DM/input.Rmd", params=list(set_title=paste(folder, "map")), output_file=paste(folder, "map.html"))

# change lr, ldl, hdl, tc
fm <- update(fm, ~. -rcs(creatinine, 4) -rcs(dbp, 4) -rcs(map, 4) -rcs(hdl, 4) -rcs(ldl, 4) -rcs(tc, 4))
rmarkdown::render("C:/Users/Chao/Documents/GitHub/DM/input.Rmd", params=list(set_title=paste(folder, "lr")), output_file=paste(folder, "lr.html"))

fm <- update(fm, ~. -rcs(creatinine, 4) -rcs(dbp, 4) -rcs(map, 4) -rcs(lr, 4) -rcs(ldl, 4) -rcs(tc, 4))
rmarkdown::render("C:/Users/Chao/Documents/GitHub/DM/input.Rmd", params=list(set_title=paste(folder, "hdl")), output_file=paste(folder, "hdl.html"))

fm <- update(fm, ~. -rcs(creatinine, 4) -rcs(dbp, 4) -rcs(map, 4) -rcs(lr) -rcs(hdl, 4) -rcs(tc, 4))
rmarkdown::render("C:/Users/Chao/Documents/GitHub/DM/input.Rmd", params=list(set_title=paste(folder, "hdl")), output_file=paste(folder, "ldl.html"))

fm <- update(fm, ~. -rcs(creatinine, 4) -rcs(dbp, 4) -rcs(map, 4) -rcs(lr) -rcs(hdl, 4) -rcs(ldl, 4))
rmarkdown::render("C:/Users/Chao/Documents/GitHub/DM/input.Rmd", params=list(set_title=paste(folder, "tc")), output_file=paste(folder, "tc.html"))


# Interactions (full model)
fm <- as.formula(Surv(years, event) ~ rcs(age, 4)*rcs(duration, 4) + rcs(log.urine_acr, 4)*rcs(age, 4) + rcs(bmi, 4)*rcs(age, 4) + rcs(haemoglobin, 4)*rcs(age, 4) + rcs(wbc, 4)*rcs(age, 4) + rcs(pulse, 4)*rcs(age, 4) + rcs(lr, 4)*rcs(age, 4) + rcs(hba1c, 4)*rcs(age, 4) + rcs(log.creatinine, 4)*rcs(age, 4) + rcs(map, 4)*rcs(age, 4) + rcs(triglyceride, 4)*rcs(age, 4) + smoking + af + cancer + pad + ckd + meds + complications + female + chd + rcs(sbp, 4)*rcs(age, 4) + rcs(dbp, 4)*rcs(age, 4) + rcs(ldl, 4)*rcs(age, 4) + rcs(hdl, 4)*rcs(age, 4) + rcs(tc, 4)*rcs(age, 4) + rcs(log.egfr_chinese, 4))

# Interactions (simplified model)
fm <- as.formula(Surv(years, event) ~ rcs(age, 4)*rcs(duration, 4) + rcs(age, 4)*rcs(log.urine_acr, 4) + rcs(age, 4)*rcs(haemoglobin, 4) + rcs(age, 4)*female)

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

