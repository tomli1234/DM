# Mortality
load("E:/vfm/mortality_t2_5.Rdata")
d$event <- as.numeric(d$event)-1
# Since d$event has two status, i.e. 1 (censored indication) and 2 (event indication), transform it to 0 and 1 status as usual.
d$years <- d$time / 365.25
units(d$years) <- "year"
label(d$years) <- "Survival Time"

vars <- c("rcs(age, nk) + rcs(duration, nk) + rcs(log.urine_acr, nk) + rcs(bmi, nk) + rcs(haemoglobin, nk) + rcs(wbc, nk) + rcs(pulse, nk) + rcs(lr, nk) + rcs(hba1c, nk) + rcs(log.creatinine, nk) + rcs(map, nk) + rcs(triglyceride, nk) + smoking + af + cancer + pad + ckd + meds + complications + female + chd + stroke + rcs(sbp, nk) + rcs(dbp, nk) + rcs(ldl, nk) + rcs(hdl, nk) + rcs(tc, nk) + rcs(log.egfr_chinese, nk)")
vars <- gsub("nk", 4, vars)
fm <- as.formula(paste("Surv(years, event) ~", vars))

rmarkdown::render("C:/Users/Chao/Documents/GitHub/DM/input.Rmd", params=list(set_title="Mortality confirmed T2DM only"), output_file="C:/Users/Chao/Documents/GitHub/DM/mortality_t2.html")

# change egfr_chinese for creatinine
# change sbp, dbp for map
# change ldl, hdl, tc for lr
 


fm <- as.formula(Surv(years, event) ~ rcs(age, 4)*rcs(duration, 4) + rcs(log.urine_acr, 4)*rcs(age, 4) + rcs(bmi, 4)*rcs(age, 4) + rcs(haemoglobin, 4)*rcs(age, 4) + rcs(wbc, 4)*rcs(age, 4) + rcs(pulse, 4)*rcs(age, 4) + rcs(lr, 4)*rcs(age, 4) + rcs(hba1c, 4)*rcs(age, 4) + rcs(log.creatinine, 4)*rcs(age, 4) + rcs(map, 4)*rcs(age, 4) + rcs(triglyceride, 4)*rcs(age, 4) + smoking + af + cancer + pad + ckd + meds + complications + female + chd + stroke + rcs(sbp, 4)*rcs(age, 4) + rcs(dbp, 4)*rcs(age, 4) + rcs(ldl, 4)*rcs(age, 4) + rcs(hdl, 4)*rcs(age, 4) + rcs(tc, 4)*rcs(age, 4) + rcs(log.egfr_chinese, 4))
rmarkdown::render("C:/Users/Chao/Documents/GitHub/DM/input.Rmd", params=list(set_title="Mortality confirmed T2DM only (interaction)"), output_file="C:/Users/Chao/Documents/GitHub/DM/mortality_t2_interaction.html")


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


#reformulate(c("rcs(age, 4)*", "b*f"), fS[[2]])
fm <- update(fm, ~. +rcs(age, 4)*terms(fm))
term_remove_index <- grep(term_remove, attr(fm_terms, 'term.labels'))
fm_approx <- drop.terms(fm_terms, term_remove_index, keep.response = TRUE)
vars <- gsub("+", "*rcs(age, 4) +", vars)

fm <- as.formula(paste("Surv(years, event) ~", vars))
fm <- update(fm, ~. -rcs(age, 4)-rcs(duration, 4) + rcs(age, 4)*rcs(duration, 4))
model_interaction <- fit.mult.impute(fm, fitter = cph, xtrans = imp, data = d, x=TRUE, y=TRUE, surv=TRUE, time.inc=5)


load("E:/vfm/mortality_20170320_5.Rdata")
rmarkdown::render("C:/Users/Chao/Desktop/input.Rmd", params=list(set_title="Mortality all DM"), output_file="mortality_all_dm.html")






# CHD
vars <- c("rcs(age, nk) + rcs(duration, nk) + rcs(log.urine_acr, nk) + rcs(bmi, nk) + rcs(haemoglobin, nk) + rcs(wbc, nk) + rcs(pulse, nk) + rcs(lr, nk) + rcs(hba1c, nk) + rcs(log.creatinine, nk) + rcs(map, nk) + rcs(triglyceride, nk) + smoking + af + cancer + pad + ckd + meds + complications + female + stroke + rcs(sbp, nk) + rcs(dbp, nk) + rcs(ldl, nk) + rcs(hdl, nk) + rcs(tc, nk) + rcs(log.egfr_chinese, nk)")

load("E:/vfm/chd_t2_5.Rdata")
rmarkdown::render("C:/Users/Chao/Desktop/input.Rmd", params=list(set_title="CHD confirmed T2DM only"), output_file="chd_t2.html")

load("E:/vfm/chd_20170320_5.Rdata")
rmarkdown::render("C:/Users/Chao/Desktop/input.Rmd", params=list(set_title="CHD all DM"), output_file="chd_all_dm.html")

# Stroke
vars <- c("rcs(age, nk) + rcs(duration, nk) + rcs(log.urine_acr, nk) + rcs(bmi, nk) + rcs(haemoglobin, nk) + rcs(wbc, nk) + rcs(pulse, nk) + rcs(lr, nk) + rcs(hba1c, nk) + rcs(log.creatinine, nk) + rcs(map, nk) + rcs(triglyceride, nk) + smoking + af + cancer + pad + ckd + meds + complications + female + chd + rcs(sbp, nk) + rcs(dbp, nk) + rcs(ldl, nk) + rcs(hdl, nk) + rcs(tc, nk) + rcs(log.egfr_chinese, nk)")

load("E:/vfm/stroke_t2_3.Rdata")
rmarkdown::render("C:/Users/Chao/Desktop/input.Rmd", params=list(set_title="Stroke confirmed T2DM only (combined haemorrhagic & occulsive)"), output_file="stroke_t2_h_o.html")

load("E:/vfm/stroke_h_t2_3.Rdata")
rmarkdown::render("C:/Users/Chao/Desktop/input.Rmd", params=list(set_title="Stroke confirmed T2DM only (haemorrhagic)"), output_file="stroke_t2_h.html")

load("E:/vfm/stroke_o_t2_3.Rdata")
rmarkdown::render("C:/Users/Chao/Desktop/input.Rmd", params=list(set_title="Stroke confirmed T2DM only (occulsive)"), output_file="stroke_t2_o.html")

load("E:/vfm/stroke_20170320_5.Rdata")
rmarkdown::render("C:/Users/Chao/Desktop/input.Rmd", params=list(set_title="Stroke all DM (combined haemorrhagic & occulsive)"), output_file="stroke_all_dm_h_o.html")



