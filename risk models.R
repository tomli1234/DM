# Mortality
vars <- c("rcs(age, nk) + rcs(duration, nk) + rcs(log.urine_acr, nk) + rcs(bmi, nk) + rcs(haemoglobin, nk) + rcs(wbc, nk) + rcs(pulse, nk) + rcs(lr, nk) + rcs(hba1c, nk) + rcs(log.creatinine, nk) + rcs(map, nk) + rcs(triglyceride, nk) + smoking + af + cancer + pad + ckd + meds + complications + female + chd + stroke + rcs(sbp, nk) + rcs(dbp, nk) + rcs(ldl, nk) + rcs(hdl, nk) + rcs(tc, nk) + rcs(log.egfr_chinese, nk)")

load("E:/vfm/mortality_t2_5.Rdata")
rmarkdown::render("C:/Users/Chao/Desktop/input.Rmd", params=list(set_title="Mortality confirmed T2DM only"), output_file="mortality_t2.html")

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



