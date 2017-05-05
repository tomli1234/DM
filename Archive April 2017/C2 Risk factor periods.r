rm(list=ls())
require(data.table)

###############################################################
"THREE COHORTS" # PRE-2006, 2006-08, 2009-11, and entire
###############################################################

###############################################################
# Valuation: Characteristics
# Patients: ALL DM including type 1, age 18 or 20???
###############################################################
setwd("E:/")
load("vfm/cohort.RData")
d <- cohort[, c("serial_no", "cohort", "female", "age_entry", "duration")]
str(d)

# Observation periods (3 x 3 years)
###################################################################
table(d$cohort, exclude = NULL)
p3 <- d
p3$ob_period <- "p3"
p2 <- d
p2$ob_period <- "p2"
d <- split(d, d$cohort)
p1 <- do.call(rbind, d[c(1,2)])
p1$ob_period <- "p1"


FUN.lastvalue <- function(x, risk_factor){
	p <- data.table(x)
	p$yr <- format(p$ref_date, "%Y")
	p$ob_period <- cut(as.numeric(p$yr), c(2006, 2009, 2012, 2015), right = FALSE, labels = c("p1", "p2", "p3"))
	p <- p[,`:=`(latest=max(ref_date), count=.N), by=c("serial_no", "ob_period")]
	# count for regression dilution?
	p <- p[latest==ref_date,]
	# if multiple matches -> use mean vales
	p <- p[, list(ref_date=max(ref_date), value=mean(eval(risk_factor)), count=max(count)), by=c("serial_no", "ob_period")]
	p <- data.frame(p[, list(serial_no, ob_period, value)])
	names(p)[(names(p) %in% "value")] <- deparse(substitute(x))
	split(p, p$ob_period)
	}
	
load("Rdata/bmi_clean.RData")
bmi <- FUN.lastvalue(bmi, quote(bmi))

load("Rdata/hba1c_clean.RData") 
hba1c <- FUN.lastvalue(hba1c, quote(hba1c))

load("Rdata/bp_clean.RData") 
map <- FUN.lastvalue(map, quote(map))

load("Rdata/cholesterol_clean.RData") 
lr <- FUN.lastvalue(lr, quote(lr))

load("Rdata/creatinine_clean.RData") 
creatinine$log.creatinine <- log(creatinine$creatinine)
log.creatinine <- creatinine[, c("serial_no", "ref_date", "log.creatinine")]
log.creatinine <- FUN.lastvalue(log.creatinine, quote(log.creatinine))

load("Rdata/urine_acr_clean.RData") 
urine_acr$log.urine_acr <- log(urine_acr$urine_acr)
log.urine_acr <- urine_acr[, c("serial_no", "ref_date", "log.urine_acr")]
log.urine_acr <- FUN.lastvalue(log.urine_acr, quote(log.urine_acr))

load("Rdata/haemoglobin_clean.RData")
haemoglobin <- FUN.lastvalue(haemoglobin, quote(haemoglobin))

load("Rdata/triglyceride_clean.RData")
triglyceride <- FUN.lastvalue(triglyceride, quote(triglyceride))

load("Rdata/pulse_clean.RData")
pulse <- FUN.lastvalue(pulse, quote(pulse))

load("Rdata/wbc_clean.RData")
wbc <- FUN.lastvalue(wbc, quote(wbc))

# other risk engines
sbp <- FUN.lastvalue(sbp, quote(sbp))
# egrf
# non hdl

b1 <- Reduce(function(x, y) merge(x, y, all=T, by=c("serial_no", "ob_period")), list(bmi[[1]], hba1c[[1]], map[[1]], lr[[1]], log.creatinine[[1]], log.urine_acr[[1]], haemoglobin[[1]], triglyceride[[1]], pulse[[1]], wbc[[1]], sbp[[1]]))

b2 <- Reduce(function(x, y) merge(x, y, all=T, by=c("serial_no", "ob_period")), list(bmi[[2]], hba1c[[2]], map[[2]], lr[[2]], log.creatinine[[2]], log.urine_acr[[2]], haemoglobin[[2]], triglyceride[[2]], pulse[[2]], wbc[[2]], sbp[[2]]))

b3 <- Reduce(function(x, y) merge(x, y, all=T, by=c("serial_no", "ob_period")), list(bmi[[3]], hba1c[[3]], map[[3]], lr[[3]], log.creatinine[[3]], log.urine_acr[[3]], haemoglobin[[3]], triglyceride[[3]], pulse[[3]], wbc[[3]], sbp[[3]]))

p1 <- merge(p1, b1, all.x=T, by=c("serial_no", "ob_period"))
p2 <- merge(p2, b2, all.x=T, by=c("serial_no", "ob_period"))
p3 <- merge(p3, b3, all.x=T, by=c("serial_no", "ob_period"))

FUN.categorical <- function(x){
	p <- x
	p$value <- TRUE
	p1 <- p[p$adate<as.Date("2009-01-01"), c("serial_no", "value")]
	p2 <- p[p$adate<as.Date("2012-01-01"), c("serial_no", "value")]
	p3 <- p[p$adate<as.Date("2015-01-01"), c("serial_no", "value")]
	names(p1)[(names(p1) %in% "value")] <- deparse(substitute(x))
	names(p2)[(names(p2) %in% "value")] <- deparse(substitute(x))
	names(p3)[(names(p3) %in% "value")] <- deparse(substitute(x))
	list(p1, p2, p3)
}

load("Rdata/af.RData")
af <- FUN.categorical(af)

load("Rdata/ckd.RData")
ckd <- FUN.categorical(ckd)

load("Rdata/pad.RData")
pad <- FUN.categorical(pad)

load("Rdata/dm_complications.RData")
complications <- FUN.categorical(dm.complications)

load("Rdata/cancer.RData")
cancer <- FUN.categorical(cancer)
	
load("Rdata/dm_meds.RData")

# Diabetes medications for imputation/CUHK mortality 
#####################################################
load("Rdata/dm_meds.Rdata")
DT <- data.table(dm.meds)
# insulin (earlist date = assume permanent once started)
insulin <- DT[insulin==TRUE, list(insulin.date = min(disp_date)), serial_no]
insulin$yr <- format(insulin$insulin.date, "%Y")
insulin$p1 <- ifelse(insulin$yr < 2009, TRUE, FALSE)
insulin$p2 <- ifelse(insulin$yr < 2012, TRUE, FALSE)
insulin$p3 <- ifelse(insulin$yr < 2015, TRUE, FALSE)
library(reshape2)
insulin <- melt(insulin, id.vars=c("serial_no"), variable.name="ob_period", measure.vars=c("p1", "p2", "p3" ), value.name="insulin")
setorder(insulin, serial_no, ob_period)
insulin

# last 6 months of non-insulin prescriptions
DT <- DT[insulin==FALSE]
DT <- DT[(disp_date>as.Date("2008-06-01") & disp_date<as.Date("2009-01-01")) | (disp_date>as.Date("2011-06-01") & disp_date<as.Date("2012-01-01")) | (disp_date>as.Date("2014-06-01") & disp_date<as.Date("2015-01-01"))]

# DT <- data.table(merge(d[, c("serial_no", "entry.date")], dm.meds, by = "serial_no"))

# rename item codes
# non-insulin (HA table 2)
DT$item_cd <- gsub("ACAR01|ACAR02", "ACARBOSE", DT$item_cd)
DT$item_cd <- gsub("EXEN01|EXEN02|EXEN03|EXEN04|S00623|S00624", "EXENATIDE", DT$item_cd)
DT$item_cd <- gsub("GLIB01|GLIB02", "GLIBENCLAMIDE", DT$item_cd)
DT$item_cd <- gsub("GLIC01|GLIC02|GLIC03|S00089", "GLICLAZIDE", DT$item_cd)
DT$item_cd <- gsub("GLIM01|GLIM02|GLIM03|S00248|S00249", "GLIMEPIRIDE", DT$item_cd)
DT$item_cd <- gsub("GLIP01|GLIP02", "GLIPIZIDE", DT$item_cd)
DT$item_cd <- gsub("GLUC01|GLUC37", "GLUCAGON", DT$item_cd)
DT$item_cd <- gsub("LINA01|S00893", "LINAGLIPTIN", DT$item_cd)
DT$item_cd <- gsub("METF01|METF02|METF03|METF04|S00466", "METFORMIN", DT$item_cd)
DT$item_cd <- gsub("PIOG01|PIOG02|S00016|S00599", "PIOGLITAZONE", DT$item_cd)
DT$item_cd <- gsub("S00790|SAXA01", "SAXAGLIPTIN", DT$item_cd)
DT$item_cd <- gsub("S00600|S00602|S00675|SITA02|SITA03|SITA06", "SITAGLIPTIN", DT$item_cd)
DT$item_cd <- gsub("TOLB01", "TOLBUTAMIDE", DT$item_cd)
DT$item_cd <- gsub("S00731|VILD01", "VILDAGLIPTIN", DT$item_cd)

DT$yr <- format(DT$disp_date, "%Y")
DT$ob_period <- cut(as.numeric(DT$yr), c(2006, 2009, 2012, 2015), right = FALSE, labels = c("p1", "p2", "p3"))

# no of non-insulin agents
non.insulin <- DT[, list(other.meds = .N), by = list(serial_no, ob_period)]

insulin <- data.frame(insulin)	
non.insulin <- data.frame(non.insulin)
meds <- merge(insulin[, c("serial_no", "ob_period", "insulin")], non.insulin[, c("serial_no", "ob_period", "other.meds")], all = T, by = c("serial_no", "ob_period"))
meds$meds <- ifelse(meds$other.meds>=3, "3+", meds$other.meds)
meds$meds[meds$insulin == TRUE] <- "insulin"
rm(dm.meds, DT)

meds <- meds[, c("serial_no", "ob_period", "meds")]
meds <- split(meds, meds$ob_period)

c1 <- Reduce(function(x, y) merge(x, y, all=T, by=c("serial_no")), list(af[[1]], ckd[[1]], pad[[1]], complications[[1]], cancer[[1]], meds[[1]]))
c2 <- Reduce(function(x, y) merge(x, y, all=T, by=c("serial_no")), list(af[[2]], ckd[[2]], pad[[2]], complications[[2]], cancer[[2]], meds[[2]]))
c3 <- Reduce(function(x, y) merge(x, y, all=T, by=c("serial_no")), list(af[[3]], ckd[[3]], pad[[3]], complications[[3]], cancer[[3]], meds[[3]]))

p1 <- merge(p1, c1, all.x=T, by=c("serial_no", "ob_period"))
p2 <- merge(p2, c2, all.x=T, by=c("serial_no", "ob_period"))
p3 <- merge(p3, c3, all.x=T, by=c("serial_no", "ob_period"))

# recode missing categorical as zero
categorical <- names(c1)[-1]

p1[, categorical][is.na(p1[, categorical])] <- 0
p2[, categorical][is.na(p2[, categorical])] <- 0
p3[, categorical][is.na(p3[, categorical])] <- 0

save(p1, p2, p3, file="vfm/periods.Rdata")

load("vfm/periods.Rdata")
#######################################################
# Multiple Imputation using aregImpute
require(Hmisc);
'TIME SERIES DATA: Age and duration, Impute combined dataset together???'

p1 <- subset(p1, select=-c(cohort, ob_period))
p2 <- subset(p2, select=-c(cohort, ob_period))
p3 <- subset(p3, select=-c(cohort, ob_period))

fm <- as.formula(paste(" ~ ", paste(names(p1), collapse = "+")))
fm <- update(fm, ~. -serial_no ) 
strt <- proc.time()
impute1 <- aregImpute(fm, n.impute = 3, data = p1, nk = 4)
impute2 <- aregImpute(fm, n.impute = 3, data = p2, nk = 4)
impute3 <- aregImpute(fm, n.impute = 3, data = p3, nk = 4)
print(proc.time() - strt)

save(impute1, impute2, impute3, file="vfm/periods_impute.Rdata")

imp1 <- impute.transcan(impute1, imputation=1, data=p1, list.out=TRUE, pr=FALSE, check=FALSE) 
imp2 <- impute.transcan(impute1, imputation=1, data=p2, list.out=TRUE, pr=FALSE, check=FALSE) 
imp3 <- impute.transcan(impute1, imputation=1, data=p3, list.out=TRUE, pr=FALSE, check=FALSE) 


p1[names(imp1)] <- imp1

# predicted scores
p1$pred <- (predict(model, newdata=p1))
p2$pred <- (predict(model, newdata=p2))
p3$pred <- (predict(model, newdata=p3))

# predicted survival at 5-years
library(pec)
p1$pred <- 1 - predictSurvProb(model, p1, times=5)
p2$pred <- 1 - predictSurvProb(model, p2, times=5)
p3$pred <- 1 - predictSurvProb(model, p3, times=5)


# Absolute Risk Table (requires actual age & duration)




