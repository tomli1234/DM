#----------------------------------------------
# Model analysis for Chd
# Zhiqiang, March 2017
#-----------------------------------------------

####load library
library(rms); 
library(Hmisc)
library(pec)
library(prodlim)
library(reshape2)
rm(list=ls())

#load data
load("chd_20170220.Rdata")
imp$formula; imp$call
d$event <- as.numeric(d$event)
d$years <- d$time / 365.25
units(d$years) <- "year"
label(d$years) <- "Survival Time"

#imputing missing data
completed <- d
imputed <- impute.transcan(imp, imputation=1, data=d, list.out=TRUE, pr=FALSE, check=FALSE) 
completed[names(imputed)] <- imputed
#Since d$event and completed$event have two status, i.e., 
#1 (censored indication) and 2 (event indication), 
#transform it to 0 and 1 status as usual.
completed$event <- completed$event-1
d$event <- d$event-1
dd <<- datadist(completed); options(datadist ="dd")

#20 predictors for chd, no interaction terms are conisdered 
nk <- 4
vars <- c("rcs(age, nk) + rcs(duration, nk) + rcs(log.urine_acr, nk) + rcs(bmi, nk) + rcs(haemoglobin, nk) + rcs(wbc, nk) + rcs(pulse, nk) + rcs(lr, nk) + rcs(hba1c, nk) + rcs(log.creatinine, nk) + rcs(map, nk) + rcs(triglyceride, nk) +smoking + af + cancer + pad + ckd + meds + complications + female")
vars <- gsub("nk", nk, vars)
fm <- as.formula(paste("Surv(years, event) ~", vars))
chd <- fit.mult.impute(fm, fitter = cph, xtrans = imp, data = d, x=TRUE, y=TRUE, surv=TRUE, time.inc=5)

#estimation curve for variable
ggplot(Predict(chd), sepdiscrete ='vertical', nlevels=4, vnames ='names')
#importance of each variable
plot(anova(chd))
#estimation result
chd

#5-year risk prediction
#if we want to predict 10-year risk, change times=5 to times=10.
completed$pred <- 1 - predictSurvProb(chd, completed, times=5)
table(completed$event)
#event after 5 years is setted to censored status? If so, run following code 
#completed$event[completed$event == 1 & completed$years > 5] <- 0
#table(completed$event)
#if above annotation code is executed, c-statistics for each group is slightly large.

#C-statistics for male group
cstat_male <- function(completed){
  score_result <- cbind(completed[completed$female==F, "pred"],completed[completed$female==F, "years"], completed[completed$female==F, "event"],completed[completed$female==F, "pred"])
  cstat=survConcordance(Surv(score_result[,2], score_result[,3])~score_result[,1])$concordance
  cstat
}
cstat_male(completed)

#C-statistic for female group
cstat_female<-function(completed){
  score_result <- cbind(completed[completed$female==T, "pred"], completed[completed$female==T, "years"], completed[completed$female==T, "event"],completed[completed$female==T, "pred"])
  cstat <- survConcordance(Surv(score_result[,2], score_result[,3])~ score_result[,1])$concordance
  cstat
}
cstat_female(completed)

#C-statistic for overall
cstat<-function(completed){
  score_result <- cbind(completed$pred, completed$years, completed$event, completed$pred)
  cstats <- survConcordance(Surv(score_result[,2], score_result[,3])~ score_result[,1])$concordance
  cstats
}
cstat(completed)

#predict 10 group's risk probabilities
#cut data.frame into 10 group
FUN.deciles_mean_pred <- function(data) {
  data <- data[order(data$pred),]
  data$group <- cut2(data$pred, g=10)
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

#multiple plot function
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL){
  library(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#predict survival probability and observed probability for overall
dat <-FUN.deciles_mean_pred(completed)
dat.melt <- melt(dat, id.vars = 'Deciles', variable.name = 'Group', value.name = 'Risk')
plot_all<-function(dat.melt){
  ggplot(dat.melt, aes(x = Deciles, y = Risk*100, color = Group)) + geom_point(size = 4) +  xlab("Tenth of predicted risk") + ylab("5 year risk (%)") + theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) + scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks=seq(0, 100, 10)) + scale_color_manual(values = c("observed" = 'skyblue1','pred' = 'darkblue'), labels=c("observed" = 'Observed','pred' = 'Predicted')) + theme(legend.position="bottom", legend.title=element_blank(), legend.text = element_text(size = 12, face = "bold"))
}
plot_all(dat.melt)

#predict survival probability and observed probability for male group
dat <- FUN.deciles_mean_pred(completed[completed$female==0,])
dat.melt <- melt(dat, id.vars = 'Deciles', variable.name = 'Group', value.name = 'Risk')
plot_male<-function(dat.melt){
  p1 <- ggplot(dat.melt, aes(x = Deciles, y = Risk*100, color = Group)) + geom_point(size = 4) +  xlab("Tenth of predicted risk") + ylab("5 year risk (%)") + theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) + scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks=seq(0, 100, 10)) + scale_color_manual(values = c("observed" = 'skyblue1','pred' = 'darkblue'), labels=c("observed" = 'Observed','pred' = 'Predicted')) + theme(legend.position="bottom", legend.title=element_blank(), legend.text = element_text(size = 12, face = "bold")) + ggtitle("Male") 
}
p1<-plot_male(dat.melt)

#predict survival probability and observed probability for female group
dat <- FUN.deciles_mean_pred(completed[completed$female==1,])
dat.melt <- melt(dat, id.vars = 'Deciles', variable.name = 'Group', value.name = 'Risk')
plot_female<-function(dat.melt){
  p2 <- ggplot(dat.melt, aes(x = Deciles, y = Risk*100, color = Group)) + geom_point(size = 4) +  xlab("Tenth of predicted risk") + ylab("5 year risk (%)") + theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) + scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks=seq(0, 100, 10)) + scale_color_manual(values = c("observed" = 'pink','pred' = 'darkred'), labels=c("observed" = 'Observed','pred' = 'Predicted')) + theme(legend.position="bottom", legend.title=element_blank(), legend.text = element_text(size = 12, face = "bold")) + ggtitle("Female") 
}
p2<-plot_female(dat.melt)
multiplot(p1, p2, cols = 1)