library(doParallel)

# Model approximation----------------------------------------------------------------
FUN.approx <- function(mod, completedata, fm) {
      lp <- predict(mod)
      a <- ols(as.formula(paste0("lp ~ ", paste(labels(terms(fm)), collapse = "+"))), sigma=1, data=completedata)
      s <- fastbw(a, aic = 10000000000)
      betas <- s$Coefficients
      X <- model.matrix(fm, data = completedata)
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
      
      r2 <- c(1, r2)
      predictor <- row.names(s[[1]])
      frac <-  c(1, frac)
      
      # Select variables that makes up 90% of R2
      cutoff <- max(which(r2 > 0.9))
      term_remove <- paste0(predictor[1: (cutoff - 1)], collapse = '|')
      fm_terms <- terms(fm)
      term_remove_index <- grep(term_remove, attr(fm_terms, 'term.labels'))
      fm_approx <- drop.terms(fm_terms, term_remove_index, keep.response = TRUE)
      
      # Fit the approximate model
      approx_model <- cph(fm_approx, data=completedata,
                          x=TRUE, y=TRUE, surv=TRUE, time.inc=5)
      
      return(list(r2 = r2,
                  frac = frac,
                  predictor = predictor,
                  approx_model = approx_model))
}

# Calibration------------------------------------------------------------------------
# Plot
plot_cal_dots <- function(dat.melt, colours = c('skyblue1','darkblue')){
      ggplot(dat.melt, aes(x = Deciles, y = Risk*100, color = Group)) + 
            geom_point(size = 4) +  
            xlab("Tenth of predicted risk") + 
            ylab("5 year risk (%)") + 
            theme(axis.title.x = element_text(size=14), 
                  axis.title.y = element_text(size=14), 
                  axis.text.y = element_text(size=12)) + 
            scale_x_continuous(breaks = NULL) + 
            scale_y_continuous(breaks=seq(0, 100, 10)) + 
            scale_color_manual(values = c("observed" = colours[1],'pred' = colours[2]),
                               labels=c("observed" = 'Observed','pred' = 'Predicted')) +
            theme(legend.position="bottom", 
                  legend.title=element_blank(), 
                  legend.text = element_text(size = 12, face = "bold"))
}

plot_cal_bar <- function(dat.melt, colours = c('skyblue1','darkblue')){
      ggplot(dat.melt, aes(x = factor(Deciles), y = Risk*100, fill = Group)) + 
            geom_bar(stat = 'identity', position = 'dodge', width = 0.7) +  
            xlab("Tenth of predicted risk") + 
            ylab("5 year risk (%)") + 
            theme_classic()+
            theme(axis.title.x = element_text(size=14), 
                  axis.title.y = element_text(size=14), 
                  axis.text.y = element_text(size=12)) + 
            scale_x_discrete() + 
            scale_y_continuous(breaks=seq(0, 100, 10)) + 
            scale_fill_manual(values = c("observed" = colours[1],'pred' = colours[2]),
                              labels=c("observed" = 'Observed','pred' = 'Predicted')) +
            theme(legend.position="bottom", 
                  legend.title=element_blank(), 
                  legend.text = element_text(size = 12, face = "bold"),
                  panel.grid.major.y = element_line( size=.1, color="grey"))
      
}

# Calibration for full model (using rms:calibrate)
Calibrate_full <- function (model, bootstrap) {
      cal.decile <- calibrate(model, u=5, B=bootstrap, cmethod="KM", m=(nrow(completed)/10))
      dat.melt <- melt(cal.decile[,c('mean.predicted', 'KM.corrected')])
      names(dat.melt) <- c('Deciles','Group','Risk')
      dat.melt[,'Risk'] <- 1 - dat.melt[,'Risk']
      levels(dat.melt[,'Group']) <- c('pred','observed')
      dat.melt$Deciles <- 11 - dat.melt$Deciles
      p <- plot_cal_bar(dat.melt) 
      return(p)
}

# Discrimnation-----------------------------------------------------------------------------
# C-index for full model (using rms::validate)
Cindex_full <- function (model, bootstrap) {
      set.seed (10)
      v <- validate (model, B=bootstrap)
      print(v)
      Dxy = v[rownames(v)=="Dxy", colnames(v)=="index.corrected"]
      c.stat <- Dxy/2+0.5
      print(paste("C-statistic =", round(c.stat, 4)))
      c.stat
}

# Bootstrap Validation for the approximate model
valid_approx <- function(number_bootstrap, fm = fm, completedata = completed){
      bootstrap <- function(){
            library(rms)
            library(Hmisc)
            library(pec)
            library(prodlim)
            library(reshape2)
            
            # Bootstrap sampling
            boostrap_sample <- completedata[sample(1:nrow(completedata),
                                                nrow(completedata),
                                                replace = TRUE), ]
            boostrap_fullmodel <- cph(fm, data=boostrap_sample,
                                      x=TRUE, y=TRUE, surv=TRUE, time.inc=5)
            approx <- FUN.approx(boostrap_fullmodel, boostrap_sample, fm)
            approx_model <- approx$approx_model
            
            # Calculate the C-index
            completedata$pred <- predict(approx_model, newdata=completedata)
            score_result <- data.frame(pred=completedata$pred,
                                       years=completedata$years,
                                       event=completedata$event)
            c_index <- c(survConcordance(Surv(years, event) ~ pred, score_result)$concordance)
            
            # Calibration
            pred <- predictSurvProb(approx_model, completedata, times=5)
            S <- Surv(completedata$years, completedata$event)
            groupkm(pred, S, g=10, u=5, pl=TRUE)
            
            completedata$pred <- 1 - predictSurvProb(approx_model, completedata, times=5)
            dat <- FUN.deciles_mean_pred(completedata)
            dat.melt <- melt(dat, id.vars = 'Deciles', variable.name = 'Group', value.name = 'Risk')
            # p1 <- plot_cal_bar(dat.melt)
            
            return(list(c_index = c_index,
                        calibration = dat.melt)) 
      }
      
      completedata <<- completedata
      
      # Use doParallel to speed up the process ---------------------------
      no_cores <- detectCores() - 1
      registerDoParallel(cores = no_cores)
      cl <- makeCluster(no_cores, type="PSOCK")
      clusterExport(cl, list("completedata","fm","FUN.approx","FUN.deciles_mean_pred","plot_cal_bar"))
      boots_results <- parLapply(cl, 1:number_bootstrap, function(x) bootstrap())
      stopCluster(cl)
      
      return(boots_results)
      
}



# cut data.frame into 10 group
FUN.deciles_mean_pred <- function(data) {
      data <- data[order(data$pred),]
      data$group <- cut2(data$pred, g = 10)
      data$event <- as.numeric(data$event)
      levels(data$group) <- c(1:10)
      observed <- rep(NA, 10)
      for (i in 1:10) {
            km <- survfit(Surv(years, event) ~ 1, data = data[data$group == i, ])
            survest <- stepfun(km$time, c(1, km$surv))
            observed[i] <- 1-survest(5) # risk at 5-years
      }
      mean_pred <- tapply(data$pred, data$group, mean)
      data.frame(observed, pred = mean_pred, Deciles=1:10)
}

# multiple plot function
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






