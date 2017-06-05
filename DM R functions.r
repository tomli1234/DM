library(doParallel)

# C-index for full model
FunDiscrimination <- function (model, bootstrap) {
      set.seed (10)
      v <- validate (model, B=bootstrap)
      print(v)
      Dxy = v[rownames(v)=="Dxy", colnames(v)=="index.corrected"]
      c.stat <- Dxy/2+0.5
      print(paste("C-statistic =", round(c.stat, 4)))
      c.stat
}

# cut data.frame into 10 group
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

# Calibration plot
plot_all_dots <- function(dat.melt, colours = c('skyblue1','darkblue')){
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


# Model approximation
FUN.approx <- function(mod, completedata) {
      lp <- predict(mod)
      a <- ols(as.formula(paste0("lp ~ ", paste(labels(terms(fm)), collapse = "+"))), sigma=1, data=completedata)
      s <- fastbw(a, aic=1000000)
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
      list(r2=c(1, r2), frac=c(1, frac), predictor=row.names(s[[1]]))
}

# Bootstrap Validation for the approximate model
Cindex_approx <- function(number_bootstrap, formula = fm, data = completed){
      bootstrap <- function(){
            library(rms)
            library(Hmisc)
            library(pec)
            library(prodlim)
            library(reshape2)
            
            # Bootstrap sampling
            boostrap_sample <- completed[sample(1:nrow(completed),
                                                nrow(completed),
                                                replace = TRUE), ]
            boostrap_fullmodel <- cph(fm, data=boostrap_sample,
                                      x=TRUE, y=TRUE, surv=TRUE, time.inc=5)
            r <- FUN.approx(boostrap_fullmodel, boostrap_sample)

            # Select variables that makes up 90% of R2
            cutoff <- max(which(r[[1]]>0.9))
            term_remove <- paste0(r[[3]][1:cutoff], collapse = '|')
            fm_terms <- terms(fm)
            term_remove_index <- grep(term_remove, attr(fm_terms, 'term.labels'))
            fm_approx <- drop.terms(fm_terms, term_remove_index, keep.response = TRUE)

            # Fit the approximate model
            approx_model <- cph(fm_approx, data=boostrap_sample,
                                x=TRUE, y=TRUE, surv=TRUE, time.inc=5)

            # Calculate the C-index
            completed$pred <- predict(approx_model, newdata=completed)
            score_result <- data.frame(pred=completed$pred,
                                       years=completed$years,
                                       event=completed$event)
            c_index <- c(survConcordance(Surv(years, event) ~ pred, score_result)$concordance)
            return(c_index)
      }
      
      # Use doParallel to speed up the process ---------------------------
      no_cores <- detectCores() - 1
      registerDoParallel(cores = no_cores)
      cl <- makeCluster(no_cores, type="PSOCK")
      clusterExport(cl, list("completed","fm","FUN.approx"))
      c_index <- parLapply(cl, 1:number_bootstrap, function(x) bootstrap())
      stopCluster(cl)
      
      return(do.call(rbind,c_index))

}


