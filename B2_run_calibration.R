library(reshape2)

# predict 10 group's risk probabilities
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

# predict survival probability and observed probability for overall
dat <-FUN.deciles_mean_pred(completed)
dat.melt <- melt(dat, id.vars = 'Deciles', variable.name = 'Group', value.name = 'Risk')
plot_all<-function(dat.melt){
  ggplot(dat.melt, aes(x = Deciles, y = Risk*100, color = Group)) + geom_point(size = 4) +  xlab("Tenth of predicted risk") + ylab("5 year risk (%)") + theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) + scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks=seq(0, 100, 10)) + scale_color_manual(values = c("observed" = 'skyblue1','pred' = 'darkblue'), labels=c("observed" = 'Observed','pred' = 'Predicted')) + theme(legend.position="bottom", legend.title=element_blank(), legend.text = element_text(size = 12, face = "bold"))
}
plot_all(dat.melt)

# predict survival probability and observed probability for male group
dat <- FUN.deciles_mean_pred(completed[completed$female==F,])
dat.melt <- melt(dat, id.vars = 'Deciles', variable.name = 'Group', value.name = 'Risk')
plot_male<-function(dat.melt){
  p1 <- ggplot(dat.melt, aes(x = Deciles, y = Risk*100, color = Group)) + geom_point(size = 4) +  xlab("Tenth of predicted risk") + ylab("5 year risk (%)") + theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) + scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks=seq(0, 100, 10)) + scale_color_manual(values = c("observed" = 'skyblue1','pred' = 'darkblue'), labels=c("observed" = 'Observed','pred' = 'Predicted')) + theme(legend.position="bottom", legend.title=element_blank(), legend.text = element_text(size = 12, face = "bold")) + ggtitle("Male") 
}
p1 <- plot_male(dat.melt)

# predict survival probability and observed probability for female group
dat <- FUN.deciles_mean_pred(completed[completed$female==T,])
dat.melt <- melt(dat, id.vars = 'Deciles', variable.name = 'Group', value.name = 'Risk')
plot_female<-function(dat.melt){
  p2 <- ggplot(dat.melt, aes(x = Deciles, y = Risk*100, color = Group)) + geom_point(size = 4) +  xlab("Tenth of predicted risk") + ylab("5 year risk (%)") + theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) + scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks=seq(0, 100, 10)) + scale_color_manual(values = c("observed" = 'pink','pred' = 'darkred'), labels=c("observed" = 'Observed','pred' = 'Predicted')) + theme(legend.position="bottom", legend.title=element_blank(), legend.text = element_text(size = 12, face = "bold")) + ggtitle("Female") 
  p2 <- ggplot(dat.melt, aes(x = Deciles, y = Risk*100, color = Group)) + geom_point(size = 4) +  xlab("Tenth of predicted risk") + ylab("5 year risk (%)") + theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) + scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks=seq(0, 100, 10)) + scale_color_manual(values = c("observed" = 'pink','pred' = 'darkred'), labels=c("observed" = 'Observed','pred' = 'Predicted')) + theme(legend.position="bottom", legend.title=element_blank(), legend.text = element_text(size = 12, face = "bold")) + ggtitle("Female") 
}
p2 <- plot_female(dat.melt)
multiplot(p1, p2, cols = 1)