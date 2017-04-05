# A script for smoothing price using GDP deflator
# Assume that the % change is the same for the two quantities

rm(list=ls())

library(tidyverse)
library(reshape2)
library(parallel)


setwd("C:\\Users\\tomli\\Desktop\\DM")


GDP_deflator <- read.csv('GDP_deflator_2014_byYear.csv') %>%
      mutate(Time = 1:nrow(.)) 

HA_cost <- read.csv('HA cost.csv') %>%
      rename('2003' = X2003,
             '2013' = X2013) %>%
      melt(id = 'Service') %>%
      rename(Year = variable,
             cost = value)




# Calculate difference in slope
diff_perc_change <- function(x, y, weight) {
      sum1 <- sum(weight[-length(x)]^0.1 * abs(diff(x)/x[-length(x)] - diff(y)/y[-length(y)]))
      sum2 <- sum(weight[-1]^0.1 * abs(diff(x)/x[-1] - diff(y)/y[-1]))
      mean(c(sum1, sum2))
}

# Algorithm to impute missing value based on the slopes-------------------------
smooth <- function(service) {
      library(tidyverse)
      
      x <- GDP_deflator$GDP_deflator_2014
      y <- numeric(nrow(GDP_deflator))
      y[c(43,53)] <- HA_cost %>%
                        filter(Service == service) %>%
                        select(cost) %>%
                        unlist()
      # 
      # missing <- (1:length(y))[-c(43, 53)]
      # y[missing] <- NA
      # y[missing] <- mean(y, na.rm=TRUE)
      # 
      # not_missing <- (1:length(y))[-missing]
      # distance <- apply(sapply(missing, function(x) x - not_missing), 2, function(x) min(abs(x)))
      # missing <- missing[order(distance)]
      # 
      # weight <- 1/(1+apply(sapply(1:length(x), function(x) x - not_missing), 2, function(x) min(abs(x))))
      # 
      # for(k in 1:100) {
      #       for(i in missing) {
      #             iter <- 1
      #             while(iter < 30) {
      #                   old_coint <- diff_perc_change(x, y, weight)
      #                   y_i_old <- y[i]
      #                   y[i] <- rnorm(n = 1, mean = y[i], sd = mean(y)/10)
      #                   eps <- diff_perc_change(x, y, weight) - old_coint
      #                   if(eps > 0) {
      #                         y[i] <- y_i_old
      #                   }
      #                   iter <- iter + 1
      #             }      
      #       } 
      # }      
      
      d1 <- diff(x[c(43, 53)])
      d2 <- diff(y[c(43, 53)])
      y0 <- y[43]
      y <- d2/d1 * (x - x[43]) + y0
      
      output <- data.frame(y) %>%
            mutate(Year = GDP_deflator$Year,
                   Service = service) %>%
            rename(cost = y)
      
      return(output)
}

# no_cores <- detectCores() - 1
# cl <- makeCluster(no_cores)
# clusterExport(cl, list("smooth","HA_cost","GDP_deflator","diff_perc_change"),envir = .GlobalEnv)
# smooth_cost <- parLapply(cl, unique(HA_cost$Service), smooth)
# smooth_cost <- Reduce(rbind, smooth_cost)



# Figure to compare the two trends---------------------------------------------------------
figure <- function(service) {
      y <- smooth_cost %>%
                  filter(Service == service) %>%
                  select(cost) %>%
                  as.matrix()
      plot(GDP_deflator_2014 ~ Year, data = GDP_deflator, type= 'l', xaxt="n", ylab = 'GDP', xlab = 'Year', main = service)
      legend(x = 1970, y = 90, legend = c('GDP deflator',paste0('Estimated cost for ', service)), 
             lty = c(1, 0), 
             pch = c(-1, 1),
             col = c(1,2))
      par(new= TRUE)
      plot(y = y, x = GDP_deflator$Year, cex=1, yaxt="n",xaxt="n", ylab = "", xlab = '', col = 2)
      abline(v = c(2003, 2013))
      axis(side = 4, col = 2)
      mtext(side = 1, at = c(1960, 1980, 2003, 2013), text = c(1960, 1980, 2003, 2013), padj = 1)
}

# write.csv(smooth_cost, 'smooth_cost.csv')
