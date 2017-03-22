# A script for smoothing price using GDP deflator
# Assume that the % change is the same for the two quantities

rm(list=ls())

library(tidyverse)

setwd("C:\\Users\\tomli\\Desktop\\DM")

GDP_inflator <- read.csv('GDP_deflator_2014.csv') %>%
      mutate(Year = rep(1973:2016, each = 4),
             Quarter = rep(1:4, length(unique(Year))),
             perc_change = as.numeric(paste0(perc_change)),
             Time = 1:nrow(.)) 

plot(GDP_inflator_2014 ~ Time, data = GDP_inflator, type = 'l')

# Cointegrated time series----------------------------
diff_perc_change <- function(x, y) {
      sum1 <- sum(weight[-length(x)]^2 * abs(diff(x)/x[-length(x)] - diff(y)/y[-length(y)]))
      sum2 <- sum(weight[-1]^2 * abs(diff(x)/x[-1] - diff(y)/y[-1]))
      mean(c(sum1, sum2))
}

x <- GDP_inflator$GDP_inflator_2014
y <- x - rnorm(n = length(x), mean = 10, sd = 0.1)
plot(x, type= 'l', ylim = c(0, 120))
points(y, cex=0.5)

# missing <- sample(1:length(y), length(y)-20)
missing <- (1:length(y))[-c(1,50,100,150)]
y[missing] <- NA
y[missing] <- mean(y, na.rm=TRUE)

not_missing <- (1:length(y))[-missing]
distance <- apply(sapply(missing, function(x) x - not_missing), 2, function(x) min(abs(x)))
missing <- missing[order(distance)]

weight <- 1/(1+apply(sapply(1:length(x), function(x) x - not_missing), 2, function(x) min(abs(x))))

for(k in 1:30) {
      for(i in missing) {
            iter <- 1
            while(iter < 10) {
                  old_coint <- diff_perc_change(x, y)
                  y_i_old <- y[i]
                  y[i] <- rnorm(n = 1, mean = y[i], sd = 5)
                  eps <- diff_perc_change(x, y) - old_coint
                  if(eps > 0) {
                        y[i] <- y_i_old
                  }
                  iter <- iter + 1
            }      
      } 
      plot(x, type= 'l', ylim = c(0, 120))
      points(y, cex=0.5)
}      
