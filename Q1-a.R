library(jsonlite)
library(ggplot2)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(MCMCpack)

#Load Dataset 
data <- read.csv("C:/Users/Anny/Desktop/ASM/simpsons.csv")
desc(data)
#Pre-processing of Data
#select data
data <- data[c('Season','Rating')]
#filter data
q1data <- data %>% filter(Season == 2 | Season == 6)
head(q1data,10)
#   Season Rating
#1       2    8.3
#2       2    8.2
#3       2    8.3
#4       2    8.1
#5       2    7.4
#6       2    8.0
#7       2    7.7
#8       2    8.4
#9       2    8.1
#10      2    7.8
q1data$Season <- factor(q1data$Season)
summary(q1data)
#Season     Rating     
#2:22   Min.   :5.800  
#6:25   1st Qu.:8.000  
#       Median :8.200  
#       Mean   :8.221
#       3rd Qu.:8.600  
#        Max.   :9.200 
#q1data$Rating <- factor(q1data$Rating)
#summary(q1data)
#Season     Rating  
#2:22   8.3    : 7  
#6:25   8.1    : 6  
#       8.6    : 5  
#       8.2    : 4  
#       7.9    : 3  
#       8      : 3  
#       (Other):19 

season2 <- data %>% filter(Season == 2)
data$Season <- factor(data$Season)
summary(season2)
# Season       Rating     
#2      :22   Min.   :7.400  
#1      : 0   1st Qu.:7.825  
#3      : 0   Median :8.050  
#4      : 0   Mean   :8.045  
#5      : 0   3rd Qu.:8.300  
#6      : 0   Max.   :8.800  
season6 <- data %>% filter(Season == 6)
data$Season <- factor(data$Season)
summary(season6)
#plot
#ggplot(q1data) + geom_boxplot(aes(Season, Rating, fill = Season)) + geom_jitter(aes(Season, Rating, shape = Season))+ scale_fill_manual(values=c("green", "yellow"))
boxplot( Rating ~ Season, q1data, col = c("yellow", "green"))
#Comparing the means in a Bayesian model
compare_gibbs <- function(y, ind, mu0 = 50, tau0 = 1/400, del0 = 0, gamma0 = 1/400, a0 = 1, b0 = 50, maxiter = 5000)
{
  y1 <- y[ind == 2]
  y2 <- y[ind == 6]
  
  n1 <- length(y1) 
  n2 <- length(y2)
  
  ##### starting values
  mu <- (mean(y1) + mean(y2)) / 2
  del <- (mean(y1) - mean(y2)) / 2
  
  mat_store <- matrix(0, nrow = maxiter, ncol = 3)
  #####
  
  ##### Gibbs sampler
  an <- a0 + (n1 + n2)/2
  
  for(s in 1 : maxiter) 
  {
    
    ##update tau
    bn <- b0 + 0.5 * (sum((y1 - mu - del) ^ 2) + sum((y2 - mu + del) ^ 2))
    tau <- rgamma(1, an, bn)
    ##
    
    ##update mu
    taun <-  tau0 + tau * (n1 + n2)
    mun <- (tau0 * mu0 + tau * (sum(y1 - del) + sum(y2 + del))) / taun
    mu <- rnorm(1, mun, sqrt(1/taun))
    ##
    
    ##update del
    gamman <-  gamma0 + tau*(n1 + n2)
    deln <- ( del0 * gamma0 + tau * (sum(y1 - mu) - sum(y2 - mu))) / gamman
    del<-rnorm(1, deln, sqrt(1/gamman))
    ##
    
    ## store parameter values
    mat_store[s, ] <- c(mu, del, tau)
  }
  colnames(mat_store) <- c("mu", "del", "tau")
  return(mat_store)
}

library('MCMCpack')
fit_model <- compare_gibbs(q1data$Rating, as.factor(q1data$Season))
pdf('R.pdf')
plot(as.mcmc(fit_model))
dev.off()
raftery.diag(as.mcmc(fit_model))
apply(fit_model, 2, mean)
#        mu        del        tau 
#8.2224775 -0.1629151  0.4146883 
apply(fit_model, 2, sd)
#  mu        del        tau 
#0.23136137 0.23391632 0.08564099 
mean(1/sqrt(fit_model[, 3]))
#1.576394
sd(1/sqrt(fit_model[, 3]))
#0.1648255
apply(fit_model, 2, function(x) quantile(x, c(0.05, 0.95)))
#       mu        del       tau
#5%  7.828663 -0.5556151 0.2828140
#95% 8.604860  0.2078148 0.5652526
#Making posterior params
rating2_simu <- rnorm(5000, fit_model[, 1] + fit_model[, 2], sd = 1/sqrt(fit_model[, 3]))
rating6_simu <- rnorm(5000, fit_model[, 1] - fit_model[, 2], sd = 1/sqrt(fit_model[, 3]))
mean(rating2_simu < rating6_simu)
pdf('1.pdf')
ggplot(data.frame(rating2_simu, rating6_simu)) + geom_point(aes(rating2_simu, rating6_simu), alpha = 0.3) + geom_abline(slope = 1, intercept = 0)
dev.off()
pdf('2.pdf')
ggplot(data.frame(rating_simu_diff = rating6_simu - rating2_simu), aes(x=rating_simu_diff)) + stat_bin(aes(rating_simu_diff)) +geom_histogram(color="black", fill="white")
dev.off()