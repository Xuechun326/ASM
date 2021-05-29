library(ggplot2)
library(MCMCpack)
library(RColorBrewer)
#Load Dataset 
data <- read.csv("C:/Users/Anny/Desktop/ASM/simpsons.csv")
#Pre-processing of Data
#select data
data <- data[c('Season','Rating')]
dim(data)
head(data)
data$Season <- factor(data$Season)
nlevels(data$Season)
#pdf('R1.pdf')

#ggplot(data) + geom_boxplot(aes(x = reorder(Season, Rating, median), Rating, 
#                              fill = reorder(Season, Rating, median)), show.legend=FALSE)

#dev.off()
#pdf('R2.pdf')
#ggplot(data, aes(x = reorder(Season, Season, length))) + stat_count()
#dev.off()
#pdf('R3.pdf')
#ggplot(data, aes(Rating)) + stat_bin()
#dev.off()
pdf('R5.pdf')
ggplot(data.frame(size = tapply(data$Rating, data$Season, length), 
                  mean_score = tapply(data$Rating, data$Season, mean)), 
       aes(size, mean_score)) + geom_point()
dev.off()
# build model
compare_m_gibbs <- function(y, ind, maxiter = 5000)
{
  
  ### weakly informative priors
  a0 <- 1/2 ; b0 <- 50 ## tau_w hyperparameters
  eta0 <-1/2 ; t0 <- 50 ## tau_b hyperparameters
  mu0<-50 ; gamma0 <- 1/25
  ###
  
  ### starting values
  m <- nlevels(ind)
  ybar <- theta <- tapply(y, ind, mean)
  tau_w <- mean(1 / tapply(y, ind, var)) ##within group precision
  mu <- mean(theta)
  tau_b <-var(theta) ##between group precision
  n_m <- tapply(y, ind, length)
  an <- a0 + sum(n_m)/2
  ###
  
  ### setup MCMC
  theta_mat <- matrix(0, nrow=maxiter, ncol=m)
  mat_store <- matrix(0, nrow=maxiter, ncol=3)
  ###
  
  ### MCMC algorithm
  for(s in 1:maxiter) 
  {
    
    # sample new values of the thetas
    for(j in 1:m) 
    {
      taun <- n_m[j] * tau_w + tau_b
      thetan <- (ybar[j] * n_m[j] * tau_w + mu * tau_b) / taun
      theta[j]<-rnorm(1, thetan, 1/sqrt(taun))
    }
    
    #sample new value of tau_w
    ss <- 0
    for(j in 1:m){
      ss <- ss + sum((y[ind == j] - theta[j])^2)
    }
    bn <- b0 + ss/2
    tau_w <- rgamma(1, an, bn)
    
    #sample a new value of mu
    gammam <- m * tau_b + gamma0
    mum <- (mean(theta) * m * tau_b + mu0 * gamma0) / gammam
    mu <- rnorm(1, mum, 1/ sqrt(gammam)) 
    
    # sample a new value of tau_b
    etam <- eta0 + m/2
    tm <- t0 + sum((theta - mu)^2) / 2
    tau_b <- rgamma(1, etam, tm)
    
    #store results
    theta_mat[s,] <- theta
    mat_store[s, ] <- c(mu, tau_w, tau_b)
  }
  colnames(mat_store) <- c("mu", "tau_w", "tau_b")
  return(list(params = mat_store, theta = theta_mat))
}
fit2 <- compare_m_gibbs(data$Rating, data$Season)
# reformat samples for ggplot
theta_df <- data.frame(Rating = as.numeric(fit2$theta), 
                       Season = rep(1:ncol(fit2$theta), each = nrow(fit2$theta))) 

theta_med <- apply(theta_df, 2, mean) ## get basic posterior summary
sort(theta_med, decreasing = TRUE) ## 
pdf('R4.pdf')
ggplot(theta_df) + geom_boxplot(aes(x = reorder(Season, Rating, median), Rating, 
                                fill = reorder(Season, Rating, median)), show.legend=FALSE)

dev.off()

