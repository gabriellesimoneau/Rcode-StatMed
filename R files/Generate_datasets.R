##############################################################
##                                                          ##
##       GENERATE DATASETS FOR ALL SIMULATIONS              ##
##                                                          ##
##----------------------------------------------------------##  
##                                                          ## 
## This code produces the simulated datasets following the  ##
## data generating mechanisms presented in the main paper   ##
## Section 4.1.                                             ##
##                                                          ## 
## Input                                                    ##    
##   - outpath: destination path for the simulated datasets ##
##                                                          ##
## Output                                                   ##
##   - RData files, each with 1000 simulated datasets       ##
##                                                          ## 
## Date: March 5, 2019                                      ##
##############################################################
rm(list = ls())
expit <- function(x) exp(x) / (1 + exp(x))
# outpath <- 

#### --------- Unknown true error distribution --------- ####
## Log-normal or Weibull survival times
## Linear or nonlinear treatment-free model
## Sample sizes: 100, 300, 500, 1000, 5000, 10000

## Log-normal survival times, linear treatment-free
theta1 <- c(4.7, 1.5, -0.8, 0.1, 0.1)
for(n in c(100,300,500,1000,5000,10000)){
  datasets <- vector("list", length = 1000)
  for(i in 1:1000){
    X1 <- runif(n, 0.1, 1.29) 
    X12 <- rbinom(n, 1, 0.4)
    A1 <- rbinom(n, 1, expit(2*X1 - 1))
    delta <- rbinom(n, 1, expit(3*X12 + 0.1))
    logT <- theta1[1] + theta1[2]*X1[delta == 1] + theta1[3]*X12[delta == 1] + theta1[4]*A1[delta == 1] + theta1[5]*A1[delta == 1]*X1[delta == 1] + rnorm(sum(delta), sd = 0.3)
    C <- rexp(n - sum(delta), rate = 1/300)
    Y <- rep(NA, n)
    Y[delta == 1] <- exp(logT)
    Y[delta == 0] <- C
    datasets[[i]] <- data.frame(X1, X12, A1, delta, Y)
  }
  save.image(paste(outpath, "/datasets", n, "_normal_linear.RData", sep = ""))
}

## Weibull survival times, linear treatment-free
for(n in c(100,300,500,1000,5000,10000)){
  datasets <- vector("list", length = 1000)
  for(i in 1:1000){
    X1 <- runif(n, 0.1, 1.29) 
    X12 <- rbinom(n, 1, 0.4)
    A1 <- rbinom(n, 1, expit(2*X1 - 1))
    delta <- rbinom(n, 1, expit(3*X12 + 0.1))
    logT <- theta1[1] + theta1[2]*X1[delta == 1] + theta1[3]*X12[delta == 1] + theta1[4]*A1[delta == 1] + theta1[5]*A1[delta == 1]*X1[delta == 1] + log(rweibull(sum(delta), shape = 4, scale = 1))
    
    C <- rexp(n - sum(delta), rate = 1/300)
    Y <- rep(NA, n)
    Y[delta == 1] <- exp(logT)
    Y[delta == 0] <- C
    
    datasets[[i]] <- data.frame(X1, X12, A1, delta, Y)
  }
  save.image(paste(outpath, "/datasets", n, "_weibull_linear.RData", sep = ""))
}

## Log-normal survival times, nonlinear treatment-free
theta1 <- c(4.7, 3, -0.9, 0.05, 0.1, 0.1)
for(n in c(100,300,500,1000,5000,10000)){
  datasets <- vector("list", length = 1000)
  for(i in 1:1000){
    X1 <- runif(n, 0.1, 1.29) 
    X14 <- X1^4
    X12 <- rbinom(n, 1, 0.4)
    A1 <- rbinom(n, 1, expit(2*X1 - 1))
    delta <- rbinom(n, 1, expit(3*X12 + 0.1))
    logT <- theta1[1] + theta1[2]*X1[delta == 1] + theta1[3]*X14[delta == 1] + theta1[4]*X12[delta == 1] + theta1[5]*A1[delta == 1] + theta1[6]*A1[delta == 1]*X1[delta == 1] + rnorm(sum(delta), sd = 0.3)
    
    C <- rexp(n - sum(delta), rate = 1/300)
    Y <- rep(NA, n)
    Y[delta == 1] <- exp(logT)
    Y[delta == 0] <- C
    
    datasets[[i]] <- data.frame(X1, X14, X12, A1, delta, Y)
  }
  save.image(paste(outpath, "/datasets", n, "_normal_nonlinear.RData", sep = ""))
}

## Weibull survival times, nonlinear treatment-free
for(n in c(100,300,500,1000,5000,10000)){
  datasets <- vector("list", length = 1000)
  for(i in 1:1000){
    X1 <- runif(n, 0.1, 1.29) 
    X14 <- X1^4
    X12 <- rbinom(n, 1, 0.4)
    A1 <- rbinom(n, 1, expit(2*X1 - 1))
    delta <- rbinom(n, 1, expit(3*X12 + 0.1))
    logT <- theta1[1] + theta1[2]*X1[delta == 1] + theta1[3]*X14[delta == 1] + theta1[4]*X12[delta == 1] + theta1[5]*A1[delta == 1] + theta1[6]*A1[delta == 1]*X1[delta == 1] + log(rweibull(sum(delta), shape = 4, scale = 1))
    
    C <- rexp(n - sum(delta), rate = 1/300)
    Y <- rep(NA, n)
    Y[delta == 1] <- exp(logT)
    Y[delta == 0] <- C
    
    datasets[[i]] <- data.frame(X1, X14, X12, A1, delta, Y)
  }
  save.image(paste(outpath, "/datasets", n, "_weibull_nonlinear.RData", sep = ""))
}



#### ------------- Regular/nonregular scenarios ------------- ####
## 8 regular/near nonregular/nonregular scenarios
## Sample sizes: 300, 500, 1000, 10000

theta1 <- c(6.3, 0.5, -0.01, 0.1, 0.1)
theta2 <- c(4, 1.1, 0.01, -0.2, 0.1, NA, NA)
b2 <- matrix(NA, nrow = 8, ncol = 2)
b2[1,] <- c(0, 0); b2[2,] <- c(0.01, 0); b2[3,] <- c(-0.5,0.5); b2[4,] <- c(-0.5, 0.49)
b2[5,] <- c(-0.2, 0.2); b2[6,] <- c(-0.2, 0.19); b2[7,] <- c(-0.9); b2[8,] <- c(0.2, -0.7)
p2 <- rep(0.5, 9)
p2[1] <- p2[2] <- 0.3; p2[5] <- 0.25

for(sc in 1:8){
  theta2[c(6,7)] <- b2[sc,]
  for(n in c(300,500,1000,10000)){
    datasets <- vector("list", length = 1000)
    for(i in 1:1000){
      t1ok <- 1
      while(t1ok){
        X11 <- runif(n, 0.1, 1.29) 
        X12 <- rbinom(n, 1, 0.4)
        A1 <- rbinom(n, 1, expit(2*X11 - 1))
        X21 <- runif(n, 0.9, 2)
        X22 <- rbinom(n, 1, p2[sc]) 
        X23 <- X21^3
        A2 <- rbinom(n, size = 1, prob = expit(-2*X21 + 2.8))
        delta <- rbinom(n, 1, expit(3*X12 + 0.1))
        eta2 <- rep(1, n) # everybody enters the second stage
        
        logT2 <- theta2[1] + theta2[2]*X21[eta2 == 1 & delta == 1] + theta2[3]*X22[eta2 == 1 & delta == 1] + theta2[4]*X23[eta2 == 1 & delta == 1] + theta2[5]*X11[eta2 == 1 & delta == 1] + theta2[6]*A2[eta2 == 1 & delta == 1] + theta2[7]*A2[eta2 == 1 & delta == 1]*X22[eta2 == 1 & delta == 1] + rnorm(sum(eta2*delta), sd = 0.3)
        trueA2opt <- ifelse(theta2[6] + theta2[7]*X22[eta2 == 1 & delta == 1] > 0, 1, 0)
        logT2opt <- logT2 + (trueA2opt - A2[eta2 == 1 & delta == 1])*(theta2[6] + theta2[7]*X22[eta2 == 1 & delta == 1])
        
        logT <- theta1[1] + theta1[2]*X11 + theta1[3]*X12 + theta1[4]*A1 + theta1[5]*A1*X11 + rnorm(n,sd = 0.3) 
        T1 <- exp(logT[eta2 == 1 & delta == 1]) - exp(logT2opt)
        if(min(T1) > 0) t1ok <- 0
      }
      C <- rexp(n - sum(delta), rate = 1/300)
      Y2 <- rep(NA, n)
      Y1 <- rep(NA, n)
      eta2d0 <- eta2[delta == 0]
      C1 <- rep(NA, length(C))
      C2 <- rep(NA, length(C))
      for(j in 1:length(C))
      {
        if(eta2d0[j] == 0){ # they should all be censored in second stage
          C1[j] <-  C[j]
          C2[j] <- 0
        }else{
          C1[j] <- runif(1, 0, C[j])
          C2[j] <- C[j] - C1[j]
        }
      }
      Y2[delta == 0] <- C2
      Y1[delta == 0] <- C1
      Y1[delta == 1 & eta2 == 1] <- T1
      Y1[delta == 1 & eta2 == 0] <- exp(logT[delta == 1 & eta2 == 0])
      Y2[delta == 1 & eta2 == 0] <- 0
      Y2[delta == 1 & eta2 == 1] <- exp(logT2)
      
      datasets[[i]] <- data.frame(X11, X12, A1, X21, X22, X23, A2, delta, Y1, Y2, eta2)
    }
    save.image(paste(outpath, "/datasets", n, "_sc", sc,".RData", sep = ""))
  }
}


