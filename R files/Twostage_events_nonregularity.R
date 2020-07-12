##############################################################
##                                                          ##
##       TWO-STAGE DTR WITH EVENTS IN BOTH STAGES           ##
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
## Date: Jul 12, 2020                                      ##
##############################################################
rm(list = ls())
expit <- function(x) exp(x) / (1 + exp(x))
# inpath <- 
# outpath <- 

#### ------------- Regular/nonregular scenarios ------------- ####
## 8 regular/near nonregular/nonregular scenarios, with events in either stages
## Sample sizes: 300, 1000

theta1 <- c(6.3, 0.5, -0.01, 0.1, 0.1)
theta2 <- c(4, 1.1, 0.01, -0.2, 0.1, NA, NA)
b2 <- matrix(NA, nrow = 8, ncol = 2)
b2[1,] <- c(0, 0); b2[2,] <- c(0.01, 0); b2[3,] <- c(-0.5,0.5); b2[4,] <- c(-0.5, 0.49)
b2[5,] <- c(-0.2, 0.2); b2[6,] <- c(-0.2, 0.19); b2[7,] <- c(-0.9); b2[8,] <- c(0.2, -0.7)
p2 <- rep(0.5, 9)
p2[1] <- p2[2] <- 0.3; p2[5] <- 0.25

for(sc in 1:8){
  theta2[c(6,7)] <- b2[sc,]
  for(n in c(300,1000)){
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
        delta <- rbinom(n, size = 1, expit(3*X12 + 0.1))
        eta2 <- rbinom(n, 1, prob = 0.8)
        delta2 <- delta[eta2 == 1]
        
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
    save.image(paste(outpath, "/datasets", n, "_twoev_sc", sc,".RData", sep = ""))
  }
}

#### ---- 8 scenarios, 2 stage DTR with Log-normal survival times, linear treatment-free ---- ####
## Corresponding scenario ID in Table 1 main paper
## Sample sizes: 300, 1000
theta1 <- c(6.3, 0.5, -0.01, 0.1, 0.1)

for(sc in c(1,2,7)){
  for(n in c(300,1000)){
    load(paste(inpath, "/datasets", n, "_twoev_sc", sc, ".RData", sep = ""))
    cov <- width <- matrix(0, nrow = 1000, ncol = 10)
    comptime <- matrix(NA, nrow = 1000, ncol = 15)
    for(i in 1:1000){
      if(i%%100==0){
        print(paste("Scenario ",sc,", n=",n,", i=",i,sep=""))
      }
      original.data <- datasets[[i]]
      # Adjusted
      t0adj <- proc.time()
      model1 <- DWSurv(time = list(~Y1, ~Y2), blip.mod = list(~ X11, ~ X22), treat.mod = list(A1 ~ X11, A2 ~ X21), tf.mod = list( ~ X11 + X12, ~ X21 + X22 + X23 + X11), cens.mod = list(delta ~ X12, delta ~ X12), var.estim = "asymptotic", asymp.opt = "adjusted", data = original.data)
      t1adj <- proc.time()
      # Naive
      t0nai <- proc.time()
      model2 <- DWSurv(time = list(~Y1, ~Y2), blip.mod = list(~ X11, ~ X22), treat.mod = list(A1 ~ X11, A2 ~ X21), tf.mod = list( ~ X11 + X12, ~ X21 + X22 + X23 + X11), cens.mod = list(delta ~ X12, delta ~ X12), var.estim = "asymptotic", asymp.opt = "naive", data = original.data)
      t1nai <- proc.time()
      # Standard bootstrap
      t0stan <- proc.time()
      model3 <- DWSurv(time = list(~Y1, ~Y2), blip.mod = list(~ X11, ~ X22), treat.mod = list(A1 ~ X11, A2 ~ X21), tf.mod = list( ~ X11 + X12, ~ X21 + X22 + X23 + X11), cens.mod = list(delta ~ X12, delta ~ X12), var.estim = "bootstrap", boot.opt = "standard", B = 1000, data = original.data)
      t1stan <- proc.time()
      # Parametric bootstrap with empirical residual distribution
      t0para1 <- proc.time()
      model4 <- DWSurv(time = list(~Y1, ~Y2), blip.mod = list(~ X11, ~ X22), treat.mod = list(A1 ~ X11, A2 ~ X21), tf.mod = list( ~ X11 + X12, ~ X21 + X22 + X23 + X11), cens.mod = list(delta ~ X12, delta ~ X12), var.estim = "bootstrap", boot.opt = "empirical", B = 1000, data = original.data)
      t1para1 <- proc.time()
      # Parametric bootstrap with normal residual distribution
      t0para2 <- proc.time()
      model5 <- DWSurv(time = list(~Y1, ~Y2), blip.mod = list(~ X11, ~ X22), treat.mod = list(A1 ~ X11, A2 ~ X21), tf.mod = list( ~ X11 + X12, ~ X21 + X22 + X23 + X11), cens.mod = list(delta ~ X12, delta ~ X12), var.estim = "bootstrap", boot.opt = "normal", B = 1000, data = original.data)
      t1para2 <- proc.time()
      # Coverage
      if(confint(model1)$`stage1`[1,1] <= theta1[4] & confint(model1)$`stage1`[1,2] >= theta1[4]) cov[i,1] <- 1
      if(confint(model1)$`stage1`[2,1] <= theta1[5] & confint(model1)$`stage1`[2,2] >= theta1[5]) cov[i,2] <- 1
      if(confint(model2)$`stage1`[1,1] <= theta1[4] & confint(model2)$`stage1`[1,2] >= theta1[4]) cov[i,3] <- 1
      if(confint(model2)$`stage1`[2,1] <= theta1[5] & confint(model2)$`stage1`[2,2] >= theta1[5]) cov[i,4] <- 1
      if(confint(model3, type = "percentile")$`stage1`[1,1] <= theta1[4] & confint(model3, type = "percentile")$`stage1`[1,2] >= theta1[4]) cov[i,5] <- 1
      if(confint(model3, type = "percentile")$`stage1`[2,1] <= theta1[5] & confint(model3, type = "percentile")$`stage1`[2,2] >= theta1[5]) cov[i,6] <- 1  
      if(confint(model4, type = "percentile")$`stage1`[1,1] <= theta1[4] & confint(model4, type = "percentile")$`stage1`[1,2] >= theta1[4]) cov[i,7] <- 1
      if(confint(model4, type = "percentile")$`stage1`[2,1] <= theta1[5] & confint(model4, type = "percentile")$`stage1`[2,2] >= theta1[5]) cov[i,8] <- 1  
      if(confint(model5, type = "percentile")$`stage1`[1,1] <= theta1[4] & confint(model5, type = "percentile")$`stage1`[1,2] >= theta1[4]) cov[i,9] <- 1
      if(confint(model5, type = "percentile")$`stage1`[2,1] <= theta1[5] & confint(model5, type = "percentile")$`stage1`[2,2] >= theta1[5]) cov[i,10] <- 1 
      # Width
      width[i,1] <- confint(model1)$`stage1`[1,2] - confint(model1)$`stage1`[1,1]
      width[i,2] <- confint(model1)$`stage1`[2,2] - confint(model1)$`stage1`[2,1]
      width[i,3] <- confint(model2)$`stage1`[1,2] - confint(model2)$`stage1`[1,1]
      width[i,4] <- confint(model2)$`stage1`[2,2] - confint(model2)$`stage1`[2,1]
      width[i,5] <- confint(model3, type = "percentile")$`stage1`[1,2] - confint(model3, type = "percentile")$`stage1`[1,1]
      width[i,6] <- confint(model3, type = "percentile")$`stage1`[2,2] - confint(model3, type = "percentile")$`stage1`[2,1]
      width[i,7] <- confint(model4, type = "percentile")$`stage1`[1,2] - confint(model4, type = "percentile")$`stage1`[1,1]
      width[i,8] <- confint(model4, type = "percentile")$`stage1`[2,2] - confint(model4, type = "percentile")$`stage1`[2,1]
      width[i,9] <- confint(model5, type = "percentile")$`stage1`[1,2] - confint(model5, type = "percentile")$`stage1`[1,1]
      width[i,10] <- confint(model5, type = "percentile")$`stage1`[2,2] - confint(model5, type = "percentile")$`stage1`[2,1]
      # Computational time
      comptime[i,1:3] <- as.numeric(t1adj[1:3]) - as.numeric(t0adj[1:3])
      comptime[i,4:6] <- as.numeric(t1nai[1:3]) - as.numeric(t0nai[1:3])
      comptime[i,7:9] <- as.numeric(t1stan[1:3]) - as.numeric(t0stan[1:3])
      comptime[i,10:12] <- as.numeric(t1para1[1:3]) - as.numeric(t0para1[1:3])
      comptime[i,13:15] <- as.numeric(t1para2[1:3]) - as.numeric(t0para2[1:3])
    } 
    dest_cov <- paste(outpath, "/nonreg_twoev_sc", sc, "_", n, "_cov.txt", sep = "")
    dest_width <- paste(outpath, "/nonreg_twoev_sc", sc, "_", n, "_width.txt", sep = "")
    dest_comptime <- paste(outpath, "/nonreg_twoev_sc", sc, "_", n, "_comptime.txt", sep = "")
    write.table(cov, dest_cov, col.names = F, row.names = F)
    write.table(width, dest_width, col.names = F, row.names = F)
    write.table(comptime, dest_comptime, col.names = F, row.names = F)
  }
}


