##############################################################
##                                                          ##
##              MODEL MISSPECIFICATION                      ##
##                                                          ##
##----------------------------------------------------------##  
##                                                          ##  
## This code produces the txt files necessary to reproduce  ##
## the results presented in the mainpaper Section 4.3 and   ##
## in the supplementary material regarding the simulation   ##
## studies on model misspecification                        ##
##                                                          ## 
##                                                          ## 
## Input                                                    ##    
##   - inpath: path containing the simulated datasets       ##
##   - outpath: destination path for results .txt files     ##
##                                                          ##
## Output                                                   ##        
##   - three txt files per simulation setting               ##
##                                                          ##
## Date: February 12, 2019                                  ##
##############################################################

rm(list = ls())
expit <- function(x) exp(x) / (1 + exp(x))
library(DTRreg)
# inpath <- 
# outpath <- 

#### ---- (cens1) censoring misspecified, Log-normal survival times, linear treatment-free ---- ####
## True censoring model: delta ~ X12
## Misspecified censoring model: delta ~ 1
## Sample sizes: 100, 300, 500, 1000, 5000, 10000
theta1 <- c(4.7, 1.5, -0.8, 0.1, 0.1)
for(n in c(100,300,500,1000,5000,10000)){
  load(paste(inpath, "/datasets", n, "_normal_linear.RData", sep = ""))
  cov <- width <- matrix(0, nrow = 1000, ncol = 10)
  comptime <- matrix(NA, nrow = 1000, ncol = 15)
  for(i in 1:1000){
    original.data <- datasets[[i]]
    # Adjusted
    t0adj <- proc.time()
    model1 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X12), cens.mod = list(delta ~ 1), var.estim = "asymptotic", asymp.opt = "adjusted", data = original.data)
    t1adj <- proc.time()
    # Naive
    t0nai <- proc.time()
    model2 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X12), cens.mod = list(delta ~ 1), var.estim = "asymptotic", asymp.opt = "naive", data = original.data)
    t1nai <- proc.time()
    # Standard boostrap
    t0stan <- proc.time()
    model3 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X12), cens.mod = list(delta ~ 1), var.estim = "bootstrap", boot.opt = "standard", B = 1000, data = original.data)
    t1stan <- proc.time()
    # Parametric bootstrap with empirical residual distribution
    t0para1 <- proc.time()
    model4 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X12), cens.mod = list(delta ~ 1), var.estim = "bootstrap", boot.opt = "empirical", B = 1000, data = original.data)
    t1para1 <- proc.time()
    # Parametric bootstrap with normal residual distribution
    t0para2 <- proc.time()
    model5 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X12), cens.mod = list(delta ~ 1), var.estim = "bootstrap", boot.opt = "normal", B = 1000, data = original.data)
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
  dest_cov <- paste(outpath, "/misspec_cens1_", n, "_cov.txt", sep = "")
  dest_width <- paste(outpath, "/misspec_cens1_", n, "_width.txt", sep = "")
  dest_comptime <- paste(outpath, "/misspec_cens1_", n, "_comptime.txt", sep = "")
  write.table(cov, dest_cov, col.names = F, row.names = F)
  write.table(width, dest_width, col.names = F, row.names = F)
  write.table(comptime, dest_comptime, col.names = F, row.names = F)
}

#### ---- (cens2) censoring misspecified, Log-normal survival times, nonlinear treatment-free ---- ####
## True censoring model: delta ~ X12
## Misspecified censoring model: delta ~ X1
## Sample sizes: 100, 300, 500, 1000, 5000, 10000
for(n in c(100,300,500,1000,5000,10000)){
  load(paste(inpath, "/datasets", n, "_normal_linear.RData", sep = ""))
  cov <- width <- matrix(0, nrow = 1000, ncol = 10)
  comptime <- matrix(NA, nrow = 1000, ncol = 15)
  for(i in 1:1000){
    original.data <- datasets[[i]]
    # Adjusted
    t0adj <- proc.time()
    model1 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X12), cens.mod = list(delta ~ X1), var.estim = "asymptotic", asymp.opt = "adjusted", data = original.data)
    t1adj <- proc.time()
    # Naive
    t0nai <- proc.time()
    model2 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X12), cens.mod = list(delta ~ X1), var.estim = "asymptotic", asymp.opt = "naive", data = original.data)
    t1nai <- proc.time()
    # Standard boostrap
    t0stan <- proc.time()
    model3 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X12), cens.mod = list(delta ~ X1), var.estim = "bootstrap", boot.opt = "standard", B = 1000, data = original.data)
    t1stan <- proc.time()
    # Parametric bootstrap with empirical residual distribution
    t0para1 <- proc.time()
    model4 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X12), cens.mod = list(delta ~ X1), var.estim = "bootstrap", boot.opt = "empirical", B = 1000, data = original.data)
    t1para1 <- proc.time()
    # Parametric bootstrap with normal residual distribution
    t0para2 <- proc.time()
    model5 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X12), cens.mod = list(delta ~ X1), var.estim = "bootstrap", boot.opt = "normal", B = 1000, data = original.data)
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
  dest_cov <- paste(outpath, "/misspec_cens2_", n, "_cov.txt", sep = "")
  dest_width <- paste(outpath, "/misspec_cens2_", n, "_width.txt", sep = "")
  dest_comptime <- paste(outpath, "/misspec_cens2_", n, "_comptime.txt", sep = "")
  write.table(cov, dest_cov, col.names = F, row.names = F)
  write.table(width, dest_width, col.names = F, row.names = F)
  write.table(comptime, dest_comptime, col.names = F, row.names = F)
}

#### ---- (trmt1) treatment misspecified, Log-normal survival times, linear treatment-free ---- ####
## True treatment model: A1 ~ X1
## Misspecified treatment model: A1 ~ 1
## Sample sizes: 100, 300, 500, 1000, 5000, 10000
for(n in c(100,300,500,1000,5000,10000)){
  load(paste(inpath, "/datasets", n, "_normal_linear.RData", sep = ""))
  cov <- width <- matrix(0, nrow = 1000, ncol = 10)
  comptime <- matrix(NA, nrow = 1000, ncol = 15)
  for(i in 1:1000){
    original.data <- datasets[[i]]
    # Adjusted
    t0adj <- proc.time()
    model1 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ 1), tf.mod = list( ~ X1 + X12), cens.mod = list(delta ~ X12), var.estim = "asymptotic", asymp.opt = "adjusted", data = original.data)
    t1adj <- proc.time()
    # Naive
    t0nai <- proc.time()
    model2 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ 1), tf.mod = list( ~ X1 + X12), cens.mod = list(delta ~ X12), var.estim = "asymptotic", asymp.opt = "naive", data = original.data)
    t1nai <- proc.time()
    # Standard boostrap
    t0stan <- proc.time()
    model3 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ 1), tf.mod = list( ~ X1 + X12), cens.mod = list(delta ~ X12), var.estim = "bootstrap", boot.opt = "standard", B = 1000, data = original.data)
    t1stan <- proc.time()
    # Parametric bootstrap with empirical residual distribution
    t0para1 <- proc.time()
    model4 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ 1), tf.mod = list( ~ X1 + X12), cens.mod = list(delta ~ X12), var.estim = "bootstrap", boot.opt = "empirical", B = 1000, data = original.data)
    t1para1 <- proc.time()
    # Parametric bootstrap with normal residual distribution
    t0para2 <- proc.time()
    model5 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ 1), tf.mod = list( ~ X1 + X12), cens.mod = list(delta ~ X12), var.estim = "bootstrap", boot.opt = "normal", B = 1000, data = original.data)
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
  dest_cov <- paste(outpath, "/misspec_trmt1_", n, "_cov.txt", sep = "")
  dest_width <- paste(outpath, "/misspec_trmt1_", n, "_width.txt", sep = "")
  dest_comptime <- paste(outpath, "/misspec_trmt1_", n, "_comptime.txt", sep = "")
  write.table(cov, dest_cov, col.names = F, row.names = F)
  write.table(width, dest_width, col.names = F, row.names = F)
  write.table(comptime, dest_comptime, col.names = F, row.names = F)
}

#### ---- (trmt2) treatment misspecified, Log-normal survival times, nonlinear treatment-free ---- ####
## True treatment model: A1 ~ X1
## Misspecified treatment model: A1 ~ X12
## Sample sizes: 100, 300, 500, 1000, 5000, 10000
for(n in c(100,300,500,1000,5000,10000)){
  load(paste(inpath, "/datasets", n, "_normal_linear.RData", sep = ""))
  cov <- width <- matrix(0, nrow = 1000, ncol = 10)
  comptime <- matrix(NA, nrow = 1000, ncol = 15)
  for(i in 1:1000){
    original.data <- datasets[[i]]
    # Adjusted
    t0adj <- proc.time()
    model1 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X12), tf.mod = list( ~ X1 + X12), cens.mod = list(delta ~ X12), var.estim = "asymptotic", asymp.opt = "adjusted", data = original.data)
    t1adj <- proc.time()
    # Naive
    t0nai <- proc.time()
    model2 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X12), tf.mod = list( ~ X1 + X12), cens.mod = list(delta ~ X12), var.estim = "asymptotic", asymp.opt = "naive", data = original.data)
    t1nai <- proc.time()
    # Standard boostrap
    t0stan <- proc.time()
    model3 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X12), tf.mod = list( ~ X1 + X12), cens.mod = list(delta ~ X12), var.estim = "bootstrap", boot.opt = "standard", B = 1000, data = original.data)
    t1stan <- proc.time()
    # Parametric bootstrap with empirical residual distribution
    t0para1 <- proc.time()
    model4 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X12), tf.mod = list( ~ X1 + X12), cens.mod = list(delta ~ X12), var.estim = "bootstrap", boot.opt = "empirical", B = 1000, data = original.data)
    t1para1 <- proc.time()
    # Parametric bootstrap with normal residual distribution
    t0para2 <- proc.time()
    model5 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X12), tf.mod = list( ~ X1 + X12), cens.mod = list(delta ~ X12), var.estim = "bootstrap", boot.opt = "normal", B = 1000, data = original.data)
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
  dest_cov <- paste(outpath, "/misspec_trmt2_", n, "_cov.txt", sep = "")
  dest_width <- paste(outpath, "/misspec_trmt2_", n, "_width.txt", sep = "")
  dest_comptime <- paste(outpath, "/misspec_trmt2_", n, "_comptime.txt", sep = "")
  write.table(cov, dest_cov, col.names = F, row.names = F)
  write.table(width, dest_width, col.names = F, row.names = F)
  write.table(comptime, dest_comptime, col.names = F, row.names = F)
}

#### ---- (tf1) treatment-free misspecified, Log-normal survival times, linear treatment-free ---- ####
## True treatment-free model: ~ X1 + X12
## Misspecified treatment-free model: ~ X1
## Sample sizes: 100, 300, 500, 1000, 5000, 10000
for(n in c(100,300,500,1000,5000,10000)){
  load(paste(inpath, "/datasets", n, "_normal_linear.RData", sep = ""))
  cov <- width <- matrix(0, nrow = 1000, ncol = 10)
  comptime <- matrix(NA, nrow = 1000, ncol = 15)
  for(i in 1:1000){
    original.data <- datasets[[i]]
    # Adjusted
    t0adj <- proc.time()
    model1 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1), cens.mod = list(delta ~ X12), var.estim = "asymptotic", asymp.opt = "adjusted", data = original.data)
    t1adj <- proc.time()
    # Naive
    t0nai <- proc.time()
    model2 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1), cens.mod = list(delta ~ X12), var.estim = "asymptotic", asymp.opt = "naive", data = original.data)
    t1nai <- proc.time()
    # Standard boostrap
    t0stan <- proc.time()
    model3 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1), cens.mod = list(delta ~ X12), var.estim = "bootstrap", boot.opt = "standard", B = 1000, data = original.data)
    t1stan <- proc.time()
    # Parametric bootstrap with empirical residual distribution
    t0para1 <- proc.time()
    model4 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1), cens.mod = list(delta ~ X12), var.estim = "bootstrap", boot.opt = "empirical", B = 1000, data = original.data)
    t1para1 <- proc.time()
    # Parametric bootstrap with normal residual distribution
    t0para2 <- proc.time()
    model5 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1), cens.mod = list(delta ~ X12), var.estim = "bootstrap", boot.opt = "normal", B = 1000, data = original.data)
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
  dest_cov <- paste(outpath, "/misspec_tf1_", n, "_cov.txt", sep = "")
  dest_width <- paste(outpath, "/misspec_tf1_", n, "_width.txt", sep = "")
  dest_comptime <- paste(outpath, "/misspec_tf1_", n, "_comptime.txt", sep = "")
  write.table(cov, dest_cov, col.names = F, row.names = F)
  write.table(width, dest_width, col.names = F, row.names = F)
  write.table(comptime, dest_comptime, col.names = F, row.names = F)
}

#### ---- (tf2) treatment-free misspecified, Log-normal survival times, nonlinear treatment-free ---- ####
## True treatment-free model: ~ X1 + X12 + X14
## Misspecified treatment-free model: ~ X1 + X12
## Sample sizes: 100, 300, 500, 1000, 5000, 10000
theta1 <- c(4.7, 3, -0.9, 0.05, 0.1, 0.1)
for(n in c(100,300,500,1000,5000,10000)){
  load(paste(inpath, "/datasets", n, "_normal_nonlinear.RData", sep = ""))
  cov <- width <- matrix(0, nrow = 1000, ncol = 10)
  comptime <- matrix(NA, nrow = 1000, ncol = 15)
  for(i in 1:1000){
    original.data <- datasets[[i]]
    # Adjusted
    t0adj <- proc.time()
    model1 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X12), cens.mod = list(delta ~ X12), var.estim = "asymptotic", asymp.opt = "adjusted", data = original.data)
    t1adj <- proc.time()
    # Naive
    t0nai <- proc.time()
    model2 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X12), cens.mod = list(delta ~ X12), var.estim = "asymptotic", asymp.opt = "naive", data = original.data)
    t1nai <- proc.time()
    # Standard boostrap
    t0stan <- proc.time()
    model3 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X12), cens.mod = list(delta ~ X12), var.estim = "bootstrap", boot.opt = "standard", B = 1000, data = original.data)
    t1stan <- proc.time()
    # Parametric bootstrap with empirical residual distribution
    t0para1 <- proc.time()
    model4 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X12), cens.mod = list(delta ~ X12), var.estim = "bootstrap", boot.opt = "empirical", B = 1000, data = original.data)
    t1para1 <- proc.time()
    # Parametric bootstrap with normal residual distribution
    t0para2 <- proc.time()
    model5 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X12), cens.mod = list(delta ~ X12), var.estim = "bootstrap", boot.opt = "normal", B = 1000, data = original.data)
    t1para2 <- proc.time()
    # Coverage
    if(confint(model1)$`stage1`[1,1] <= theta1[5] & confint(model1)$`stage1`[1,2] >= theta1[5]) cov[i,1] <- 1
    if(confint(model1)$`stage1`[2,1] <= theta1[6] & confint(model1)$`stage1`[2,2] >= theta1[6]) cov[i,2] <- 1
    if(confint(model2)$`stage1`[1,1] <= theta1[5] & confint(model2)$`stage1`[1,2] >= theta1[5]) cov[i,3] <- 1
    if(confint(model2)$`stage1`[2,1] <= theta1[6] & confint(model2)$`stage1`[2,2] >= theta1[6]) cov[i,4] <- 1
    if(confint(model3, type = "percentile")$`stage1`[1,1] <= theta1[5] & confint(model3, type = "percentile")$`stage1`[1,2] >= theta1[5]) cov[i,5] <- 1
    if(confint(model3, type = "percentile")$`stage1`[2,1] <= theta1[6] & confint(model3, type = "percentile")$`stage1`[2,2] >= theta1[6]) cov[i,6] <- 1  
    if(confint(model4, type = "percentile")$`stage1`[1,1] <= theta1[5] & confint(model4, type = "percentile")$`stage1`[1,2] >= theta1[5]) cov[i,7] <- 1
    if(confint(model4, type = "percentile")$`stage1`[2,1] <= theta1[6] & confint(model4, type = "percentile")$`stage1`[2,2] >= theta1[6]) cov[i,8] <- 1  
    if(confint(model5, type = "percentile")$`stage1`[1,1] <= theta1[5] & confint(model5, type = "percentile")$`stage1`[1,2] >= theta1[5]) cov[i,9] <- 1
    if(confint(model5, type = "percentile")$`stage1`[2,1] <= theta1[6] & confint(model5, type = "percentile")$`stage1`[2,2] >= theta1[6]) cov[i,10] <- 1 
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
  dest_cov <- paste(outpath, "/misspec_tf2_", n, "_cov.txt", sep = "")
  dest_width <- paste(outpath, "/misspec_tf2_", n, "_width.txt", sep = "")
  dest_comptime <- paste(outpath, "/misspec_tf2_", n, "_comptime.txt", sep = "")
  write.table(cov, dest_cov, col.names = F, row.names = F)
  write.table(width, dest_width, col.names = F, row.names = F)
  write.table(comptime, dest_comptime, col.names = F, row.names = F)
}