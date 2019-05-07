##############################################################
##                                                          ##
##            UNKNOWN ERROR DISTRIBUTION                    ##
##                                                          ##
##----------------------------------------------------------##  
##                                                          ##  
## This code produces the txt files necessary to reproduce  ##
## the results presented in the mainpaper Section 4.2 and   ##
## in the supplementary material regarding the simulation   ##
## studies on unknown survival time distribution            ##
## (Log-normal or Weibull).                                 ##
##                                                          ## 
##                                                          ## 
## Input                                                    ##    
##   - inpath: path of the simulated datasets               ##
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

#### ---- Log-normal survival times, linear treatment-free ---- ####
theta1 <- c(4.7, 1.5, -0.8, 0.1, 0.1)
for(n in c(100,300,500,1000,5000,10000)){
  load(paste(inpath, "/datasets", n, "_normal_linear.RData", sep = ""))
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
  dest_cov <- paste(outpath, "/error_linearnormal_", n, "_cov.txt", sep = "")
  dest_width <- paste(outpath, "/error_linearnormal_", n, "_width.txt", sep = "")
  dest_comptime <- paste(outpath, "/error_linearnormal_", n, "_comptime.txt", sep = "")
  write.table(cov, dest_cov, col.names = F, row.names = F)
  write.table(width, dest_width, col.names = F, row.names = F)
  write.table(comptime, dest_comptime, col.names = F, row.names = F)
}

#### ---- Weibull survival times, linear treatment-free ---- ####
for(n in c(100,300,500,1000,5000,10000)){
  load(paste(inpath, "/datasets", n, "_weibull_linear.RData", sep = ""))
  cov <- width <- matrix(0, nrow = 1000, ncol = 10)
  comptime <- matrix(NA, nrow = 1000, ncol = 15)
  for(i in 1:1000){
    original.data <- datasets[[i]]
    # Asymptotic
    t0adj <- proc.time()
    model1 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X12), cens.mod = list(delta ~ X12), var.estim = "asymptotic", asymp.opt = "adjusted", data = original.data)
    t1adj <- proc.time()
    # Finite
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
  dest_cov <- paste(outpath, "/error_linearweibull_", n, "_cov.txt", sep = "")
  dest_width <- paste(outpath, "/error_linearweibull_", n, "_width.txt", sep = "")
  dest_comptime <- paste(outpath, "/error_linearweibull_", n, "_comptime.txt", sep = "")
  write.table(cov, dest_cov, col.names = F, row.names = F)
  write.table(width, dest_width, col.names = F, row.names = F)
  write.table(comptime, dest_comptime, col.names = F, row.names = F)
}

#### ---- Log-normal survival times, nonlinear treatment-free ---- ####
theta1 <- c(4.7, 3, -0.9, 0.05, 0.1, 0.1)
for(n in c(100,300,500,1000,5000,10000)){
  load(paste(inpath, "/datasets", n, "_normal_nonlinear.RData", sep = ""))
  cov <- width <- matrix(0, nrow = 1000, ncol = 10)
  comptime <- matrix(NA, nrow = 1000, ncol = 15)
  for(i in 1:1000){
    original.data <- datasets[[i]]
    # Asymptotic
    t0adj <- proc.time()
    model1 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X14 + X12), cens.mod = list(delta ~ X12), var.estim = "asymptotic", asymp.opt = "adjusted", data = original.data)
    t1adj <- proc.time()
    # Finite
    t0nai <- proc.time()
    model2 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X14 + X12), cens.mod = list(delta ~ X12), var.estim = "asymptotic", asymp.opt = "naive", data = original.data)
    t1nai <- proc.time()
    # Standard boostrap
    t0stan <- proc.time()
    model3 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X14 + X12), cens.mod = list(delta ~ X12), var.estim = "bootstrap", boot.opt = "standard", B = 1000, data = original.data)
    t1stan <- proc.time()
    # Parametric bootstrap with empirical residual distribution
    t0para1 <- proc.time()
    model4 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X14 + X12), cens.mod = list(delta ~ X12), var.estim = "bootstrap", boot.opt = "empirical", B = 1000, data = original.data)
    t1para1 <- proc.time()
    # Parametric bootstrap with normal residual distribution
    t0para2 <- proc.time()
    model5 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X14 + X12), cens.mod = list(delta ~ X12), var.estim = "bootstrap", boot.opt = "normal", B = 1000, data = original.data)
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
  dest_cov <- paste(outpath, "/error_nonlinearnormal_", n, "_cov.txt", sep = "")
  dest_width <- paste(outpath, "/error_nonlinearnormal_", n, "_width.txt", sep = "")
  dest_comptime <- paste(outpath, "/error_nonlinearnormal_", n, "_comptime.txt", sep = "")
  write.table(cov, dest_cov, col.names = F, row.names = F)
  write.table(width, dest_width, col.names = F, row.names = F)
  write.table(comptime, dest_comptime, col.names = F, row.names = F)
}

#### ---- Weibull survival times, nonlinear treatment-free ---- ####
for(n in c(100,300,500,1000,5000,10000)){
  load(paste(inpath, "/datasets", n, "_weibull_nonlinear.RData", sep = ""))
  cov <- width <- matrix(0, nrow = 1000, ncol = 10)
  comptime <- matrix(NA, nrow = 1000, ncol = 15)
  for(i in 1:1000){
    original.data <- datasets[[i]]
    # Asymptotic
    t0adj <- proc.time()
    model1 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X14 + X12), cens.mod = list(delta ~ X12), var.estim = "asymptotic", asymp.opt = "adjusted", data = original.data)
    t1adj <- proc.time()
    # Finite
    t0nai <- proc.time()
    model2 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X14 + X12), cens.mod = list(delta ~ X12), var.estim = "asymptotic", asymp.opt = "naive", data = original.data)
    t1nai <- proc.time()
    # Standard boostrap
    t0stan <- proc.time()
    model3 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X14 + X12), cens.mod = list(delta ~ X12), var.estim = "bootstrap", boot.opt = "standard", B = 1000, data = original.data)
    t1stan <- proc.time()
    # Parametric bootstrap with empirical residual distribution
    t0para1 <- proc.time()
    model4 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X14 + X12), cens.mod = list(delta ~ X12), var.estim = "bootstrap", boot.opt = "empirical", B = 1000, data = original.data)
    t1para1 <- proc.time()
    # Parametric bootstrap with normal residual distribution
    t0para2 <- proc.time()
    model5 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X14 + X12), cens.mod = list(delta ~ X12), var.estim = "bootstrap", boot.opt = "normal", B = 1000, data = original.data)
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
  dest_cov <- paste(outpath, "/error_nonlinearweibull_", n, "_cov.txt", sep = "")
  dest_width <- paste(outpath, "/error_nonlinearweibull_", n, "_width.txt", sep = "")
  dest_comptime <- paste(outpath, "/error_nonlinearweibull_", n, "_comptime.txt", sep = "")
  write.table(cov, dest_cov, col.names = F, row.names = F)
  write.table(width, dest_width, col.names = F, row.names = F)
  write.table(comptime, dest_comptime, col.names = F, row.names = F)
}
