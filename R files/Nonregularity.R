##############################################################
##                                                          ##
##                   NONREGULARITY                          ##
##                                                          ##
##----------------------------------------------------------##  
##                                                          ##  
## This code produces the txt files necessary to reproduce  ##
## the results presented in the mainpaper Section 4.4 and   ##
## in the supplementary material regarding the simulation   ##
## studies on nonregularity                                 ##
##                                                          ## 
##                                                          ## 
## Input                                                    ##    
##   - inpath: path of the simulated datasets               ##
##   - outpath: destination path for results .txt files     ##
##   - dtrpath: path with DTRreg_v1.4.R                     ##
##                                                          ##
## Output                                                   ##        
##   - three txt files per simulation setting               ##
##                                                          ##
## Date: May 07, 2019                                       ##
##############################################################
rm(list = ls())
expit <- function(x) exp(x) / (1 + exp(x))
# dtrpath <- 
# inpath <- 
# outpath <- 
source(paste(dtrpath, "/DTRreg_v1.4.R"))

#### ---- 8 scenarios, 2 stage DTR with Log-normal survival times, linear treatment-free ---- ####
## Corresponding scenario ID in Table 1 main paper
## Sample sizes: 300, 500, 1000, 10000
theta1 <- c(6.3, 0.5, -0.01, 0.1, 0.1)
for(sc in 1:8){
  for(n in c(300,500,1000,10000)){
    load(paste(inpath, "/datasets", n, "_sc", sc, ".RData", sep = ""))
    cov <- width <- matrix(0, nrow = 1000, ncol = 10)
    comptime <- matrix(NA, nrow = 1000, ncol = 15)
    for(i in 1:1000){
      original.data <- datasets[[i]]
      # Adjusted
      t0adj <- proc.time()
      model1 <- dWSurv(time = list(~Y1, ~Y2), blip.mod = list(~ X11, ~ X22), treat.mod = list(A1 ~ X11, A2 ~ X21), tf.mod = list( ~ X11 + X12, ~ X21 + X22 + X23 + X11), cens.mod = list(delta ~ X12, delta ~ X12), var.estim = "asymptotic", asymp.opt = "adjusted", data = original.data)
      t1adj <- proc.time()
      # Naive
      t0nai <- proc.time()
      model2 <- dWSurv(time = list(~Y1, ~Y2), blip.mod = list(~ X11, ~ X22), treat.mod = list(A1 ~ X11, A2 ~ X21), tf.mod = list( ~ X11 + X12, ~ X21 + X22 + X23 + X11), cens.mod = list(delta ~ X12, delta ~ X12), var.estim = "asymptotic", asymp.opt = "naive", data = original.data)
      t1nai <- proc.time()
      # Standard bootstrap
      t0stan <- proc.time()
      model3 <- dWSurv(time = list(~Y1, ~Y2), blip.mod = list(~ X11, ~ X22), treat.mod = list(A1 ~ X11, A2 ~ X21), tf.mod = list( ~ X11 + X12, ~ X21 + X22 + X23 + X11), cens.mod = list(delta ~ X12, delta ~ X12), var.estim = "bootstrap", boot.opt = "standard", B = 1000, data = original.data)
      t1stan <- proc.time()
      # Parametric bootstrap with empirical residual distribution
      t0para1 <- proc.time()
      model4 <- dWSurv(time = list(~Y1, ~Y2), blip.mod = list(~ X11, ~ X22), treat.mod = list(A1 ~ X11, A2 ~ X21), tf.mod = list( ~ X11 + X12, ~ X21 + X22 + X23 + X11), cens.mod = list(delta ~ X12, delta ~ X12), var.estim = "bootstrap", boot.opt = "empirical", B = 1000, data = original.data)
      t1para1 <- proc.time()
      # Parametric bootstrap with normal residual distribution
      t0para2 <- proc.time()
      model5 <- dWSurv(time = list(~Y1, ~Y2), blip.mod = list(~ X11, ~ X22), treat.mod = list(A1 ~ X11, A2 ~ X21), tf.mod = list( ~ X11 + X12, ~ X21 + X22 + X23 + X11), cens.mod = list(delta ~ X12, delta ~ X12), var.estim = "bootstrap", boot.opt = "normal", B = 1000, data = original.data)
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
    dest_cov <- paste(outpath, "/nonreg_sc", sc, "_", n, "_cov.txt", sep = "")
    dest_width <- paste(outpath, "/nonreg_sc", sc, "_", n, "_width.txt", sep = "")
    dest_comptime <- paste(outpath, "/nonreg_sc", sc, "_", n, "_comptime.txt", sep = "")
    write.table(cov, dest_cov, col.names = F, row.names = F)
    write.table(width, dest_width, col.names = F, row.names = F)
    write.table(comptime, dest_comptime, col.names = F, row.names = F)
  }
}