##############################################################
##                                                          ##
##  SUPPLEMENTARY MATERIAL: DETAIL ASYMPTOTIC VARIANCE      ##
##                                                          ##
##----------------------------------------------------------##  
##                                                          ## 
## This code reproduces section 5 in the supplementary      ##
## material                                                 ##
##                                                          ## 
## Input                                                    ##    
##   - inpath: path containing the simulated datasets       ##
##   - outpath: destination path for figures                ##
##   - respath: path containing the results .txt files      ##
##                                                          ##
## Output                                                   ##
##   - PNG figures                                          ##
##                                                          ## 
## Date: March 5, 2019                                      ##
##############################################################
rm(list = ls())
expit <- function(x) exp(x) / (1 + exp(x))
library(DTRreg)
# inpath <- 
# outpath <- 
# respath <- 

#### ---- Run the simulations with n=1000, and extract estimated psi11 and standard errors ---- ####
n <- 1000
linear_correct_psi <- nonlinear_correct_psi <- matrix(NA, nrow = 1000, ncol = 2)
linear_incorrect_psi <- nonlinear_incorrect_psi <- matrix(NA, nrow = 1000, ncol = 2)
linear_correct_sepsi <- nonlinear_correct_sepsi <- matrix(NA, nrow = 1000, ncol = 4)
linear_incorrect_sepsi <- nonlinear_incorrect_sepsi <- matrix(NA, nrow = 1000, ncol = 4)
## linear treatment-free (1)
load(paste(inpath, "/datasets", n, "_normal_linear.RData", sep = ""))
for(i in 1:1000){
  original.data <- datasets[[i]]
  # asymptotic adjusted
  model1 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X12), cens.mod = list(delta ~ X12), var.estim = "asymptotic", asymp.opt = "adjusted", data = original.data)
  # asymptotic naive
  model2 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X12), cens.mod = list(delta ~ X12), var.estim = "asymptotic", asymp.opt = "naive", data = original.data)
  # Misspecified, asymptotic adjusted
  model3 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1), cens.mod = list(delta ~ X12), var.estim = "asymptotic", asymp.opt = "adjusted", data = original.data)
  # Misspecified, asymptotic naive
  model4 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1), cens.mod = list(delta ~ X12), var.estim = "asymptotic", asymp.opt = "naive", data = original.data)
  # psi estimates
  linear_correct_psi[i,] <- model1$psi[[1]]
  linear_incorrect_psi[i,] <- model3$psi[[1]]
  # se estimates
  linear_correct_sepsi[i,] <- c(sqrt(diag(model1$covmat[[1]])), sqrt(diag(model2$covmat[[1]])))
  linear_incorrect_sepsi[i,] <- c(sqrt(diag(model3$covmat[[1]])), sqrt(diag(model4$covmat[[1]])))
}
## nonlinear treatment-free (1)
load(paste(inpath, "/datasets", n, "_normal_nonlinear.RData", sep = ""))
for(i in 1:1000){
  original.data <- datasets[[i]]
  # asymptotic adjusted
  model1 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X14 + X12), cens.mod = list(delta ~ X12), var.estim = "asymptotic", asymp.opt = "adjusted", data = original.data)
  # asymptotic naive
  model2 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X14 + X12), cens.mod = list(delta ~ X12), var.estim = "asymptotic", asymp.opt = "naive", data = original.data)
  # Misspecified, asymptotic adjusted
  model3 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X12), cens.mod = list(delta ~ X12), var.estim = "asymptotic", asymp.opt = "adjusted", data = original.data)
  # Misspecified, asymptotic naive
  model4 <- dWSurv(time = list(~Y), blip.mod = list(~ X1), treat.mod = list(A1 ~ X1), tf.mod = list( ~ X1 + X12), cens.mod = list(delta ~ X12), var.estim = "asymptotic", asymp.opt = "naive", data = original.data)
  # psi estimates
  nonlinear_correct_psi[i,] <- model1$psi[[1]]
  nonlinear_incorrect_psi[i,] <- model3$psi[[1]]
  # se estimates
  nonlinear_correct_sepsi[i,] <- c(sqrt(diag(model1$covmat[[1]])), sqrt(diag(model2$covmat[[1]])))
  nonlinear_incorrect_sepsi[i,] <- c(sqrt(diag(model3$covmat[[1]])), sqrt(diag(model4$covmat[[1]])))
}
save.image(paste(respath, "/asymptotic_details.RData", sep = ""))

#### --------------- Distribution of psi11 linear/nonlinear correct/incorrect ----------- ####
cex.title <- 2
cex.lab <- 1.5
cex.axistitle <- 1.6
xrange <- range(c(linear_correct_psi[,2], nonlinear_correct_psi[,2],linear_incorrect_psi[,2], nonlinear_incorrect_psi[,2]))
png(paste(outpath, "/Figure6_Supp.png", sep = ""), width = 1355, height = 303)
par(mfrow = c(1,4), mar = c(4.1, 2.1, 3.1, 1.1))
hist(linear_correct_psi[,2], main = "Linear treatment-free", ylim = c(0,250), xlim = xrange, ylab = "", xaxt="n", breaks = seq(-0.30,0.55,by=0.05), xlab = expression(hat(psi)[11]), cex.axis = cex.lab, cex.lab = cex.axistitle, cex.main = cex.title)
axis(1, at=seq(-0.30,0.55,by=0.2), labels=seq(-0.30,0.55,by=0.2), cex.axis = cex.lab)
abline(v = 0.1, lwd = 3)
hist(linear_incorrect_psi[,2], main = "Misspecified linear treatment-free (1)", ylim = c(0,250), xlim = xrange, ylab = "", xaxt="n", breaks = seq(-0.30,0.55,by=0.05), xlab = expression(hat(psi)[11]), cex.axis = cex.lab, cex.lab = cex.axistitle, cex.main = cex.title)
axis(1, at=seq(-0.30,0.55,by=0.2), labels=seq(-0.30,0.55,by=0.2), cex.axis = cex.lab)
abline(v = 0.1, lwd = 3)
hist(nonlinear_correct_psi[,2], main = "Nonlinear treatment-free", ylim = c(0,250), xlim = xrange, xaxt="n", ylab = "", breaks = seq(-0.30,0.55,by=0.05), xlab = expression(hat(psi)[11]), cex.axis = cex.lab, cex.lab = cex.axistitle, cex.main = cex.title)
axis(1, at=seq(-0.30,0.55,by=0.2), labels=seq(-0.30,0.55,by=0.2), cex.axis = cex.lab)
abline(v = 0.1, lwd = 3)
hist(nonlinear_incorrect_psi[,2], main = "Misspecified nonlinear treatment-free (2)", ylim = c(0,250), xlim = xrange, xaxt="n", ylab = "", breaks = seq(-0.35,0.55,by=0.05), xlab = expression(hat(psi)[11]), cex.axis = cex.lab, cex.lab = cex.axistitle, cex.main = cex.title)
axis(1, at=seq(-0.35,0.55,by=0.2), labels=seq(-0.35,0.55,by=0.2), cex.axis = cex.lab)
abline(v = 0.1, lwd = 3)
dev.off()

#### --------------- Distribution of se psi11 linear/nonlinear correct/incorrect ----------- ####
seMC_linear_correct <- sd(linear_correct_psi[,2])
seMC_linear_incorrect <- sd(linear_incorrect_psi[,2])
seMC_nonlinear_correct <- sd(nonlinear_correct_psi[,2])
seMC_nonlinear_incorrect <- sd(nonlinear_incorrect_psi[,2])

png(paste(outpath, "/Figure7_Supp.png", sep = ""), width = 1355, height = 606)

par(mfrow = c(2,4), mar = c(2.1, 5.1, 3.1, 1.1))
hist(linear_correct_sepsi[,2], ylim = c(0,400), breaks = seq(0.06,0.115,by=0.005), main = "Linear treatment-free", xlab = "", ylab = "Asymptotic (adjusted)", cex.axis = cex.lab, cex.lab = cex.axistitle, cex.main = cex.title)
abline(v = seMC_linear_correct, lwd = 3)
par(mar = c(2.1, 2.1, 3.1, 1.1))
hist(linear_incorrect_sepsi[,2], main = "Misspecified linear treatment-free (1)", xlab = "", ylab = "", cex.axis = cex.lab, cex.lab = cex.axistitle, cex.main = cex.title)
abline(v = seMC_linear_incorrect, lwd = 3)
hist(nonlinear_correct_sepsi[,2], ylim = c(0,400), breaks = seq(0.06,0.11,by=0.005), main = "Nonlinear treatment-free", ylab = "", xlab = "", cex.axis = cex.lab, cex.lab = cex.axistitle, cex.main = cex.title)
abline(v = seMC_nonlinear_correct, lwd = 3)
hist(nonlinear_incorrect_sepsi[,2], breaks = seq(0.09,0.18,by=0.005), main = "Misspecified nonlinear treatment-free (2)", ylab = "", xlab = "", cex.axis = cex.lab, cex.lab = cex.axistitle, cex.main = cex.title)
abline(v = seMC_nonlinear_incorrect, lwd = 3)

par(mar = c(4.1, 5.1, 1.1, 1.1))
hist(linear_correct_sepsi[,4], ylim = c(0,400), breaks = seq(0.06,0.11,by=0.005), xlab = expression(hat(psi)[11]), ylab = "Asymptotic (naive)", cex.axis = cex.lab, cex.lab = cex.axistitle, main = "")
abline(v = seMC_linear_correct, lwd = 3)
par(mar = c(4.1, 2.1, 1.1, 1.1))
hist(linear_incorrect_sepsi[,4], ylim = c(0,400), xlab = expression(hat(psi)[11]), cex.axis = cex.lab, cex.lab = cex.axistitle, main = "")
abline(v = seMC_linear_incorrect, lwd = 3)
hist(nonlinear_correct_sepsi[,4], ylim = c(0,400), breaks = seq(0.06,0.11,by=0.005), xlab = expression(hat(psi)[11]), cex.axis = cex.lab, cex.lab = cex.axistitle, main = "")
abline(v = seMC_nonlinear_correct, lwd = 3)
hist(nonlinear_incorrect_sepsi[,4], breaks = seq(0.09,0.175,by=0.005), xlab = expression(hat(psi)[11]), cex.axis = cex.lab, cex.lab = cex.axistitle, main = "")
abline(v = seMC_nonlinear_incorrect, lwd = 3)
dev.off()


