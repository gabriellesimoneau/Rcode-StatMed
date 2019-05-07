##############################################################
##                                                          ##
##     REPRODUCE FIGURES SUPPLEMENTARY MATERIAL             ##
##                                                          ##
##----------------------------------------------------------##  
##                                                          ## 
## This code reproduces the figures in the supplementary    ##
## material using the precompiled results located in        ##
## subfolder "Results'                                      ##
##                                                          ## 
## Input                                                    ##    
##   - inpath: path containing results .txt files           ##
##   - outpath: destination path for figures                ##
##                                                          ##
## Output                                                   ##
##   - PNG figures                                          ##
##                                                          ## 
## Date: March 5, 2019                                      ##
##############################################################

rm(list = ls())
# inpath <- 
# outpath <- 

#### --------- Figure 1: coverage for psi10, unknown error distribution, linear treatment-free  --------- ####
# import data and extract coverage
norm100 <- read.table(paste(inpath,"/error_linearnormal_100_cov.txt", sep = ""))
norm300 <- read.table(paste(inpath,"/error_linearnormal_300_cov.txt", sep = ""))
norm500 <- read.table(paste(inpath,"/error_linearnormal_500_cov.txt", sep = ""))
norm1000 <- read.table(paste(inpath,"/error_linearnormal_1000_cov.txt", sep = ""))
norm5000 <- read.table(paste(inpath,"/error_linearnormal_5000_cov.txt", sep = ""))
norm10000 <- read.table(paste(inpath,"/error_linearnormal_10000_cov.txt", sep = ""))
wei100 <- read.table(paste(inpath,"/error_linearweibull_100_cov.txt", sep = ""))
wei300 <- read.table(paste(inpath,"/error_linearweibull_300_cov.txt", sep = ""))
wei500 <- read.table(paste(inpath,"/error_linearweibull_500_cov.txt", sep = ""))
wei1000 <- read.table(paste(inpath,"/error_linearweibull_1000_cov.txt", sep = ""))
wei5000 <- read.table(paste(inpath,"/error_linearweibull_5000_cov.txt", sep = ""))
wei10000 <- read.table(paste(inpath,"/error_linearweibull_10000_cov.txt", sep = ""))
# extract coverage for psi11 only
keep <- seq(1,9,by=2)
norm100_cov <- apply(norm100[,keep], 2, mean)
norm300_cov <- apply(norm300[,keep], 2, mean)
norm500_cov <- apply(norm500[,keep], 2, mean)
norm1000_cov <- apply(norm1000[,keep], 2, mean)
norm5000_cov <- apply(norm5000[,keep], 2, mean)
norm10000_cov <- apply(norm10000[,keep], 2, mean)
wei100_cov <- apply(wei100[,keep], 2, mean)
wei300_cov <- apply(wei300[,keep], 2, mean)
wei500_cov <- apply(wei500[,keep], 2, mean)
wei1000_cov <- apply(wei1000[,keep], 2, mean)
wei5000_cov <- apply(wei5000[,keep], 2, mean)
wei10000_cov <- apply(wei10000[,keep], 2, mean)
# prepare plot
data <- expand.grid(error = c("Normal", "Weibull"), method = c("Asymptotic (adjusted)", "Asymptotic (naive)", "Standard bootstrap", "Parametric bootstrap 1", "Parametric bootstrap 2"), n = c(100, 300, 500, 1000, 5000, 10000))
data$n.coord <- rep(c(10, 20, 30, 40, 50, 60), each = 10)
data$coverage <- rep(NA, nrow(data))
data$coverage[which(data$error == "Normal")] <- c(norm100_cov, norm300_cov, norm500_cov, norm1000_cov, norm5000_cov, norm10000_cov)
data$coverage[which(data$error == "Weibull")] <-  c(wei100_cov, wei300_cov, wei500_cov, wei1000_cov, wei5000_cov, wei10000_cov)
cex.points <- 2.1
cex.lwd <- 1.9
cex.xaxislab <- cex.yaxislab <- 1.6
cex.title <- 2
cex.axistitle <- 1.7
png(paste(outpath, "/Figure1_Supp.png", sep = ""), width = 1136, height = 434)
mat <- matrix(c(2,1,3,1), nrow = 2)
layout(mat, c(1,1), c(8,1))
par(xpd = TRUE, mar = c(0,0,0,0))
plot(0,type='n',axes=FALSE,ann=FALSE, xlim = c(0,10))
legend(0,1,bty='n',c("Asymptotic (adjusted)"), lty = c("solid"), pch = c(15), cex = cex.axistitle-0.1, lwd = cex.lwd, pt.cex = cex.points, seg.len = 3)
legend(2.15,1,bty='n',c("Asymptotic (naive)"), lty = c("twodash"), pch = c(7), cex = cex.axistitle-0.1, lwd = cex.lwd, pt.cex = cex.points, seg.len = 3)
legend(4.05,1,bty='n',c("Standard bootstrap"), lty = c("longdash"), pch = c(16), cex = cex.axistitle-0.1, lwd = cex.lwd, pt.cex = cex.points, seg.len = 3)
legend(6.05,1,bty='n',c("Parametric bootstrap 1"), lty = c("dotted"), pch = c(10), cex = cex.axistitle-0.1, lwd = cex.lwd, pt.cex = cex.points, seg.len = 3)
legend(8.3,1,bty='n',c("Parametric bootstrap 2"), lty = c("dotdash"), pch = c(17), cex = cex.axistitle-0.1, lwd = cex.lwd, pt.cex = cex.points, seg.len = 3)
par(mar = c(4.1, 5.1, 3.1, 1.1), xpd = FALSE)
# normal
plot(x = data$n.coord[which(data$method == "Asymptotic (adjusted)" & data$error == "Normal")], y = data$coverage[which(data$method == "Asymptotic (adjusted)" & data$error == "Normal")], type = "l", ylim = c(0.9, 1.00), xlab = "Sample size", ylab = expression(paste("Coverage ", psi[10])), xaxt = "n", yaxt = "n", lty = "solid", lwd = cex.lwd, main = "Log-normal survival times", cex.lab = cex.axistitle, cex.main = cex.title)
axis(2, at = c(0.9, 0.95, 1.00), labels = c("0.90", "0.95", "1"), cex.axis = cex.yaxislab)
axis(1, at = c(10, 20, 30, 40, 50, 60), labels = c("100", "300", "500", "1000", "5000", "10,000"), cex.axis = cex.xaxislab)
abline(0.95, 0)
abline(0.95 - 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
abline(0.95 + 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
points(x = data$n.coord[which(data$method == "Asymptotic (adjusted)" & data$error == "Normal")], y = data$coverage[which(data$method == "Asymptotic (adjusted)" & data$error == "Normal")], pch = 15, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Asymptotic (naive)" & data$error == "Normal")], y = data$coverage[which(data$method == "Asymptotic (naive)" & data$error == "Normal")], lty = "twodash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Asymptotic (naive)" & data$error == "Normal")], y = data$coverage[which(data$method == "Asymptotic (naive)" & data$error == "Normal")], pch = 7, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Standard bootstrap" & data$error == "Normal")], y = data$coverage[which(data$method == "Standard bootstrap" & data$error == "Normal")], lty = "longdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Standard bootstrap" & data$error == "Normal")], y = data$coverage[which(data$method == "Standard bootstrap" & data$error == "Normal")], pch = 16, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 1" & data$error == "Normal")], y = data$coverage[which(data$method == "Parametric bootstrap 1" & data$error == "Normal")], lty = "dotted", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 1" & data$error == "Normal")], y = data$coverage[which(data$method == "Parametric bootstrap 1" & data$error == "Normal")], pch = 10, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 2" & data$error == "Normal")], y = data$coverage[which(data$method == "Parametric bootstrap 2" & data$error == "Normal")], lty = "dotdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 2" & data$error == "Normal")], y = data$coverage[which(data$method == "Parametric bootstrap 2" & data$error == "Normal")], pch = 17, cex = cex.points)
# weibull 
par(mar = c(4.1, 4.1, 3.1, 1.1), xpd = FALSE)
plot(x = data$n.coord[which(data$method == "Asymptotic (adjusted)" & data$error == "Weibull")], y = data$coverage[which(data$method == "Asymptotic (adjusted)" & data$error == "Weibull")], type = "l", ylim = c(0.9, 1.00), xlab = "Sample size", ylab = "", xaxt = "n", yaxt = "n", lty = "solid", lwd = cex.lwd, main = "Weibull survival times", cex.lab = cex.axistitle, cex.main = cex.title)
axis(2, at = c(0.9, 0.95, 1.00), labels = c("0.90", "0.95", "1"), cex.axis = cex.yaxislab)
axis(1, at = c(10, 20, 30, 40, 50, 60), labels = c("100", "300", "500", "1000", "5000", "10,000"), cex.axis = cex.xaxislab)
abline(0.95, 0)
abline(0.95 - 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
abline(0.95 + 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
points(x = data$n.coord[which(data$method == "Asymptotic (adjusted)" & data$error == "Weibull")], y = data$coverage[which(data$method == "Asymptotic (adjusted)" & data$error == "Weibull")], pch = 15, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Asymptotic (naive)" & data$error == "Weibull")], y = data$coverage[which(data$method == "Asymptotic (naive)" & data$error == "Weibull")], lty = "twodash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Asymptotic (naive)" & data$error == "Weibull")], y = data$coverage[which(data$method == "Asymptotic (naive)" & data$error == "Weibull")], pch = 7, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Standard bootstrap" & data$error == "Weibull")], y = data$coverage[which(data$method == "Standard bootstrap" & data$error == "Weibull")], lty = "longdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Standard bootstrap" & data$error == "Weibull")], y = data$coverage[which(data$method == "Standard bootstrap" & data$error == "Weibull")], pch = 16, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 1" & data$error == "Weibull")], y = data$coverage[which(data$method == "Parametric bootstrap 1" & data$error == "Weibull")], lty = "dotted", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 1" & data$error == "Weibull")], y = data$coverage[which(data$method == "Parametric bootstrap 1" & data$error == "Weibull")], pch = 10, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 2" & data$error == "Weibull")], y = data$coverage[which(data$method == "Parametric bootstrap 2" & data$error == "Weibull")], lty = "dotdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 2" & data$error == "Weibull")], y = data$coverage[which(data$method == "Parametric bootstrap 2" & data$error == "Weibull")], pch = 17, cex = cex.points)
dev.off()


#### --------- Figure 2: coverage for psi10, model misspecified, outcome 1/2  --------- ####
#import data
outcome1001 <- read.table(paste(inpath,"/misspec_tf1_100_cov.txt", sep = ""))
outcome5001 <- read.table(paste(inpath,"/misspec_tf1_500_cov.txt", sep = ""))
outcome10001 <- read.table(paste(inpath,"/misspec_tf1_1000_cov.txt", sep = ""))
outcome100001 <- read.table(paste(inpath,"/misspec_tf1_10000_cov.txt", sep = ""))
outcome1002 <- read.table(paste(inpath,"/misspec_tf2_100_cov.txt", sep = ""))
outcome5002 <- read.table(paste(inpath,"/misspec_tf2_500_cov.txt", sep = ""))
outcome10002 <- read.table(paste(inpath,"/misspec_tf2_1000_cov.txt", sep = ""))
outcome100002 <- read.table(paste(inpath,"/misspec_tf2_10000_cov.txt", sep = ""))
keep <- seq(1,9,by=2)
o1001 <- apply(outcome1001[,keep], 2, mean)
o5001 <- apply(outcome5001[,keep], 2, mean)
o10001 <- apply(outcome10001[,keep], 2, mean)
o100001 <- apply(outcome100001[,keep], 2, mean)
o1002 <- apply(outcome1002[,keep], 2, mean)
o5002 <- apply(outcome5002[,keep], 2, mean)
o10002 <- apply(outcome10002[,keep], 2, mean)
o100002 <- apply(outcome100002[,keep], 2, mean)
dat1 <- expand.grid(model = c("tf1", "tf2"), method = c("Asymptotic (adjusted)", "Asymptotic (naive)", "Standard bootstrap", "Parametric bootstrap 1", "Parametric bootstrap 2"), n = c(100, 500, 1000, 10000))
dat1$n.coord <- rep(c(10, 26.6, 42.6, 60), each = 10)
dat1$coverage <- rep(NA, nrow(dat1))
dat1$coverage[which(dat1$model == "tf1")] <- c(o1001, o5001, o10001, o100001)
dat1$coverage[which(dat1$model == "tf2")] <- c(o1002, o5002, o10002, o100002)
dat1$psi <- rep("psi10", nrow(dat1))
keep <- seq(2,10,by=2)
o1001 <- apply(outcome1001[,keep], 2, mean)
o5001 <- apply(outcome5001[,keep], 2, mean)
o10001 <- apply(outcome10001[,keep], 2, mean)
o100001 <- apply(outcome100001[,keep], 2, mean)
o1002 <- apply(outcome1002[,keep], 2, mean)
o5002 <- apply(outcome5002[,keep], 2, mean)
o10002 <- apply(outcome10002[,keep], 2, mean)
o100002 <- apply(outcome100002[,keep], 2, mean)
dat2 <- expand.grid(model = c("tf1", "tf2"), method = c("Asymptotic (adjusted)", "Asymptotic (naive)", "Standard bootstrap", "Parametric bootstrap 1", "Parametric bootstrap 2"), n = c(100, 500, 1000, 10000))
dat2$n.coord <- rep(c(10, 26.6, 42.6, 60), each = 10)
dat2$coverage <- rep(NA, nrow(dat2))
dat2$coverage[which(dat2$model == "tf1")] <- c(o1001, o5001, o10001, o100001)
dat2$coverage[which(dat2$model == "tf2")] <- c(o1002, o5002, o10002, o100002)
dat2$psi <- rep("psi11", nrow(dat2))
dat <- rbind(dat1, dat2)
# prepare plot
cex.points <- 2.5
cex.lwd <- 1.9
cex.xaxislab <- cex.yaxislab <- 1.7
cex.title <- 2.5
cex.axistitle <- 2
png(paste(outpath, "/Figure2_Supp.png", sep = ""), width = 1136, height = 700)
mat <- matrix(c(2,4,1,3,5,1), nrow = 3)
layout(mat, c(1,1,1), c(8,8,1))
par(xpd = TRUE, mar = c(0,0,0,0))
plot(0,type='n',axes=FALSE,ann=FALSE, xlim = c(0,10))
legend(0,1,bty='n',c("Asymptotic (adjusted)"), lty = c("solid"), pch = c(15), cex = cex.axistitle, lwd = cex.lwd, pt.cex = cex.points, seg.len = 3)
legend(2.15,1,bty='n',c("Asymptotic (naive)"), lty = c("twodash"), pch = c(7), cex = cex.axistitle, lwd = cex.lwd, pt.cex = cex.points, seg.len = 3)
legend(4.05,1,bty='n',c("Standard bootstrap"), lty = c("longdash"), pch = c(16), cex = cex.axistitle, lwd = cex.lwd, pt.cex = cex.points, seg.len = 3)
legend(6.05,1,bty='n',c("Parametric bootstrap 1"), lty = c("dotted"), pch = c(10), cex = cex.axistitle, lwd = cex.lwd, pt.cex = cex.points, seg.len = 3)
legend(8.3,1,bty='n',c("Parametric bootstrap 2"), lty = c("dotdash"), pch = c(17), cex = cex.axistitle, lwd = cex.lwd, pt.cex = cex.points, seg.len = 3)
# psi10 out1
data <- dat[which(dat$model == "tf1" & dat$psi == "psi10"),]
par(mar = c(2.1, 5.1, 3.1, 1.1), xpd = FALSE)
plot(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], type = "l", ylim = c(0.8, 1.00), xlab = "", ylab = expression(paste("Coverage ", psi[10])), xaxt = "n", yaxt = "n", lty = "solid", lwd = cex.lwd, main = "Treatment-free misspecified (1)", cex.lab = cex.axistitle, cex.main = cex.title)
axis(2, at = c(0.8, 0.85, 0.9, 0.95, 1.00), labels = c("0.80", "0.85", "0.90", "0.95", "1"), cex.axis = cex.yaxislab)
axis(1, at = c(10, 26.6, 42.6, 60), labels = c("100", "500", "1000", "10,000"), cex.axis = cex.xaxislab)
abline(0.95, 0)
abline(0.95 - 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
abline(0.95 + 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
points(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], pch = 15, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], lty = "twodash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], pch = 7, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], lty = "longdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], pch = 16, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], lty = "dotted", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], pch = 10, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], lty = "dotdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], pch = 17, cex = cex.points)
# psi10 out2
data <- dat[which(dat$model == "tf2" & dat$psi == "psi10"),]
par(mar = c(2.1, 3.1, 3.1, 1.1), xpd = FALSE)
plot(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], type = "l", ylim = c(0.9, 1.00), xlab = "", ylab = "", xaxt = "n", yaxt = "n", lty = "solid", lwd = cex.lwd, main = "Treatment-free misspecified (2)", cex.lab = cex.axistitle, cex.main = cex.title)
axis(2, at = c(0.9, 0.95, 1.00), labels = c("0.90", "0.95", "1"), cex.axis = cex.yaxislab)
axis(1, at = c(10, 26.6, 42.6, 60), labels = c("100", "500", "1000", "10,000"), cex.axis = cex.xaxislab)
abline(0.95, 0)
abline(0.95 - 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
abline(0.95 + 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
points(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], pch = 15, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], lty = "twodash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], pch = 7, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], lty = "longdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], pch = 16, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], lty = "dotted", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], pch = 10, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], lty = "dotdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], pch = 17, cex = cex.points)
# psi11 out1
data <- dat[which(dat$model == "tf1" & dat$psi == "psi11"),]
par(mar = c(4.1, 5.1, 1.1, 1.1), xpd = FALSE)
plot(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], type = "l", ylim = c(0.8, 1.00), xlab = "Sample size", ylab = expression(paste("Coverage ", psi[11])), xaxt = "n", yaxt = "n", lty = "solid", lwd = cex.lwd, cex.lab = cex.axistitle)
axis(2, at = c(0.8, 0.85, 0.9, 0.95, 1.00), labels = c("0.80", "0.85", "0.90", "0.95", "1"), cex.axis = cex.yaxislab)
axis(1, at = c(10, 26.6, 42.6, 60), labels = c("100", "500", "1000", "10,000"), cex.axis = cex.xaxislab)
abline(0.95, 0)
abline(0.95 - 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
abline(0.95 + 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
points(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], pch = 15, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], lty = "twodash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], pch = 7, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], lty = "longdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], pch = 16, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], lty = "dotted", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], pch = 10, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], lty = "dotdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], pch = 17, cex = cex.points)
# psi11 out2
data <- dat[which(dat$model == "tf2" & dat$psi == "psi11"),]
par(mar = c(4.1, 3.1, 1.1, 1.1), xpd = FALSE)
plot(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], type = "l", ylim = c(0.9, 1.00), xlab = "Sample size", ylab = "", xaxt = "n", yaxt = "n", lty = "solid", lwd = cex.lwd, cex.lab = cex.axistitle)
axis(2, at = c(0.9, 0.95, 1.00), labels = c("0.90", "0.95", "1"), cex.axis = cex.yaxislab)
axis(1, at = c(10, 26.6, 42.6, 60), labels = c("100", "500", "1000", "10,000"), cex.axis = cex.xaxislab)
abline(0.95, 0)
abline(0.95 - 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
abline(0.95 + 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
points(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], pch = 15, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], lty = "twodash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], pch = 7, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], lty = "longdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], pch = 16, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], lty = "dotted", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], pch = 10, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], lty = "dotdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], pch = 17, cex = cex.points)
dev.off()

#### --------- Figure 3: coverage for psi10, model misspecified, treatment 1/2  --------- ####
#import data
outcome1001 <- read.table(paste(inpath,"/misspec_trmt1_100_cov.txt", sep = ""))
outcome5001 <- read.table(paste(inpath,"/misspec_trmt1_500_cov.txt", sep = ""))
outcome10001 <- read.table(paste(inpath,"/misspec_trmt1_1000_cov.txt", sep = ""))
outcome100001 <- read.table(paste(inpath,"/misspec_trmt1_10000_cov.txt", sep = ""))
outcome1002 <- read.table(paste(inpath,"/misspec_trmt2_100_cov.txt", sep = ""))
outcome5002 <- read.table(paste(inpath,"/misspec_trmt2_500_cov.txt", sep = ""))
outcome10002 <- read.table(paste(inpath,"/misspec_trmt2_1000_cov.txt", sep = ""))
outcome100002 <- read.table(paste(inpath,"/misspec_trmt2_10000_cov.txt", sep = ""))
keep <- seq(1,9,by=2)
o1001 <- apply(outcome1001[,keep], 2, mean)
o5001 <- apply(outcome5001[,keep], 2, mean)
o10001 <- apply(outcome10001[,keep], 2, mean)
o100001 <- apply(outcome100001[,keep], 2, mean)
o1002 <- apply(outcome1002[,keep], 2, mean)
o5002 <- apply(outcome5002[,keep], 2, mean)
o10002 <- apply(outcome10002[,keep], 2, mean)
o100002 <- apply(outcome100002[,keep], 2, mean)
dat1 <- expand.grid(model = c("tf1", "tf2"), method = c("Asymptotic (adjusted)", "Asymptotic (naive)", "Standard bootstrap", "Parametric bootstrap 1", "Parametric bootstrap 2"), n = c(100, 500, 1000, 10000))
dat1$n.coord <- rep(c(10, 26.6, 42.6, 60), each = 10)
dat1$coverage <- rep(NA, nrow(dat1))
dat1$coverage[which(dat1$model == "tf1")] <- c(o1001, o5001, o10001, o100001)
dat1$coverage[which(dat1$model == "tf2")] <- c(o1002, o5002, o10002, o100002)
dat1$psi <- rep("psi10", nrow(dat1))
keep <- seq(2,10,by=2)
o1001 <- apply(outcome1001[,keep], 2, mean)
o5001 <- apply(outcome5001[,keep], 2, mean)
o10001 <- apply(outcome10001[,keep], 2, mean)
o100001 <- apply(outcome100001[,keep], 2, mean)
o1002 <- apply(outcome1002[,keep], 2, mean)
o5002 <- apply(outcome5002[,keep], 2, mean)
o10002 <- apply(outcome10002[,keep], 2, mean)
o100002 <- apply(outcome100002[,keep], 2, mean)
dat2 <- expand.grid(model = c("tf1", "tf2"), method = c("Asymptotic (adjusted)", "Asymptotic (naive)", "Standard bootstrap", "Parametric bootstrap 1", "Parametric bootstrap 2"), n = c(100, 500, 1000, 10000))
dat2$n.coord <- rep(c(10, 26.6, 42.6, 60), each = 10)
dat2$coverage <- rep(NA, nrow(dat2))
dat2$coverage[which(dat2$model == "tf1")] <- c(o1001, o5001, o10001, o100001)
dat2$coverage[which(dat2$model == "tf2")] <- c(o1002, o5002, o10002, o100002)
dat2$psi <- rep("psi11", nrow(dat2))
dat <- rbind(dat1, dat2)
# prepare plot
cex.points <- 2.5
cex.lwd <- 1.9
cex.xaxislab <- cex.yaxislab <- 1.7
cex.title <- 2.5
cex.axistitle <- 2
png(paste(outpath, "/Figure3_Supp.png", sep = ""), width = 1136, height = 700)
mat <- matrix(c(2,4,1,3,5,1), nrow = 3)
layout(mat, c(1,1,1), c(8,8,1))
par(xpd = TRUE, mar = c(0,0,0,0))
plot(0,type='n',axes=FALSE,ann=FALSE, xlim = c(0,10))
legend(0,1,bty='n',c("Asymptotic (adjusted)"), lty = c("solid"), pch = c(15), cex = cex.axistitle, lwd = cex.lwd, pt.cex = cex.points, seg.len = 3)
legend(2.15,1,bty='n',c("Asymptotic (naive)"), lty = c("twodash"), pch = c(7), cex = cex.axistitle, lwd = cex.lwd, pt.cex = cex.points, seg.len = 3)
legend(4.05,1,bty='n',c("Standard bootstrap"), lty = c("longdash"), pch = c(16), cex = cex.axistitle, lwd = cex.lwd, pt.cex = cex.points, seg.len = 3)
legend(6.05,1,bty='n',c("Parametric bootstrap 1"), lty = c("dotted"), pch = c(10), cex = cex.axistitle, lwd = cex.lwd, pt.cex = cex.points, seg.len = 3)
legend(8.3,1,bty='n',c("Parametric bootstrap 2"), lty = c("dotdash"), pch = c(17), cex = cex.axistitle, lwd = cex.lwd, pt.cex = cex.points, seg.len = 3)
# psi10 out1
data <- dat[which(dat$model == "tf1" & dat$psi == "psi10"),]
par(mar = c(2.1, 5.1, 3.1, 1.1), xpd = FALSE)
plot(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], type = "l", ylim = c(0.9, 1.00), xlab = "", ylab = expression(paste("Coverage ", psi[10])), xaxt = "n", yaxt = "n", lty = "solid", lwd = cex.lwd, main = "Treatment misspecified (1)", cex.lab = cex.axistitle, cex.main = cex.title)
axis(2, at = c(0.9, 0.95, 1.00), labels = c("0.90", "0.95", "1"), cex.axis = cex.yaxislab)
axis(1, at = c(10, 26.6, 42.6, 60), labels = c("100", "500", "1000", "10,000"), cex.axis = cex.xaxislab)
abline(0.95, 0)
abline(0.95 - 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
abline(0.95 + 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
points(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], pch = 15, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], lty = "twodash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], pch = 7, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], lty = "longdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], pch = 16, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], lty = "dotted", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], pch = 10, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], lty = "dotdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], pch = 17, cex = cex.points)
# psi10 out2
data <- dat[which(dat$model == "tf2" & dat$psi == "psi10"),]
par(mar = c(2.1, 3.1, 3.1, 1.1), xpd = FALSE)
plot(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], type = "l", ylim = c(0.9, 1.00), xlab = "", ylab = "", xaxt = "n", yaxt = "n", lty = "solid", lwd = cex.lwd, main = "Treatment misspecified (2)", cex.lab = cex.axistitle, cex.main = cex.title)
axis(2, at = c(0.9, 0.95, 1.00), labels = c("0.90", "0.95", "1"), cex.axis = cex.yaxislab)
axis(1, at = c(10, 26.6, 42.6, 60), labels = c("100", "500", "1000", "10,000"), cex.axis = cex.xaxislab)
abline(0.95, 0)
abline(0.95 - 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
abline(0.95 + 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
points(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], pch = 15, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], lty = "twodash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], pch = 7, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], lty = "longdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], pch = 16, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], lty = "dotted", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], pch = 10, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], lty = "dotdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], pch = 17, cex = cex.points)
# psi10 out1
data <- dat[which(dat$model == "tf1" & dat$psi == "psi11"),]
par(mar = c(4.1, 5.1, 1.1, 1.1), xpd = FALSE)
plot(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], type = "l", ylim = c(0.9, 1.00), xlab = "Sample size", ylab = expression(paste("Coverage ", psi[11])), xaxt = "n", yaxt = "n", lty = "solid", lwd = cex.lwd, cex.lab = cex.axistitle)
axis(2, at = c(0.9, 0.95, 1.00), labels = c("0.90", "0.95", "1"), cex.axis = cex.yaxislab)
axis(1, at = c(10, 26.6, 42.6, 60), labels = c("100", "500", "1000", "10,000"), cex.axis = cex.xaxislab)
abline(0.95, 0)
abline(0.95 - 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
abline(0.95 + 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
points(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], pch = 15, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], lty = "twodash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], pch = 7, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], lty = "longdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], pch = 16, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], lty = "dotted", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], pch = 10, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], lty = "dotdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], pch = 17, cex = cex.points)
# psi11 out2
data <- dat[which(dat$model == "tf2" & dat$psi == "psi11"),]
par(mar = c(4.1, 3.1, 1.1, 1.1), xpd = FALSE)
plot(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], type = "l", ylim = c(0.9, 1.00), xlab = "Sample size", ylab = "", xaxt = "n", yaxt = "n", lty = "solid", lwd = cex.lwd, cex.lab = cex.axistitle)
axis(2, at = c(0.9, 0.95, 1.00), labels = c("0.90", "0.95", "1"), cex.axis = cex.yaxislab)
axis(1, at = c(10, 26.6, 42.6, 60), labels = c("100", "500", "1000", "10,000"), cex.axis = cex.xaxislab)
abline(0.95, 0)
abline(0.95 - 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
abline(0.95 + 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
points(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], pch = 15, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], lty = "twodash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], pch = 7, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], lty = "longdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], pch = 16, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], lty = "dotted", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], pch = 10, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], lty = "dotdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], pch = 17, cex = cex.points)
dev.off()


#### --------- Figure 4: coverage for psi10, model misspecified, censoring 1/2  --------- ####
#import data
outcome1001 <- read.table(paste(inpath,"/misspec_cens1_100_cov.txt", sep = ""))
outcome5001 <- read.table(paste(inpath,"/misspec_cens1_500_cov.txt", sep = ""))
outcome10001 <- read.table(paste(inpath,"/misspec_cens1_1000_cov.txt", sep = ""))
outcome100001 <- read.table(paste(inpath,"/misspec_cens1_10000_cov.txt", sep = ""))
outcome1002 <- read.table(paste(inpath,"/misspec_cens2_100_cov.txt", sep = ""))
outcome5002 <- read.table(paste(inpath,"/misspec_cens2_500_cov.txt", sep = ""))
outcome10002 <- read.table(paste(inpath,"/misspec_cens2_1000_cov.txt", sep = ""))
outcome100002 <- read.table(paste(inpath,"/misspec_cens2_10000_cov.txt", sep = ""))
keep <- seq(1,9,by=2)
o1001 <- apply(outcome1001[,keep], 2, mean)
o5001 <- apply(outcome5001[,keep], 2, mean)
o10001 <- apply(outcome10001[,keep], 2, mean)
o100001 <- apply(outcome100001[,keep], 2, mean)
o1002 <- apply(outcome1002[,keep], 2, mean)
o5002 <- apply(outcome5002[,keep], 2, mean)
o10002 <- apply(outcome10002[,keep], 2, mean)
o100002 <- apply(outcome100002[,keep], 2, mean)
dat1 <- expand.grid(model = c("tf1", "tf2"), method = c("Asymptotic (adjusted)", "Asymptotic (naive)", "Standard bootstrap", "Parametric bootstrap 1", "Parametric bootstrap 2"), n = c(100, 500, 1000, 10000))
dat1$n.coord <- rep(c(10, 26.6, 42.6, 60), each = 10)
dat1$coverage <- rep(NA, nrow(dat1))
dat1$coverage[which(dat1$model == "tf1")] <- c(o1001, o5001, o10001, o100001)
dat1$coverage[which(dat1$model == "tf2")] <- c(o1002, o5002, o10002, o100002)
dat1$psi <- rep("psi10", nrow(dat1))
keep <- seq(2,10,by=2)
o1001 <- apply(outcome1001[,keep], 2, mean)
o5001 <- apply(outcome5001[,keep], 2, mean)
o10001 <- apply(outcome10001[,keep], 2, mean)
o100001 <- apply(outcome100001[,keep], 2, mean)
o1002 <- apply(outcome1002[,keep], 2, mean)
o5002 <- apply(outcome5002[,keep], 2, mean)
o10002 <- apply(outcome10002[,keep], 2, mean)
o100002 <- apply(outcome100002[,keep], 2, mean)
dat2 <- expand.grid(model = c("tf1", "tf2"), method = c("Asymptotic (adjusted)", "Asymptotic (naive)", "Standard bootstrap", "Parametric bootstrap 1", "Parametric bootstrap 2"), n = c(100, 500, 1000, 10000))
dat2$n.coord <- rep(c(10, 26.6, 42.6, 60), each = 10)
dat2$coverage <- rep(NA, nrow(dat2))
dat2$coverage[which(dat2$model == "tf1")] <- c(o1001, o5001, o10001, o100001)
dat2$coverage[which(dat2$model == "tf2")] <- c(o1002, o5002, o10002, o100002)
dat2$psi <- rep("psi11", nrow(dat2))
dat <- rbind(dat1, dat2)
# prepare plot
cex.points <- 2.5
cex.lwd <- 1.9
cex.xaxislab <- cex.yaxislab <- 1.7
cex.title <- 2.5
cex.axistitle <- 2
png(paste(outpath, "/Figure4_Supp.png", sep = ""), width = 1136, height = 700)
mat <- matrix(c(2,4,1,3,5,1), nrow = 3)
layout(mat, c(1,1,1), c(8,8,1))
par(xpd = TRUE, mar = c(0,0,0,0))
plot(0,type='n',axes=FALSE,ann=FALSE, xlim = c(0,10))
legend(0,1,bty='n',c("Asymptotic (adjusted)"), lty = c("solid"), pch = c(15), cex = cex.axistitle, lwd = cex.lwd, pt.cex = cex.points, seg.len = 3)
legend(2.15,1,bty='n',c("Asymptotic (naive)"), lty = c("twodash"), pch = c(7), cex = cex.axistitle, lwd = cex.lwd, pt.cex = cex.points, seg.len = 3)
legend(4.05,1,bty='n',c("Standard bootstrap"), lty = c("longdash"), pch = c(16), cex = cex.axistitle, lwd = cex.lwd, pt.cex = cex.points, seg.len = 3)
legend(6.05,1,bty='n',c("Parametric bootstrap 1"), lty = c("dotted"), pch = c(10), cex = cex.axistitle, lwd = cex.lwd, pt.cex = cex.points, seg.len = 3)
legend(8.3,1,bty='n',c("Parametric bootstrap 2"), lty = c("dotdash"), pch = c(17), cex = cex.axistitle, lwd = cex.lwd, pt.cex = cex.points, seg.len = 3)
# psi10 out1
data <- dat[which(dat$model == "tf1" & dat$psi == "psi10"),]
par(mar = c(2.1, 5.1, 3.1, 1.1), xpd = FALSE)
plot(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], type = "l", ylim = c(0.9, 1.00), xlab = "", ylab = expression(paste("Coverage ", psi[10])), xaxt = "n", yaxt = "n", lty = "solid", lwd = cex.lwd, main = "Censoring misspecified (1)", cex.lab = cex.axistitle, cex.main = cex.title)
axis(2, at = c(0.9, 0.95, 1.00), labels = c("0.90", "0.95", "1"), cex.axis = cex.yaxislab)
axis(1, at = c(10, 26.6, 42.6, 60), labels = c("100", "500", "1000", "10,000"), cex.axis = cex.xaxislab)
abline(0.95, 0)
abline(0.95 - 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
abline(0.95 + 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
points(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], pch = 15, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], lty = "twodash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], pch = 7, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], lty = "longdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], pch = 16, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], lty = "dotted", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], pch = 10, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], lty = "dotdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], pch = 17, cex = cex.points)
# psi10 out2
data <- dat[which(dat$model == "tf2" & dat$psi == "psi10"),]
par(mar = c(2.1, 3.1, 3.1, 1.1), xpd = FALSE)
plot(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], type = "l", ylim = c(0.9, 1.00), xlab = "", ylab = "", xaxt = "n", yaxt = "n", lty = "solid", lwd = cex.lwd, main = "Censoring misspecified (2)", cex.lab = cex.axistitle, cex.main = cex.title)
axis(2, at = c(0.9, 0.95, 1.00), labels = c("0.90", "0.95", "1"), cex.axis = cex.yaxislab)
axis(1, at = c(10, 26.6, 42.6, 60), labels = c("100", "500", "1000", "10,000"), cex.axis = cex.xaxislab)
abline(0.95, 0)
abline(0.95 - 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
abline(0.95 + 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
points(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], pch = 15, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], lty = "twodash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], pch = 7, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], lty = "longdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], pch = 16, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], lty = "dotted", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], pch = 10, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], lty = "dotdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], pch = 17, cex = cex.points)
# psi10 out1
data <- dat[which(dat$model == "tf1" & dat$psi == "psi11"),]
par(mar = c(4.1, 5.1, 1.1, 1.1), xpd = FALSE)
plot(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], type = "l", ylim = c(0.9, 1.00), xlab = "Sample size", ylab = expression(paste("Coverage ", psi[11])), xaxt = "n", yaxt = "n", lty = "solid", lwd = cex.lwd, cex.lab = cex.axistitle)
axis(2, at = c(0.9, 0.95, 1.00), labels = c("0.90", "0.95", "1"), cex.axis = cex.yaxislab)
axis(1, at = c(10, 26.6, 42.6, 60), labels = c("100", "500", "1000", "10,000"), cex.axis = cex.xaxislab)
abline(0.95, 0)
abline(0.95 - 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
abline(0.95 + 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
points(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], pch = 15, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], lty = "twodash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], pch = 7, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], lty = "longdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], pch = 16, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], lty = "dotted", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], pch = 10, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], lty = "dotdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], pch = 17, cex = cex.points)
# psi11 out2
data <- dat[which(dat$model == "tf2" & dat$psi == "psi11"),]
par(mar = c(4.1, 3.1, 1.1, 1.1), xpd = FALSE)
plot(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], type = "l", ylim = c(0.9, 1.00), xlab = "Sample size", ylab = "", xaxt = "n", yaxt = "n", lty = "solid", lwd = cex.lwd, cex.lab = cex.axistitle)
axis(2, at = c(0.9, 0.95, 1.00), labels = c("0.90", "0.95", "1"), cex.axis = cex.yaxislab)
axis(1, at = c(10, 26.6, 42.6, 60), labels = c("100", "500", "1000", "10,000"), cex.axis = cex.xaxislab)
abline(0.95, 0)
abline(0.95 - 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
abline(0.95 + 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
points(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], pch = 15, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], lty = "twodash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], pch = 7, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], lty = "longdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], pch = 16, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], lty = "dotted", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], pch = 10, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], lty = "dotdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], pch = 17, cex = cex.points)
dev.off()


#### --------- Figure 5: coverage for psi10, all 8 nonregular scenarios  ---------####
# import and prepare data
g <- 0.15
keep <- seq(1,9,by=2)
n <- 300
sc1 <- read.table(paste(inpath,"/nonreg_sc1_", n,"_cov.txt", sep = ""))
sc2 <- read.table(paste(inpath,"/nonreg_sc2_", n,"_cov.txt", sep = ""))
sc3 <- read.table(paste(inpath,"/nonreg_sc3_", n,"_cov.txt", sep = ""))
sc4 <- read.table(paste(inpath,"/nonreg_sc4_", n,"_cov.txt", sep = ""))
sc5 <- read.table(paste(inpath,"/nonreg_sc5_", n,"_cov.txt", sep = ""))
sc6 <- read.table(paste(inpath,"/nonreg_sc6_", n,"_cov.txt", sep = ""))
sc7 <- read.table(paste(inpath,"/nonreg_sc7_", n,"_cov.txt", sep = ""))
sc8 <- read.table(paste(inpath,"/nonreg_sc8_", n,"_cov.txt", sep = ""))
cov1 <- apply(sc1, 2, sum)[keep]
cov2 <- apply(sc2, 2, sum)[keep]
cov3 <- apply(sc3, 2, sum)[keep]
cov4 <- apply(sc4, 2, sum)[keep]
cov5 <- apply(sc5, 2, sum)[keep]
cov6 <- apply(sc6, 2, sum)[keep]
cov7 <- apply(sc7, 2, sum)[keep]
cov8 <- apply(sc8, 2, sum)[keep]
coverage <- c(cov1, cov2, cov3, cov4, cov5, cov6, cov7, cov8)/1000
method <- factor(rep(c("Asymptotic (adjusted)", "Asymptotic (naive)", "Standard bootstrap", "Parametric bootstrap 1", "Parametric bootstrap 2"), 8), levels = c("Asymptotic (adjusted)", "Asymptotic (naive)", "Standard bootstrap", "Parametric bootstrap 1", "Parametric bootstrap 2"))
scenario <- factor(rep(c("Scenario 1 (Nonregular)", "Scenario 2 (Near-nonregular)", "Scenario 3 (Nonregular)", "Scenario 4 (Near-nonregular)", "Scenario 5 (Nonregular)", "Scenario 6 (Near-nonregular)", "Scenario 7 (Regular)", "Scenario 8 (Regular)"), each = 5), levels = c("Scenario 1 (Nonregular)", "Scenario 2 (Near-nonregular)", "Scenario 3 (Nonregular)", "Scenario 4 (Near-nonregular)", "Scenario 5 (Nonregular)", "Scenario 6 (Near-nonregular)", "Scenario 7 (Regular)", "Scenario 8 (Regular)"))
data300 <- data.frame(coverage, method, scenario)
data300$n <- rep(300, nrow(data300))
data300$n.coord <- rep(10, nrow(data300))
n <- 500
sc1 <- read.table(paste(inpath,"/nonreg_sc1_", n,"_cov.txt", sep = ""))
sc2 <- read.table(paste(inpath,"/nonreg_sc2_", n,"_cov.txt", sep = ""))
sc3 <- read.table(paste(inpath,"/nonreg_sc3_", n,"_cov.txt", sep = ""))
sc4 <- read.table(paste(inpath,"/nonreg_sc4_", n,"_cov.txt", sep = ""))
sc5 <- read.table(paste(inpath,"/nonreg_sc5_", n,"_cov.txt", sep = ""))
sc6 <- read.table(paste(inpath,"/nonreg_sc6_", n,"_cov.txt", sep = ""))
sc7 <- read.table(paste(inpath,"/nonreg_sc7_", n,"_cov.txt", sep = ""))
sc8 <- read.table(paste(inpath,"/nonreg_sc8_", n,"_cov.txt", sep = ""))
cov1 <- apply(sc1, 2, sum)[keep]
cov2 <- apply(sc2, 2, sum)[keep]
cov3 <- apply(sc3, 2, sum)[keep]
cov4 <- apply(sc4, 2, sum)[keep]
cov5 <- apply(sc5, 2, sum)[keep]
cov6 <- apply(sc6, 2, sum)[keep]
cov7 <- apply(sc7, 2, sum)[keep]
cov8 <- apply(sc8, 2, sum)[keep]
coverage <- c(cov1, cov2, cov3, cov4, cov5, cov6, cov7, cov8)/1000
method <- factor(rep(c("Asymptotic (adjusted)", "Asymptotic (naive)", "Standard bootstrap", "Parametric bootstrap 1", "Parametric bootstrap 2"), 8), levels = c("Asymptotic (adjusted)", "Asymptotic (naive)", "Standard bootstrap", "Parametric bootstrap 1", "Parametric bootstrap 2"))
scenario <- factor(rep(c("Scenario 1 (Nonregular)", "Scenario 2 (Near-nonregular)", "Scenario 3 (Nonregular)", "Scenario 4 (Near-nonregular)", "Scenario 5 (Nonregular)", "Scenario 6 (Near-nonregular)", "Scenario 7 (Regular)", "Scenario 8 (Regular)"), each = 5), levels = c("Scenario 1 (Nonregular)", "Scenario 2 (Near-nonregular)", "Scenario 3 (Nonregular)", "Scenario 4 (Near-nonregular)", "Scenario 5 (Nonregular)", "Scenario 6 (Near-nonregular)", "Scenario 7 (Regular)", "Scenario 8 (Regular)"))
data500 <- data.frame(coverage, method, scenario)
data500$n <- rep(500, nrow(data500))
data500$n.coord <- rep(26.6, nrow(data500))
n <- 1000
sc1 <- read.table(paste(inpath,"/nonreg_sc1_", n,"_cov.txt", sep = ""))
sc2 <- read.table(paste(inpath,"/nonreg_sc2_", n,"_cov.txt", sep = ""))
sc3 <- read.table(paste(inpath,"/nonreg_sc3_", n,"_cov.txt", sep = ""))
sc4 <- read.table(paste(inpath,"/nonreg_sc4_", n,"_cov.txt", sep = ""))
sc5 <- read.table(paste(inpath,"/nonreg_sc5_", n,"_cov.txt", sep = ""))
sc6 <- read.table(paste(inpath,"/nonreg_sc6_", n,"_cov.txt", sep = ""))
sc7 <- read.table(paste(inpath,"/nonreg_sc7_", n,"_cov.txt", sep = ""))
sc8 <- read.table(paste(inpath,"/nonreg_sc8_", n,"_cov.txt", sep = ""))
cov1 <- apply(sc1, 2, sum)[keep]
cov2 <- apply(sc2, 2, sum)[keep]
cov3 <- apply(sc3, 2, sum)[keep]
cov4 <- apply(sc4, 2, sum)[keep]
cov5 <- apply(sc5, 2, sum)[keep]
cov6 <- apply(sc6, 2, sum)[keep]
cov7 <- apply(sc7, 2, sum)[keep]
cov8 <- apply(sc8, 2, sum)[keep]
coverage <- c(cov1, cov2, cov3, cov4, cov5, cov6, cov7, cov8)/1000
method <- factor(rep(c("Asymptotic (adjusted)", "Asymptotic (naive)", "Standard bootstrap", "Parametric bootstrap 1", "Parametric bootstrap 2"), 8), levels = c("Asymptotic (adjusted)", "Asymptotic (naive)", "Standard bootstrap", "Parametric bootstrap 1", "Parametric bootstrap 2"))
scenario <- factor(rep(c("Scenario 1 (Nonregular)", "Scenario 2 (Near-nonregular)", "Scenario 3 (Nonregular)", "Scenario 4 (Near-nonregular)", "Scenario 5 (Nonregular)", "Scenario 6 (Near-nonregular)", "Scenario 7 (Regular)", "Scenario 8 (Regular)"), each = 5), levels = c("Scenario 1 (Nonregular)", "Scenario 2 (Near-nonregular)", "Scenario 3 (Nonregular)", "Scenario 4 (Near-nonregular)", "Scenario 5 (Nonregular)", "Scenario 6 (Near-nonregular)", "Scenario 7 (Regular)", "Scenario 8 (Regular)"))
data1000 <- data.frame(coverage, method, scenario)
data1000$n <- rep(1000, nrow(data1000))
data1000$n.coord <- rep(42.6, nrow(data1000))
n <- 10000
sc1 <- read.table(paste(inpath,"/nonreg_sc1_", n,"_cov.txt", sep = ""))
sc2 <- read.table(paste(inpath,"/nonreg_sc2_", n,"_cov.txt", sep = ""))
sc3 <- read.table(paste(inpath,"/nonreg_sc3_", n,"_cov.txt", sep = ""))
sc4 <- read.table(paste(inpath,"/nonreg_sc4_", n,"_cov.txt", sep = ""))
sc5 <- read.table(paste(inpath,"/nonreg_sc5_", n,"_cov.txt", sep = ""))
sc6 <- read.table(paste(inpath,"/nonreg_sc6_", n,"_cov.txt", sep = ""))
sc7 <- read.table(paste(inpath,"/nonreg_sc7_", n,"_cov.txt", sep = ""))
sc8 <- read.table(paste(inpath,"/nonreg_sc8_", n,"_cov.txt", sep = ""))
cov1 <- apply(sc1, 2, sum)[keep]
cov2 <- apply(sc2, 2, sum)[keep]
cov3 <- apply(sc3, 2, sum)[keep]
cov4 <- apply(sc4, 2, sum)[keep]
cov5 <- apply(sc5, 2, sum)[keep]
cov6 <- apply(sc6, 2, sum)[keep]
cov7 <- apply(sc7, 2, sum)[keep]
cov8 <- apply(sc8, 2, sum)[keep]
coverage <- c(cov1, cov2, cov3, cov4, cov5, cov6, cov7, cov8)/1000
method <- factor(rep(c("Asymptotic (adjusted)", "Asymptotic (naive)", "Standard bootstrap", "Parametric bootstrap 1", "Parametric bootstrap 2"), 8), levels = c("Asymptotic (adjusted)", "Asymptotic (naive)", "Standard bootstrap", "Parametric bootstrap 1", "Parametric bootstrap 2"))
scenario <- factor(rep(c("Scenario 1 (Nonregular)", "Scenario 2 (Near-nonregular)", "Scenario 3 (Nonregular)", "Scenario 4 (Near-nonregular)", "Scenario 5 (Nonregular)", "Scenario 6 (Near-nonregular)", "Scenario 7 (Regular)", "Scenario 8 (Regular)"), each = 5), levels = c("Scenario 1 (Nonregular)", "Scenario 2 (Near-nonregular)", "Scenario 3 (Nonregular)", "Scenario 4 (Near-nonregular)", "Scenario 5 (Nonregular)", "Scenario 6 (Near-nonregular)", "Scenario 7 (Regular)", "Scenario 8 (Regular)"))
data10000 <- data.frame(coverage, method, scenario)
data10000$n <- rep(10000, nrow(data10000))
data10000$n.coord <- rep(60, nrow(data10000))
dat <- rbind(data300, data500, data1000, data10000)
# prepare plot
cex.points <- 2.5
cex.lwd <- 1.9
cex.xaxislab <- cex.yaxislab <- 1.7
cex.title <- 2.5
cex.axistitle <- 2
png(paste(outpath, "/Figure5_Supp.png", sep = ""), width = 1136, height = 1200)
mat <- matrix(c(2,4,6,8,1,3,5,7,9,1), nrow = 5)
layout(mat, c(1,1,1,1,1), c(8,8,8,8,1))
par(xpd = TRUE, mar = c(0,0,0,0))
plot(0,type='n',axes=FALSE,ann=FALSE, xlim = c(0,10))
legend(0,1,bty='n',c("Asymptotic (adjusted)"), lty = c("solid"), pch = c(15), cex = cex.axistitle, lwd = cex.lwd, pt.cex = cex.points, seg.len = 3)
legend(2.15,1,bty='n',c("Asymptotic (naive)"), lty = c("twodash"), pch = c(7), cex = cex.axistitle, lwd = cex.lwd, pt.cex = cex.points, seg.len = 3)
legend(4.05,1,bty='n',c("Standard bootstrap"), lty = c("longdash"), pch = c(16), cex = cex.axistitle, lwd = cex.lwd, pt.cex = cex.points, seg.len = 3)
legend(6.05,1,bty='n',c("Parametric bootstrap 1"), lty = c("dotted"), pch = c(10), cex = cex.axistitle, lwd = cex.lwd, pt.cex = cex.points, seg.len = 3)
legend(8.3,1,bty='n',c("Parametric bootstrap 2"), lty = c("dotdash"), pch = c(17), cex = cex.axistitle, lwd = cex.lwd, pt.cex = cex.points, seg.len = 3)
# psi10 sc1
data <- dat[which(dat$scenario == "Scenario 1 (Nonregular)"),]
par(mar = c(2.1, 5.1, 3.1, 1.1), xpd = FALSE)
plot(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], type = "l", ylim = c(0.9, 1.00), xlab = "", ylab = expression(paste("Coverage ", psi[10])), xaxt = "n", yaxt = "n", lty = "solid", lwd = cex.lwd, main = "Scenario 1 (Nonregular)", cex.lab = cex.axistitle, cex.main = cex.title)
axis(2, at = c(0.9, 0.95, 1.00), labels = c("0.90", "0.95", "1"), cex.axis = cex.yaxislab)
axis(1, at = c(10, 26.6, 42.6, 60), labels = c("300", "500", "1000", "10,000"), cex.axis = cex.xaxislab)
abline(0.95, 0)
abline(0.95 - 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
abline(0.95 + 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
points(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], pch = 15, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], lty = "twodash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], pch = 7, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], lty = "longdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], pch = 16, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], lty = "dotted", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], pch = 10, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], lty = "dotdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], pch = 17, cex = cex.points)
# psi10 sc2
data <- dat[which(dat$scenario == "Scenario 2 (Near-nonregular)"),]
plot(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], type = "l", ylim = c(0.9, 1.00), xlab = "", ylab = expression(paste("Coverage ", psi[10])), xaxt = "n", yaxt = "n", lty = "solid", lwd = cex.lwd, main = "Scenario 2 (Near-nonregular)", cex.lab = cex.axistitle, cex.main = cex.title)
axis(2, at = c(0.9, 0.95, 1.00), labels = c("0.90", "0.95", "1"), cex.axis = cex.yaxislab)
axis(1, at = c(10, 26.6, 42.6, 60), labels = c("300", "500", "1000", "10,000"), cex.axis = cex.xaxislab)
abline(0.95, 0)
abline(0.95 - 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
abline(0.95 + 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
points(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], pch = 15, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], lty = "twodash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], pch = 7, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], lty = "longdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], pch = 16, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], lty = "dotted", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], pch = 10, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], lty = "dotdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], pch = 17, cex = cex.points)
# psi10 sc3
data <- dat[which(dat$scenario == "Scenario 3 (Nonregular)"),]
plot(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], type = "l", ylim = c(0.9, 1.00), xlab = "", ylab = expression(paste("Coverage ", psi[10])), xaxt = "n", yaxt = "n", lty = "solid", lwd = cex.lwd, main = "Scenario 3 (Nonregular)", cex.lab = cex.axistitle, cex.main = cex.title)
axis(2, at = c(0.9, 0.95, 1.00), labels = c("0.90", "0.95", "1"), cex.axis = cex.yaxislab)
axis(1, at = c(10, 26.6, 42.6, 60), labels = c("300", "500", "1000", "10,000"), cex.axis = cex.xaxislab)
abline(0.95, 0)
abline(0.95 - 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
abline(0.95 + 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
points(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], pch = 15, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], lty = "twodash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], pch = 7, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], lty = "longdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], pch = 16, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], lty = "dotted", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], pch = 10, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], lty = "dotdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], pch = 17, cex = cex.points)
# psi10 sc4
data <- dat[which(dat$scenario == "Scenario 4 (Near-nonregular)"),]
plot(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], type = "l", ylim = c(0.9, 1.00), xlab = "", ylab = expression(paste("Coverage ", psi[10])), xaxt = "n", yaxt = "n", lty = "solid", lwd = cex.lwd, main = "Scenario 4 (Near-nonregular)", cex.lab = cex.axistitle, cex.main = cex.title)
axis(2, at = c(0.9, 0.95, 1.00), labels = c("0.90", "0.95", "1"), cex.axis = cex.yaxislab)
axis(1, at = c(10, 26.6, 42.6, 60), labels = c("300", "500", "1000", "10,000"), cex.axis = cex.xaxislab)
abline(0.95, 0)
abline(0.95 - 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
abline(0.95 + 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
points(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], pch = 15, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], lty = "twodash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], pch = 7, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], lty = "longdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], pch = 16, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], lty = "dotted", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], pch = 10, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], lty = "dotdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], pch = 17, cex = cex.points)
# psi10 sc5
data <- dat[which(dat$scenario == "Scenario 5 (Nonregular)"),]
plot(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], type = "l", ylim = c(0.9, 1.00), xlab = "", ylab = expression(paste("Coverage ", psi[10])), xaxt = "n", yaxt = "n", lty = "solid", lwd = cex.lwd, main = "Scenario 5 (Nonregular)", cex.lab = cex.axistitle, cex.main = cex.title)
axis(2, at = c(0.9, 0.95, 1.00), labels = c("0.90", "0.95", "1"), cex.axis = cex.yaxislab)
axis(1, at = c(10, 26.6, 42.6, 60), labels = c("300", "500", "1000", "10,000"), cex.axis = cex.xaxislab)
abline(0.95, 0)
abline(0.95 - 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
abline(0.95 + 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
points(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], pch = 15, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], lty = "twodash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], pch = 7, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], lty = "longdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], pch = 16, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], lty = "dotted", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], pch = 10, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], lty = "dotdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], pch = 17, cex = cex.points)
# psi10 sc6
data <- dat[which(dat$scenario == "Scenario 6 (Near-nonregular)"),]
par(mar = c(2.1, 5.1, 3.1, 1.1), xpd = FALSE)
plot(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], type = "l", ylim = c(0.9, 1.00), xlab = "", ylab = expression(paste("Coverage ", psi[10])), xaxt = "n", yaxt = "n", lty = "solid", lwd = cex.lwd, main = "Scenario 6 (Near-nonregular)", cex.lab = cex.axistitle, cex.main = cex.title)
axis(2, at = c(0.9, 0.95, 1.00), labels = c("0.90", "0.95", "1"), cex.axis = cex.yaxislab)
axis(1, at = c(10, 26.6, 42.6, 60), labels = c("300", "500", "1000", "10,000"), cex.axis = cex.xaxislab)
abline(0.95, 0)
abline(0.95 - 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
abline(0.95 + 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
points(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], pch = 15, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], lty = "twodash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], pch = 7, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], lty = "longdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], pch = 16, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], lty = "dotted", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], pch = 10, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], lty = "dotdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], pch = 17, cex = cex.points)
# psi10 sc7
data <- dat[which(dat$scenario == "Scenario 7 (Regular)"),]
par(mar = c(4.1, 5.1, 3.1, 1.1), xpd = FALSE)
plot(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], type = "l", ylim = c(0.9, 1.00), xlab = "Sample size", ylab = expression(paste("Coverage ", psi[10])), xaxt = "n", yaxt = "n", lty = "solid", lwd = cex.lwd, main = "Scenario 7 (Regular)", cex.lab = cex.axistitle, cex.main = cex.title)
axis(2, at = c(0.9, 0.95, 1.00), labels = c("0.90", "0.95", "1"), cex.axis = cex.yaxislab)
axis(1, at = c(10, 26.6, 42.6, 60), labels = c("300", "500", "1000", "10,000"), cex.axis = cex.xaxislab)
abline(0.95, 0)
abline(0.95 - 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
abline(0.95 + 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
points(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], pch = 15, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], lty = "twodash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], pch = 7, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], lty = "longdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], pch = 16, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], lty = "dotted", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], pch = 10, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], lty = "dotdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], pch = 17, cex = cex.points)
# psi10 sc8
data <- dat[which(dat$scenario == "Scenario 8 (Regular)"),]
plot(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], type = "l", ylim = c(0.9, 1.00), xlab = "Sample size", ylab = expression(paste("Coverage ", psi[10])), xaxt = "n", yaxt = "n", lty = "solid", lwd = cex.lwd, main = "Scenario 8 (Regular)", cex.lab = cex.axistitle, cex.main = cex.title)
axis(2, at = c(0.9, 0.95, 1.00), labels = c("0.90", "0.95", "1"), cex.axis = cex.yaxislab)
axis(1, at = c(10, 26.6, 42.6, 60), labels = c("300", "500", "1000", "10,000"), cex.axis = cex.xaxislab)
abline(0.95, 0)
abline(0.95 - 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
abline(0.95 + 1.96*(0.95*0.05/1000)^(1/2), 0, lty = "dashed")
points(x = data$n.coord[which(data$method == "Asymptotic (adjusted)")], y = data$coverage[which(data$method == "Asymptotic (adjusted)")], pch = 15, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], lty = "twodash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Asymptotic (naive)")], y = data$coverage[which(data$method == "Asymptotic (naive)")], pch = 7, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], lty = "longdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Standard bootstrap")], y = data$coverage[which(data$method == "Standard bootstrap")], pch = 16, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], lty = "dotted", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 1")], y = data$coverage[which(data$method == "Parametric bootstrap 1")], pch = 10, cex = cex.points)
lines(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], lty = "dotdash", lwd = cex.lwd)
points(x = data$n.coord[which(data$method == "Parametric bootstrap 2")], y = data$coverage[which(data$method == "Parametric bootstrap 2")], pch = 17, cex = cex.points)
dev.off()

