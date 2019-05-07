##############################################################
##                                                          ##
##          REPRODUCE FIGURES 1-3 IN PAPER                  ##
##                                                          ##
##----------------------------------------------------------##  
##                                                          ## 
## This code reproduces the figures in the paper using the  ##
## precompiled results located in subfolder "Results'       ##
##                                                          ## 
## Input                                                    ##    
##   - inpath: path containing results .txt files           ##
##   - outpath: destination path for the figures            ##
##                                                          ##
## Output                                                   ##
##   - 3 PNG figures                                        ##
##                                                          ## 
## Date: February 12, 2019                                  ##
##############################################################
rm(list = ls())
# inpath <- 
# outpath <- 

#### --------- Figure 1: coverage for psi11, unknown error distribution  --------- ####
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
keep <- seq(2,10,by=2)
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
png(paste(outpath, "/Figure1.png", sep = ""), width = 1136, height = 434)
mat <- matrix(c(2,1,3,1), nrow = 2)
layout(mat, c(1,1), c(8,1))
par(xpd = TRUE, mar = c(0,0,0,0))
plot(0,type='n',axes=FALSE,ann=FALSE, xlim = c(0,10))
legend(0,1,bty='n',c("Asymptotic (adjusted)"), lty = c("solid"), pch = c(15), cex = cex.axistitle-0.1, lwd = cex.lwd, pt.cex = cex.points, seg.len = 3)
legend(2.15,1,bty='n',c("Asymptotic (naive)"), lty = c("twodash"), pch = c(7), cex = cex.axistitle-0.1, lwd = cex.lwd, pt.cex = cex.points, seg.len = 3)
legend(4.05,1,bty='n',c("Standard bootstrap"), lty = c("longdash"), pch = c(16), cex = cex.axistitle-0.1, lwd = cex.lwd, pt.cex = cex.points, seg.len = 3)
legend(6.05,1,bty='n',c("Parametric bootstrap 1"), lty = c("dotted"), pch = c(10), cex = cex.axistitle-0.1, lwd = cex.lwd, pt.cex = cex.points, seg.len = 3)
legend(8.3,1,bty='n',c("Parametric bootstrap 2"), lty = c("dotdash"), pch = c(17), cex = cex.axistitle-0.1, lwd = cex.lwd, pt.cex = cex.points, seg.len = 3)
par(mar = c(4.1, 4.1, 3.1, 1.1), xpd = FALSE)
# normal
plot(x = data$n.coord[which(data$method == "Asymptotic (adjusted)" & data$error == "Normal")], y = data$coverage[which(data$method == "Asymptotic (adjusted)" & data$error == "Normal")], type = "l", ylim = c(0.9, 1.00), xlab = "Sample size", ylab = "Coverage", xaxt = "n", yaxt = "n", lty = "solid", lwd = cex.lwd, main = "Log-normal survival times", cex.lab = cex.axistitle, cex.main = cex.title)
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
plot(x = data$n.coord[which(data$method == "Asymptotic (adjusted)" & data$error == "Weibull")], y = data$coverage[which(data$method == "Asymptotic (adjusted)" & data$error == "Weibull")], type = "l", ylim = c(0.9, 1.00), xlab = "Sample size", ylab = "Coverage", xaxt = "n", yaxt = "n", lty = "solid", lwd = cex.lwd, main = "Weibull survival times", cex.lab = cex.axistitle, cex.main = cex.title)
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


#### --------- Figure 2: coverage for psi11, model misspecification  --------- ####
#import data
loutcome1000 <- read.table(paste(inpath,"/misspec_tf1_1000_cov.txt", sep = ""))
lt1rmt1000 <- read.table(paste(inpath,"/misspec_trmt1_1000_cov.txt", sep = ""))
lc1ens1000 <- read.table(paste(inpath,"/misspec_cens1_1000_cov.txt", sep = ""))
loutcome100 <- read.table(paste(inpath,"/misspec_tf1_100_cov.txt", sep = ""))
lt1rmt100 <- read.table(paste(inpath,"/misspec_trmt1_100_cov.txt", sep = ""))
lc1ens100 <- read.table(paste(inpath,"/misspec_cens1_100_cov.txt", sep = ""))
keep <- seq(2,10,by=2)
o1000 <- apply(loutcome1000[,keep], 2, mean)
t1000 <- apply(lt1rmt1000[,keep], 2, mean)
c1000 <- apply(lc1ens1000[,keep], 2, mean)
o100 <- apply(loutcome100[,keep], 2, mean)
t100 <- apply(lt1rmt100[,keep], 2, mean)
c100 <- apply(lc1ens100[,keep], 2, mean)
# prepare data
data <- expand.grid(method = c("Asymptotic (adjusted)", "Asymptotic (naive)", "Standard bootstrap", "Parametric bootstrap 1", "Parametric bootstrap 2"), model = c("Outcome", "Treatment", "Censoring"), n = c(100, 1000))
g <- 0.15
data$y.coord <- c(rep(c(5-g,4-g,3-g,2-g,1-g), 3), rep(c(5+g,4+g,3+g,2+g,1+g), 3))
data$coverage <- c(o100,t100,c100,o1000,t1000,c1000)
lb <- 0.95 - 1.96*(0.95*0.05/1000)^(1/2)
ub <- 0.95 + 1.96*(0.95*0.05/1000)^(1/2)
# prepare plot
cex.points <- 3
cex.lwd <- 1.9
cex.xaxislab <- 1.5 
cex.yaxislab <- 1.8
cex.title <- 2
cex.axistitle <- 1.8
png(paste(outpath, "/Figure2.png", sep = ""), width = 1000, height = 350)
mat <- matrix(c(1,2,3), nrow = 1)
layout(mat, c(2.3,1,1), c(1))
par(mar = c(5.1,16.1,3.1,1.1), xpd = TRUE)
# treatment-free misspecified
dat <- data[which(data$model == "Outcome" & data$n == 100),]
plot(x = dat$coverage[which(dat$method == "Asymptotic (adjusted)")], y = dat$y.coord[which(dat$method == "Asymptotic (adjusted)")], xlim = c(0.84,0.97), ylim = c(0.8,5.2), ylab = "", yaxt='n', xaxt = 'n', xlab = "Coverage", main = "Treatment-free misspecified", cex.main = cex.title, cex.lab = cex.axistitle)
text(x = 0.835, y = c(5,4,3,2,1), pos = 2,labels = c("Asymptotic (adjusted)", "Asymptotic (naive)", "Standard bootstrap", "Parametric bootstrap 1", "Parametric bootstrap 2"), cex = cex.yaxislab)
axis(1, at = c(0.85,0.90,0.95), cex.axis = cex.xaxislab)
lines(x = c(0.835,0.975), y = c(5-g,5-g), col = "azure2", lwd = 2)
lines(x = c(0.835,0.975), y = c(4-g,4-g), col = "azure2", lwd = 2)
lines(x = c(0.835,0.975), y = c(3-g,3-g), col = "azure2", lwd = 2)
lines(x = c(0.835,0.975), y = c(2-g,2-g), col = "azure2", lwd = 2)
lines(x = c(0.835,0.975), y = c(1-g,1-g), col = "azure2", lwd = 2)
lines(x = c(0.835,0.975), y = c(5+g,5+g), col = "azure2", lwd = 2)
lines(x = c(0.835,0.975), y = c(4+g,4+g), col = "azure2", lwd = 2)
lines(x = c(0.835,0.975), y = c(3+g,3+g), col = "azure2", lwd = 2)
lines(x = c(0.835,0.975), y = c(2+g,2+g), col = "azure2", lwd = 2)
lines(x = c(0.835,0.975), y = c(1+g,1+g), col = "azure2", lwd = 2)
points(x = dat$coverage, y = dat$y.coord, pch = 7, cex = cex.points)
dat <- data[which(data$model == "Outcome" & data$n == 1000),]
points(x = dat$coverage, y = dat$y.coord, pch = 15, cex = cex.points)
lines(x = rep(lb,2), y = c(0.64,5.36), lty = "dashed")
lines(x = rep(ub,2), y = c(0.64,5.36), lty = "dashed")
lines(x = rep(0.95,2), y = c(0.64,5.36))
box(lwd = 2.5)
# treatment misspecified
par(mar = c(4.1,1.1,3.1,1.1), xpd = FALSE)
dat <- data[which(data$model == "Treatment" & data$n == 100),]
plot(x = dat$coverage[which(dat$method == "Asymptotic (adjusted)")], y = dat$y.coord[which(dat$method == "Asymptotic (adjusted)")], xlim = c(0.91,0.97), ylim = c(0.8,5.2), ylab = "", yaxt='n', xaxt = 'n', xlab = "Coverage", main = "Treatment misspecified", cex.main = cex.title, cex.lab = cex.axistitle)
axis(1, at = c(0.91,0.93,0.95,0.97), cex.axis = cex.xaxislab)
lines(x = c(0.9,0.975), y = c(5-g,5-g), col = "azure2", lwd = 2)
lines(x = c(0.9,0.975), y = c(4-g,4-g), col = "azure2", lwd = 2)
lines(x = c(0.9,0.975), y = c(3-g,3-g), col = "azure2", lwd = 2)
lines(x = c(0.9,0.975), y = c(2-g,2-g), col = "azure2", lwd = 2)
lines(x = c(0.9,0.975), y = c(1-g,1-g), col = "azure2", lwd = 2)
lines(x = c(0.9,0.975), y = c(5+g,5+g), col = "azure2", lwd = 2)
lines(x = c(0.9,0.975), y = c(4+g,4+g), col = "azure2", lwd = 2)
lines(x = c(0.9,0.975), y = c(3+g,3+g), col = "azure2", lwd = 2)
lines(x = c(0.9,0.975), y = c(2+g,2+g), col = "azure2", lwd = 2)
lines(x = c(0.9,0.975), y = c(1+g,1+g), col = "azure2", lwd = 2)
points(x = dat$coverage, y = dat$y.coord, pch = 7, cex = cex.points)
dat <- data[which(data$model == "Treatment" & data$n == 1000),]
points(x = dat$coverage, y = dat$y.coord, pch = 15, cex = cex.points)
lines(x = rep(lb,2), y = c(0.64,5.36), lty = "dashed")
lines(x = rep(ub,2), y = c(0.64,5.36), lty = "dashed")
lines(x = rep(0.95,2), y = c(0.64,5.36))
box(lwd = 2.5)
# censoring misspecified
dat <- data[which(data$model == "Censoring" & data$n == 100),]
plot(x = dat$coverage[which(dat$method == "Asymptotic (adjusted)")], y = dat$y.coord[which(dat$method == "Asymptotic (adjusted)")], xlim = c(0.91,0.97), ylim = c(0.8,5.2), ylab = "", yaxt='n', xaxt = 'n', xlab = "Coverage", main = "Censoring misspecified", cex.main = cex.title, cex.lab = cex.axistitle)
axis(1, at = c(0.91,0.93,0.95,0.97), cex.axis = cex.xaxislab)
lines(x = c(0.9,0.975), y = c(5-g,5-g), col = "azure2", lwd = 2)
lines(x = c(0.9,0.975), y = c(4-g,4-g), col = "azure2", lwd = 2)
lines(x = c(0.9,0.975), y = c(3-g,3-g), col = "azure2", lwd = 2)
lines(x = c(0.9,0.975), y = c(2-g,2-g), col = "azure2", lwd = 2)
lines(x = c(0.9,0.975), y = c(1-g,1-g), col = "azure2", lwd = 2)
lines(x = c(0.9,0.975), y = c(5+g,5+g), col = "azure2", lwd = 2)
lines(x = c(0.9,0.975), y = c(4+g,4+g), col = "azure2", lwd = 2)
lines(x = c(0.9,0.975), y = c(3+g,3+g), col = "azure2", lwd = 2)
lines(x = c(0.9,0.975), y = c(2+g,2+g), col = "azure2", lwd = 2)
lines(x = c(0.9,0.975), y = c(1+g,1+g), col = "azure2", lwd = 2)
points(x = dat$coverage, y = dat$y.coord, pch = 7, cex = cex.points)
dat <- data[which(data$model == "Censoring" & data$n == 1000),]
points(x = dat$coverage, y = dat$y.coord, pch = 15, cex = cex.points)
lines(x = rep(lb,2), y = c(0.64,5.36), lty = "dashed")
lines(x = rep(ub,2), y = c(0.64,5.36), lty = "dashed")
lines(x = rep(0.95,2), y = c(0.64,5.36))
box(lwd = 2.5)
dev.off()


#### --------- Figure 3: coverage for psi11, nonregulatrity  --------- ####
# import and prepare data
g <- 0.15
keep <- seq(2,10,by=2)
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
method <- factor(rep(c("Asymptotic (adjusted)", "Asymptotic (naive)", "Standard bootstrap", "Parametric bootstrap 1", "Parametric bootstrap 2"), 8), levels = c("Asymptotic", "Asymptotic (naive)", "Standard bootstrap", "Parametric bootstrap 1", "Parametric bootstrap 2"))
scenario <- factor(rep(c("Scenario 1 (Nonregular)", "Scenario 2 (Near-nonregular)", "Scenario 3 (Nonregular)", "Scenario 4 (Near-nonregular)", "Scenario 5 (Nonregular)", "Scenario 6 (Near-nonregular)", "Scenario 7 (Regular)", "Scenario 8 (Regular)"), each = 5), levels = c("Scenario 1 (Nonregular)", "Scenario 2 (Near-nonregular)", "Scenario 3 (Nonregular)", "Scenario 4 (Near-nonregular)", "Scenario 5 (Nonregular)", "Scenario 6 (Near-nonregular)", "Scenario 7 (Regular)", "Scenario 8 (Regular)"))
data1000 <- data.frame(coverage, method, scenario)
data1000$n <- rep(1000, nrow(data1000))
data1000$y.coord <- rep(c(5+g,4+g,3+g,2+g,1+g), 8)
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
method <- factor(rep(c("Asymptotic (adjusted)", "Asymptotic (naive)", "Standard bootstrap", "Parametric bootstrap 1", "Parametric bootstrap 2"), 8), levels = c("Asymptotic", "Asymptotic (naive)", "Standard bootstrap", "Parametric bootstrap 1", "Parametric bootstrap 2"))
scenario <- factor(rep(c("Scenario 1 (Nonregular)", "Scenario 2 (Near-nonregular)", "Scenario 3 (Nonregular)", "Scenario 4 (Near-nonregular)", "Scenario 5 (Nonregular)", "Scenario 6 (Near-nonregular)", "Scenario 7 (Regular)", "Scenario 8 (Regular)"), each = 5), levels = c("Scenario 1 (Nonregular)", "Scenario 2 (Near-nonregular)", "Scenario 3 (Nonregular)", "Scenario 4 (Near-nonregular)", "Scenario 5 (Nonregular)", "Scenario 6 (Near-nonregular)", "Scenario 7 (Regular)", "Scenario 8 (Regular)"))
data300 <- data.frame(coverage, method, scenario)
data300$n <- rep(300, nrow(data300))
data300$y.coord <- rep(c(5-g,4-g,3-g,2-g,1-g), 8)
data <- rbind(data300, data1000)
#prepare plot
lb <- 0.95 - 1.96*(0.95*0.05/1000)^(1/2)
ub <- 0.95 + 1.96*(0.95*0.05/1000)^(1/2)
cex.points <- 3.5
cex.lwd <- 1.9
cex.xaxislab <- 1.7
cex.yaxislab <- 1.8
cex.title <- 2
cex.axistitle <- 1.8
f3 <- 0.5
png(paste(outpath, "/Figure3.png", sep = ""), width = 1345, height = 577)
mat <- matrix(seq(1,8), nrow = 2, ncol = 4)
layout(mat, c(1.52,1,1,1), c(1,1))
par(mar = c(4.1,16.1,3.1,1.1), xpd = TRUE)
dat <- data[which(data$scenario == "Scenario 1 (Nonregular)" & data$n == 1000),]
plot(x = dat$coverage[which(dat$method == "Asymptotic (adjusted)")], y = dat$y.coord[which(dat$method == "Asymptotic (adjusted)")], xlim = c(0.92,0.98), ylim = c(0.8,5.2), ylab = "", yaxt='n', xaxt = 'n', xlab = "Coverage", cex.lab = cex.axistitle, main = "Scenario 1 (Nonregular)", cex.main = cex.title)
text(x = 0.918, y = c(5,4,3,2,1), pos = 2,labels = c("Asymptotic (adjusted)", "Asymptotic (naive)", "Standard bootstrap", "Parametric bootstrap 1", "Parametric bootstrap 2"), cex = cex.yaxislab)
axis(1, at = c(0.93,0.95,0.97), cex.axis = cex.xaxislab)
lines(x = c(0.917,0.982), y = c(5+g,5+g), col = "azure2", lwd = 2)
lines(x = c(0.917,0.982), y = c(4+g,4+g), col = "azure2", lwd = 2)
lines(x = c(0.917,0.982), y = c(3+g,3+g), col = "azure2", lwd = 2)
lines(x = c(0.917,0.982), y = c(2+g,2+g), col = "azure2", lwd = 2)
lines(x = c(0.917,0.982), y = c(1+g,1+g), col = "azure2", lwd = 2)
lines(x = c(0.917,0.982), y = c(5-g,5-g), col = "azure2", lwd = 2)
lines(x = c(0.917,0.982), y = c(4-g,4-g), col = "azure2", lwd = 2)
lines(x = c(0.917,0.982), y = c(3-g,3-g), col = "azure2", lwd = 2)
lines(x = c(0.917,0.982), y = c(2-g,2-g), col = "azure2", lwd = 2)
lines(x = c(0.917,0.982), y = c(1-g,1-g), col = "azure2", lwd = 2)
points(x = dat$coverage, y = dat$y.coord, pch = 15, cex = cex.points)
dat <- data[which(data$scenario == "Scenario 1 (Nonregular)" & data$n == 300),]
points(x = dat$coverage, y = dat$y.coord, pch = 10, cex = cex.points+f3)
lines(x = rep(lb,2), y = c(0.64,5.36), lty = "dashed")
lines(x = rep(ub,2), y = c(0.64,5.36), lty = "dashed")
lines(x = rep(0.95,2), y = c(0.64,5.36))
box(lwd = 2.5)

dat <- data[which(data$scenario == "Scenario 2 (Near-nonregular)" & data$n == 1000),]
plot(x = dat$coverage[which(dat$method == "Asymptotic (adjusted)")], y = dat$y.coord[which(dat$method == "Asymptotic (adjusted)")], xlim = c(0.92,0.98), ylim = c(0.8,5.2), ylab = "", yaxt='n', xaxt = 'n', xlab = "Coverage", cex.lab = cex.axistitle, main = "Scenario 2 (Near-nonregular)", cex.main = cex.title)
text(x = 0.918, y = c(5,4,3,2,1), pos = 2,labels = c("Asymptotic (adjusted)", "Asymptotic (naive)", "Standard bootstrap", "Parametric bootstrap 1", "Parametric bootstrap 2"), cex = cex.yaxislab)
axis(1, at = c(0.93,0.95,0.97), cex.axis = cex.xaxislab)
lines(x = c(0.917,0.982), y = c(5+g,5+g), col = "azure2", lwd = 2)
lines(x = c(0.917,0.982), y = c(4+g,4+g), col = "azure2", lwd = 2)
lines(x = c(0.917,0.982), y = c(3+g,3+g), col = "azure2", lwd = 2)
lines(x = c(0.917,0.982), y = c(2+g,2+g), col = "azure2", lwd = 2)
lines(x = c(0.917,0.982), y = c(1+g,1+g), col = "azure2", lwd = 2)
lines(x = c(0.917,0.982), y = c(5-g,5-g), col = "azure2", lwd = 2)
lines(x = c(0.917,0.982), y = c(4-g,4-g), col = "azure2", lwd = 2)
lines(x = c(0.917,0.982), y = c(3-g,3-g), col = "azure2", lwd = 2)
lines(x = c(0.917,0.982), y = c(2-g,2-g), col = "azure2", lwd = 2)
lines(x = c(0.917,0.982), y = c(1-g,1-g), col = "azure2", lwd = 2)
points(x = dat$coverage, y = dat$y.coord, pch = 15, cex = cex.points)
dat <- data[which(data$scenario == "Scenario 2 (Near-nonregular)" & data$n == 300),]
points(x = dat$coverage, y = dat$y.coord, pch = 10, cex = cex.points+f3)
lines(x = rep(lb,2), y = c(0.64,5.36), lty = "dashed")
lines(x = rep(ub,2), y = c(0.64,5.36), lty = "dashed")
lines(x = rep(0.95,2), y = c(0.64,5.36))
box(lwd = 2.5)

par(mar = c(4.1,1.1,3.1,1.1), xpd = FALSE)
dat <- data[which(data$scenario == "Scenario 3 (Nonregular)" & data$n == 1000),]
plot(x = dat$coverage[which(dat$method == "Asymptotic (adjusted)")], y = dat$y.coord[which(dat$method == "Asymptotic (adjusted)")], xlim = c(0.92,0.98), ylim = c(0.8,5.2), ylab = "", yaxt='n', xaxt = 'n', xlab = "Coverage", main = "Scenario 3 (Nonregular)", cex.main = cex.title, cex.lab = cex.axistitle)
axis(1, at = c(0.93,0.95,0.97), cex.axis = cex.xaxislab)
lines(x = c(0.91,0.99), y = c(5+g,5+g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(4+g,4+g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(3+g,3+g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(2+g,2+g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(1+g,1+g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(5-g,5-g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(4-g,4-g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(3-g,3-g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(2-g,2-g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(1-g,1-g), col = "azure2", lwd = 2)
points(x = dat$coverage, y = dat$y.coord, pch = 15, cex = cex.points)
dat <- data[which(data$scenario == "Scenario 3 (Nonregular)" & data$n == 300),]
points(x = dat$coverage, y = dat$y.coord, pch = 10, cex = cex.points+f3)
lines(x = rep(lb,2), y = c(0.64,5.36), lty = "dashed")
lines(x = rep(ub,2), y = c(0.64,5.36), lty = "dashed")
lines(x = rep(0.95,2), y = c(0.64,5.36))
box(lwd = 2.5)

dat <- data[which(data$scenario == "Scenario 4 (Near-nonregular)" & data$n == 1000),]
plot(x = dat$coverage[which(dat$method == "Asymptotic (adjusted)")], y = dat$y.coord[which(dat$method == "Asymptotic (adjusted)")], xlim = c(0.92,0.98), ylim = c(0.8,5.2), ylab = "", yaxt='n', xaxt = 'n', xlab = "Coverage", main = "Scenario 4 (Near-nonregular)", cex.main = cex.title, cex.lab = cex.axistitle)
axis(1, at = c(0.93,0.95,0.97), cex.axis = cex.xaxislab)
lines(x = c(0.91,0.99), y = c(5+g,5+g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(4+g,4+g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(3+g,3+g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(2+g,2+g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(1+g,1+g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(5-g,5-g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(4-g,4-g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(3-g,3-g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(2-g,2-g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(1-g,1-g), col = "azure2", lwd = 2)
points(x = dat$coverage, y = dat$y.coord, pch = 15, cex = cex.points)
dat <- data[which(data$scenario == "Scenario 4 (Near-nonregular)" & data$n == 300),]
points(x = dat$coverage, y = dat$y.coord, pch = 10, cex = cex.points+f3)
lines(x = rep(lb,2), y = c(0.64,5.36), lty = "dashed")
lines(x = rep(ub,2), y = c(0.64,5.36), lty = "dashed")
lines(x = rep(0.95,2), y = c(0.64,5.36))
box(lwd = 2.5)

dat <- data[which(data$scenario == "Scenario 5 (Nonregular)" & data$n == 1000),]
plot(x = dat$coverage[which(dat$method == "Asymptotic (adjusted)")], y = dat$y.coord[which(dat$method == "Asymptotic (adjusted)")], xlim = c(0.92,0.98), ylim = c(0.8,5.2), ylab = "", yaxt='n', xaxt = 'n', xlab = "Coverage", main = "Scenario 5 (Nonregular)", cex.main = cex.title, cex.lab = cex.axistitle)
axis(1, at = c(0.93,0.95,0.97), cex.axis = cex.xaxislab)
lines(x = c(0.91,0.99), y = c(5+g,5+g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(4+g,4+g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(3+g,3+g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(2+g,2+g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(1+g,1+g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(5-g,5-g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(4-g,4-g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(3-g,3-g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(2-g,2-g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(1-g,1-g), col = "azure2", lwd = 2)
points(x = dat$coverage, y = dat$y.coord, pch = 15, cex = cex.points)
dat <- data[which(data$scenario == "Scenario 5 (Nonregular)" & data$n == 300),]
points(x = dat$coverage, y = dat$y.coord, pch = 10, cex = cex.points+f3)
lines(x = rep(lb,2), y = c(0.64,5.36), lty = "dashed")
lines(x = rep(ub,2), y = c(0.64,5.36), lty = "dashed")
lines(x = rep(0.95,2), y = c(0.64,5.36))
box(lwd = 2.5)

dat <- data[which(data$scenario == "Scenario 6 (Near-nonregular)" & data$n == 1000),]
plot(x = dat$coverage[which(dat$method == "Asymptotic (adjusted)")], y = dat$y.coord[which(dat$method == "Asymptotic (adjusted)")], xlim = c(0.92,0.98), ylim = c(0.8,5.2), ylab = "", yaxt='n', xaxt = 'n', xlab = "Coverage", main = "Scenario 6 (Near-nonregular)", cex.main = cex.title, cex.lab = cex.axistitle)
axis(1, at = c(0.93,0.95,0.97), cex.axis = cex.xaxislab)
lines(x = c(0.91,0.99), y = c(5+g,5+g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(4+g,4+g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(3+g,3+g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(2+g,2+g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(1+g,1+g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(5-g,5-g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(4-g,4-g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(3-g,3-g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(2-g,2-g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(1-g,1-g), col = "azure2", lwd = 2)
points(x = dat$coverage, y = dat$y.coord, pch = 15, cex = cex.points)
dat <- data[which(data$scenario == "Scenario 6 (Near-nonregular)" & data$n == 300),]
points(x = dat$coverage, y = dat$y.coord, pch = 10, cex = cex.points+f3)
lines(x = rep(lb,2), y = c(0.64,5.36), lty = "dashed")
lines(x = rep(ub,2), y = c(0.64,5.36), lty = "dashed")
lines(x = rep(0.95,2), y = c(0.64,5.36))
box(lwd = 2.5)

dat <- data[which(data$scenario == "Scenario 7 (Regular)" & data$n == 1000),]
plot(x = dat$coverage[which(dat$method == "Asymptotic (adjusted)")], y = dat$y.coord[which(dat$method == "Asymptotic (adjusted)")], xlim = c(0.92,0.98), ylim = c(0.8,5.2), ylab = "", yaxt='n', xaxt = 'n', xlab = "Coverage", main = "Scenario 7 (Regular)", cex.main = cex.title, cex.lab = cex.axistitle)
axis(1, at = c(0.93,0.95,0.97), cex.axis = cex.xaxislab)
lines(x = c(0.91,0.99), y = c(5+g,5+g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(4+g,4+g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(3+g,3+g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(2+g,2+g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(1+g,1+g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(5-g,5-g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(4-g,4-g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(3-g,3-g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(2-g,2-g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(1-g,1-g), col = "azure2", lwd = 2)
points(x = dat$coverage, y = dat$y.coord, pch = 15, cex = cex.points)
dat <- data[which(data$scenario == "Scenario 7 (Regular)" & data$n == 300),]
points(x = dat$coverage, y = dat$y.coord, pch = 10, cex = cex.points+f3)
lines(x = rep(lb,2), y = c(0.64,5.36), lty = "dashed")
lines(x = rep(ub,2), y = c(0.64,5.36), lty = "dashed")
lines(x = rep(0.95,2), y = c(0.64,5.36))
box(lwd = 2.5)

dat <- data[which(data$scenario == "Scenario 8 (Regular)" & data$n == 1000),]
plot(x = dat$coverage[which(dat$method == "Asymptotic (adjusted)")], y = dat$y.coord[which(dat$method == "Asymptotic (adjusted)")], xlim = c(0.92,0.98), ylim = c(0.8,5.2), ylab = "", yaxt='n', xaxt = 'n', xlab = "Coverage", main = "Scenario 8 (Regular)", cex.main = cex.title, cex.lab = cex.axistitle)
axis(1, at = c(0.93,0.95,0.97), cex.axis = cex.xaxislab)
lines(x = c(0.91,0.99), y = c(5+g,5+g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(4+g,4+g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(3+g,3+g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(2+g,2+g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(1+g,1+g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(5-g,5-g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(4-g,4-g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(3-g,3-g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(2-g,2-g), col = "azure2", lwd = 2)
lines(x = c(0.91,0.99), y = c(1-g,1-g), col = "azure2", lwd = 2)
points(x = dat$coverage, y = dat$y.coord, pch = 15, cex = cex.points)
dat <- data[which(data$scenario == "Scenario 8 (Regular)" & data$n == 300),]
points(x = dat$coverage, y = dat$y.coord, pch = 10, cex = cex.points+f3)
lines(x = rep(lb,2), y = c(0.64,5.36), lty = "dashed")
lines(x = rep(ub,2), y = c(0.64,5.36), lty = "dashed")
lines(x = rep(0.95,2), y = c(0.64,5.36))
box(lwd = 2.5)
dev.off()


