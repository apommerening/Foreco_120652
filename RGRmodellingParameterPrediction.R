# RGR modelling. Duhnen in August 2021. Updated on 17.11.2021.
rm(list = ls())
# install.packages("Rcpp", dep = T)
# library(spatstat)
library(Rcpp) 
options(digits = 14, width = 100)

# Sys.setenv("LANGUAGE"="EN")

pc <- F
drive <- ""
if(pc == T)
 drive <- "C:"  

sourceFile <- paste(drive, "/Users/arng0001/Dropbox/Gwaith/Cyhoeddi/RGRmodel/", sep = "")
dataFile <- paste(drive, "/Users/arng0001/Dropbox/Gwaith/Cyhoeddi/RGRmodel/Data/", sep = "")
codeFile <- paste(drive, "/Users/arng0001/Dropbox/Gwaith/Cyhoeddi/RGRmodel/Code/", sep = "")
outputFile <- paste(drive, "/Users/arng0001/Dropbox/Gwaith/Cyhoeddi/RGRmodel/Results/", sep = "")

sourceCpp(paste(codeFile, "findNeighbours.cpp", sep = ""))
source(paste(codeFile, "GrowthFunctions.R", sep = ""))

#-- Read in data for regressions --#
#----------------------------------#
# rm(TreeList)
timeSeries <- "giRegresGD.txt" 
TreeList  <- read.table(paste(dataFile, timeSeries, sep = ""), header = T)
names(TreeList)
tapply(TreeList$year, TreeList$plotno, unique)
sort(table(TreeList$plotno))
no <- 410311  
TreeList  <- TreeList[TreeList$plotno == no, ]
range(TreeList$dbh)
range(TreeList$dm)
tapply(TreeList$dm, TreeList$year, unique)
range(TreeList$year)
TreeList[200 : 300, c(2 : 3, 6 : 7, 11 : 14)]
table(TreeList$Species)
table(TreeList$year)
TreeList  <- TreeList[order(TreeList$Treeno, TreeList$year, decreasing = FALSE), ]
TreeList[1 : 100, c(2, 6, 12, 15)]

par(mar = c(2, 5, 0.5, 0.5), mfrow = c(1, 1))
plot(TreeList$dbh, TreeList$RGR, las = 1, ylab = "", xlab = "", 
     cex = .9, col = "black", pch = 16, axes = FALSE)#, xlim = c(0, max(alpha)), ylim = c(0.20, 0.35))#, ylim = c(0, 1)
axis(1, lwd = 2, cex.axis = 1.8) 
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2)

# Data ellipses - slopes of semi-major axes
source(paste(codeFile, "EllipseCode.R", sep = ""))
library(MASS)
# install.packages("robustbase", dep = T)
library(robustbase) 
dm <- tapply(TreeList$dbh, TreeList$year, mean)
years <- as.numeric(row.names(table(TreeList$year)))
slope1 <- slope2 <- slope3 <- cvr <- cvs <- Dyer1 <- Dyer2 <- Dyer3 <- Dyer4 <- c()
# i <- 1
for (i in 1 : length(years)) {
  singlePlot <- TreeList[TreeList$year == years[i],  ]
  singlePlot <- singlePlot[!is.na(singlePlot$RGR),]
  sma.lm1 <- sma(RGR ~ dbh, data = singlePlot, method = "SMA")
  slope1[i] <- coef(sma.lm1)[[2]]
  fit <- fit.ellipse(singlePlot$dbh, singlePlot$RGR)
  slope2[i] <- fit$angle
  elli <- standard.ellipse(singlePlot$dbh, singlePlot$RGR, confs = NULL, steps = 2)
  slope3[i] <- elli$theta
  singlePlot$rdbh <- singlePlot$dbh / mean(singlePlot$dbh)
  sma.lm2 <- sma(RGR ~ rdbh, data = singlePlot, method = "SMA")
  Dyer1[i] <- coef(sma.lm2)[[2]]
  sma.lm3 <- rlm(RGR ~ rdbh, data = singlePlot, method = "MM")
  Dyer2[i] <- coef(sma.lm3)[[2]]
  sma.lm4 <- lmrob(RGR ~ rdbh, data = singlePlot)
  Dyer3[i] <- coef(sma.lm4)[[2]]
  sma.lm5 <- lqs(RGR ~ rdbh, data = singlePlot, method = "lts")
  Dyer4[i] <- coef(sma.lm5)[[2]]
  cvr[i] <- sd(singlePlot$RGR, na.rm = T) / mean(singlePlot$RGR, na.rm = T)
  cvs[i] <- sd(singlePlot$dbh, na.rm = T) / mean(singlePlot$dbh, na.rm = T)
  rm(singlePlot, sma.lm1, sma.lm2, sma.lm3, sma.lm4, sma.lm5, fit, elli)
}


# Calculate mean slopes
test <- data.frame(slope1, slope2, slope3)
meanSlope <- apply(test, 1, mean) 
rm(test)
test <- data.frame(Dyer1, Dyer2, Dyer3, Dyer4)
meanDyer <- apply(test, 1, mean) 
rm(test)

smean <- ksmooth(x = dm, y = meanSlope, kernel = "normal", bandwidth = 3)
smeanDyer <- ksmooth(x = dm, y = meanDyer, kernel = "normal", bandwidth = 3)
sslope1 <- ksmooth(x = dm, y = slope1, kernel = "normal", bandwidth = 3)
sslope2 <- ksmooth(x = dm, y = slope2, kernel = "normal", bandwidth = 3)
sslope3 <- ksmooth(x = dm, y = slope3, kernel = "normal", bandwidth = 3)
sDyer1 <- ksmooth(x = dm, y = Dyer1, kernel = "normal", bandwidth = 3)
sDyer2 <- ksmooth(x = dm, y = Dyer2, kernel = "normal", bandwidth = 3)
sDyer3 <- ksmooth(x = dm, y = Dyer3, kernel = "normal", bandwidth = 3)
sDyer4 <- ksmooth(x = dm, y = Dyer4, kernel = "normal", bandwidth = 3)

par(mar = c(2.0, 3.6, 0.5, 0.5))
plot(sslope1$x, sslope1$y, type = "l", axes = FALSE, ylab = "", xlab = "", main = "", 
     cex = 1.0, lwd = 2, col = "black", ylim = c(-0.4, 0.3), xlim = c(0, 50))
lines(sslope2$x, sslope2$y, lwd = 2, col = "red")
lines(sslope3$x, sslope3$y, lwd = 2, col = "blue")
lines(sDyer1$x, sDyer1$y, lwd = 2, col = "purple")
lines(sDyer2$x, sDyer2$y, lwd = 2, col = "purple", lty = 2)
lines(sDyer3$x, sDyer3$y, lwd = 2, col = "purple", lty = 3)
lines(sDyer4$x, sDyer4$y, lwd = 2, col = "purple", lty = 4)
lines(smeanDyer$x, smeanDyer$y, lwd = 2, lty = 1, col = "red")
lines(smean$x, smean$y, lwd = 2, lty = 1, col = "green")
axis(side = 1, lwd = 2, las = 1, cex.axis = 1.7)
axis(side = 2, lwd = 2, las = 1, cex.axis = 1.7)
abline(h = 0, lwd = 1, lty = 2)
box(lwd = 2)
rm(slope1, slope2, slope3, sslope1, sslope2, sslope3, smean)

# Growth dominance
gcy <- c()
years <- as.numeric(row.names(table(TreeList$year)))
for (j in 1 : length(years)) {
  TreeList2 <- TreeList[TreeList$year == years[j],  ]
  TreeList2$ba <- pi * (TreeList2$dbh / 200)^2
  TreeList2 <- TreeList2[!is.na(TreeList2$AGR.g),]
  TreeList2$AGR.g[TreeList2$AGR.g < 0] <- 0
  TreeList2 <- TreeList2[order(TreeList2$ba, decreasing = FALSE), ]
  cumBA <- cumsum(TreeList2$ba) / sum(TreeList2$ba)
  cumInc <- cumsum(TreeList2$AGR.g) / sum(TreeList2$AGR.g)
  xarea <- 0
  for(k in 2 : length(TreeList2$ba))
    xarea[k] <- (cumBA[k] - cumBA[k - 1]) * ((cumInc[k] - cumInc[k - 1]) / 2 + cumInc[k - 1])
  gcy[j] <- 1 - sum(xarea) / 0.5
  rm(TreeList2, xarea, cumBA, cumInc)
}
gcy

lossG <- function(abdn, gcy, meanDyer) {
  est <- abdn[1] * gcy / (abdn[2] + gcy)^abdn[3] # Hassell
  dev <- (est - meanDyer)^2
  return(sum(dev, na.rm = T))
}  
myPar <- c(0.06, 5.2, 1)
regresG <- optim(myPar, lossG, gcy = gcy, meanDyer = meanDyer, hessian = T, control = list(maxit = 90000, 
              reltol = .Machine$double.eps, temp = 80000, trace = 6, REPORT = 5000)) # maxit = 900000
regresG$par



pdf(file = paste(outputFile, "CorrGrowthDominanceDyer", no, ".pdf", sep = ""))
par(mar = c(2.0, 3.5, 0.5, 0.5))
plot(gcy, meanDyer, axes = FALSE, ylab = "", xlab = "", main = "", 
     cex = 1.0, pch = 16, col = "black", ylim = c(-0.5, 0.3), xlim = c(-0.2, 0.3))
curve(regresG$par[1] * x / (regresG$par[2] + x)^regresG$par[3], # Hassell
      from = min(gcy),  to = max(gcy),  lwd = 2, col = "red", add = TRUE)
axis(side = 1, lwd = 2, las = 1, cex.axis = 1.7)
axis(side = 2, lwd = 2, las = 1, cex.axis = 1.7)
abline(h = 0, lwd = 1, lty = 2)
abline(v = 0, lwd = 1, lty = 2)
box(lwd = 2)
dev.off()

range(gcy)
range(meanDyer)

cor(gcy, meanDyer)


# Grid search for alpha
alpha <- seq(0, 10, 0.05) # 0.01
distexp <- seq(0.5, 7.0, 0.05)
kk <- 30
counter <- 0
corr <- matrix(NA, nrow = length(alpha), ncol = length(distexp))
for (jj in 1 : length(distexp)) {
  for (z in 1 : length(alpha)) {
    TreeList$ui3 <- matrix(NA, nrow = length(TreeList$Treeno), ncol = length(table(TreeList$year)))
    counter = counter + 1
    cat("Simulation: ", counter, " of ", length(alpha) * length(distexp), " alpha ", alpha[z], " distexp ", distexp[jj], "\n")
    for (i in 1 : length(table(TreeList$year))) {
      TreeList1 <- TreeList[TreeList$year == as.numeric(names(table(TreeList$year)[i])),]
      dn <- findNeighbours(TreeList1$xmax[1], TreeList1$ymax[1], TreeList1$x, TreeList1$y, kk, 1) # 1 means periodic boundary conditions
      neighbours <- distances <- matrix(NA, nrow = length(TreeList1$x), ncol = kk)
      for (j in 1 : kk) {
        distances[,j] <- as.vector(dn[[j]]) # 4
        neighbours[,j] <- 1 + as.vector(dn[[kk + j]]) # 4
      }
      ui3 <- calcTrigDiff(TreeList1$dbh, as.matrix(neighbours), as.matrix(distances), kk, alpha[z], distexp[jj]) # alpha[z] # alpha = 0.5 so that exponent = 1
      df2 <- data.frame(Treeno = TreeList1$Treeno, year = TreeList1$year, ui3)
      test <- merge(TreeList, df2, by = c("Treeno", "year"), all.x = T, all.y = F)
      TreeList$ui3[,i] <- test$ui3.y
      rm(TreeList1, neighbours, distances, test, ui3, df2, dn)
    }
    ui3 <- rowSums(TreeList$ui3, na.rm = T)  
    TreeList$ui3 <- ui3
    rm(ui3)
    dummy <- TreeList
    dummy <- dummy[dummy$dbh > 0, ] 
    dummy <- dummy[c("RGR", "ui3")]
    dummy$RGR[is.infinite(dummy$RGR)] <- NA
    dummy <- dummy[dummy$RGR > 0,]
    dummy <- dummy[!is.na(dummy$RGR),]
    dummy <- na.omit(dummy)
    corr[z, jj] <- cor(dummy$ui3, dummy$RGR, use = "complete.obs")
    rm(dummy)
  }
}

xy <- which(abs(corr) == max(abs(corr), na.rm = T), arr.ind = TRUE)
max(corr, na.rm = T)

alphaDummy <- alpha
par(mar = c(2, 6.2, 0.5, 0.5), mfrow = c(1, 1))
plot(alphaDummy, corr[,xy[[2]]], las = 1, ylab = "", xlab = "", 
     cex = .9, col = "black", pch = 16, axes = FALSE)#, xlim = c(0, max(alpha)), ylim = c(0.20, 0.35))#, ylim = c(0, 1)
axis(1, lwd = 2, cex.axis = 1.8) 
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2)

(alpha <- alpha[xy[[1]]])
(exponent <- distexp[xy[[2]]])

# Dominance
TreeList$ui3 <- TreeList$dr <- matrix(NA, nrow = length(TreeList$Treeno), ncol = length(table(TreeList$year)))
# i <- 1
for (i in 1 : length(table(TreeList$year))) {
  TreeList1 <- TreeList[TreeList$year == as.numeric(names(table(TreeList$year)[i])),]
  dr <- TreeList1$dbh / mean(TreeList1$dbh, na.rm = T) # median(TreeList1$dbh, na.rm = T)
  dn <- findNeighbours(TreeList1$xmax[1], TreeList1$ymax[1], TreeList1$x, TreeList1$y, kk, 1) # 1 means periodic boundary conditions
  neighbours <- distances <- matrix(NA, nrow = length(TreeList1$x), ncol = kk)
  for (j in 1 : kk) {
      distances[,j] <- as.vector(dn[[j]]) # 4
      neighbours[,j] <- 1 + as.vector(dn[[kk + j]]) # 4
  }
  ui3 <- calcTrigDiff(TreeList1$dbh, as.matrix(neighbours), as.matrix(distances), kk, alpha, exponent) # alpha[z] # alpha = 0.5 so that exponent = 1
  df2 <- data.frame(Treeno = TreeList1$Treeno, year = TreeList1$year, ui3)
  test <- merge(TreeList, df2, by = c("Treeno", "year"), all.x = T, all.y = F)
  TreeList$ui3[,i] <- test$ui3.y
  rm(test, df2)
  df2 <- data.frame(Treeno = TreeList1$Treeno, year = TreeList1$year, dr)
  test <- merge(TreeList, df2, by = c("Treeno", "year"), all.x = T, all.y = F)
  TreeList$dr[,i] <- test$dr.y
  rm(TreeList1, neighbours, distances, ui3, df2, dr, dn, test)
}
ui3 <- rowSums(TreeList$ui3, na.rm = T)  
TreeList$ui3[1 : 10, 1 : length(table(TreeList$year))]
ui3[1 : 10]
TreeList$ui3 <- ui3
rm(ui3)
TreeList[1 : 100, c(2, 6, 12, 15, 19)]
range(TreeList$ui3)
dr <- rowSums(TreeList$dr, na.rm = T)  
TreeList$dr <- dr
rm(dr)
range(TreeList$dr)

par(mar = c(2, 4.1, 0.5, 0.5), mfrow = c(1, 1))
plot(TreeList$dbh, TreeList$ui3, las = 1, ylab = "", xlab = "", 
  cex = .9, col = "black", pch = 19, axes = FALSE, xlim = c(0, ceiling(max(TreeList$dbh, na.rm = T))))#, ylim = c(0, 1))
axis(1, lwd = 2, cex.axis = 1.8) 
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2)

par(mar = c(2, 4.1, 0.5, 0.5), mfrow = c(1, 1))
plot(TreeList$ui3, TreeList$RGR, las = 1, ylab = "", xlab = "", 
  cex = .9, col = "black", pch = 19, axes = FALSE, ylim = c(0, 0.7)) # xlim = c(0, 1),
axis(1, lwd = 2, cex.axis = 1.8) 
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2)

par(mar = c(2, 4.1, 0.5, 0.5), mfrow = c(1, 1))
plot(TreeList$dbh / TreeList$dm, TreeList$RGR, las = 1, ylab = "", xlab = "", 
     cex = .9, col = "black", pch = 19, axes = FALSE, ylim = c(0, 0.7)) # xlim = c(0, 1),
axis(1, lwd = 2, cex.axis = 1.8) 
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2)

par(mar = c(2, 4.1, 0.5, 0.5), mfrow = c(1, 1))
plot(TreeList$ui3, TreeList$dbh / TreeList$dm, las = 1, ylab = "", xlab = "", 
     cex = .9, col = "black", pch = 19, axes = FALSE, ylim = c(0, 1.7), xlim = c(0, 1))
axis(1, lwd = 2, cex.axis = 1.8) 
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2)







TreeList <- TreeList[TreeList$dbh > 0, ] 
TreeList$RGR[is.infinite(TreeList$RGR)] <- NA
TreeList$RGR[TreeList$RGR < 0] <- 0
TreeList <- TreeList[!is.na(TreeList$RGR),]
 

# Regression for median RGR trend
w <- rep(1, length(TreeList$dbh))
lossM <- function(abdn, dbh, RGR, ui3, wx) {
  est <- relativeGrowthFunction(dbh, abdn[1], abdn[2], ui3)
  dev <- 1 / wx * (est - RGR)^2
  return(sum(dev, na.rm = T) / sum(1 / wx, na.rm = T))
}  
myPar <- c(0.1, 1.0)
regresm <- optim(myPar, lossM, dbh = TreeList$dbh, RGR = TreeList$RGR, ui3 = TreeList$ui3, wx = w,
     hessian = T, control = list(maxit = 90000, reltol = .Machine$double.eps, temp = 80000, trace = 6, REPORT = 5000)) # maxit = 900000
options(scipen = 100)
regresm$par
# install.packages("quantreg", dep = T)
library(quantreg)
nlsout <- nlrq(RGR ~ relativeGrowthFunction(dbh, c1, c2, ui3), data = TreeList, start = list(c1 = round(regresm$par[1], 2), c2 =  round(regresm$par[2], 2)), 
               tau = 0.5, trace = T, method = "Nelder-Mead")
summary(nlsout)
param <- c(summary(nlsout)$coefficients[1], summary(nlsout)$coefficients[2])
par(mar = c(2, 3.1, 0.5, 0.5), mfrow = c(1, 1))
plot(TreeList$dbh, TreeList$RGR, las = 1, ylab = "", xlab = "", 
     cex = .9, col = "black", pch = 19, axes = FALSE, ylim = c(0, 2.0)) # xlim = c(0, 1),
curve(relativeGrowthFunction(x, regresm$par[1], regresm$par[2], xcomp = 0.50), from = 0, to = 50, lwd = 2, lty = 1, col = "blue", add = TRUE)
curve(relativeGrowthFunction(x, param[1], param[2], xcomp = 0.50), from = 0, to = 50, lwd = 2, lty = 1, col = "red", add = TRUE)
axis(1, lwd = 2, cex.axis = 1.8) 
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2)
yhat <- relativeGrowthFunction(TreeList$dbh, regresm$par[1], regresm$par[2], TreeList$ui3)
(Eff <- 1 - (sum((TreeList$RGR - yhat)^2, na.rm = T) / sum((TreeList$RGR - mean(TreeList$RGR, na.rm = T))^2, na.rm = T))) # Efficiency
yhat <- relativeGrowthFunction(TreeList$dbh, param[1], param[2], TreeList$ui3)
(Eff <- 1 - (sum((TreeList$RGR - yhat)^2, na.rm = T) / sum((TreeList$RGR - mean(TreeList$RGR, na.rm = T))^2, na.rm = T))) # Efficiency
error.hats <- TreeList$RGR - yhat
group <- as.numeric(cut(TreeList$dbh, 50))
y <- as.numeric(tapply(abs(error.hats), group, var, na.rm = T))
x <- as.numeric(tapply(TreeList$dbh, group, median, na.rm = T))
w <- group
w <- y[w]

# Testing global c1/c2 parameters
par(mar = c(2, 3.1, 0.5, 0.5), mfrow = c(1, 1))
curve(relativeGrowthFunction(x, regresm$par[1], regresm$par[2], xcomp = 0.5), from = 0, 50, lwd = 2, lty = 1, col = "red", 
      axes = FALSE, ylab = "", xlab = "", xlim = c(0, 60), ylim = c(0, 1.1))
curve(relativeGrowthFunction(x, regresm$par[1], regresm$par[2] / 30, xcomp = 0.5), from = 0, to = 50, lwd = 2, lty = 1, col = "blue", add = TRUE)
curve(relativeGrowthFunction(x, regresm$par[1], regresm$par[2] / 40, xcomp = 0.5), from = 0, to = 50, lwd = 2, lty = 1, col = "black", add = TRUE)
legend("topright", legend = c("0", "1 / 30", "1 / 40"),
       col = c("red", "blue", "black"), lty = 1, lwd = 2, cex = 1.8, bty = "n")  
axis(1, lwd = 2, cex.axis = 1.8) 
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2)

par(mar = c(2, 3.1, 0.5, 0.5), mfrow = c(1, 1))
curve(relativeGrowthFunction(x, regresm$par[1], regresm$par[2], xcomp = 0.5), from = 0, 50, lwd = 2, lty = 1, col = "red", 
      axes = FALSE, ylab = "", xlab = "", xlim = c(0, 60), ylim = c(0, 1.1))
curve(relativeGrowthFunction(x, regresm$par[1] * 1.5, regresm$par[2], xcomp = 0.5), from = 0, to = 50, lwd = 2, lty = 1, col = "blue", add = TRUE)
curve(relativeGrowthFunction(x, regresm$par[1] / 1.5, regresm$par[2], xcomp = 0.5), from = 0, to = 50, lwd = 2, lty = 1, col = "black", add = TRUE)
legend("topright", legend = c("0", "* 1.5", "/ 1.5"),
       col = c("red", "blue", "black"), lty = 1, lwd = 2, cex = 1.8, bty = "n")  
axis(1, lwd = 2, cex.axis = 1.8) 
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2)

# Individual-based regression: Starting weights (Run before first LS regression)
w <- rep(1, length(TreeList$dbh))
# Weights for weighted regression (Run before second LS regression and modify number in next line if necessary)
group <- as.numeric(cut(TreeList$dbh, 50))
y <- as.numeric(tapply(abs(error.hats), group, var, na.rm = T))
x <- as.numeric(tapply(TreeList$dbh, group, median, na.rm = T))
w <- group
w <- y[w]

# LS estimation
par0 <- c(0.05, 1.4, 1.05, 1.4)# , 0.3) 
myFact <- 0 # 1E-3 # 0 
myDenom <- 0
abdnD <- optim(par0, lossD, dbh = TreeList$dbh, period = TreeList$period, RGR = TreeList$RGR, dend = TreeList$d.end, xcomp = TreeList$ui3, xdr = TreeList$dbh / TreeList$dm, varr = w, # xdr = TreeList$rank, TreeList$dbh / TreeList$dm
               fact = myFact, denom = myDenom, hessian = T, control = list(maxit = 9000, reltol = .Machine$double.eps, temp = 80000, trace = 6, REPORT = 5000)) # maxit = 900000
options(scipen = 100)
abdnD$par
4 + any(abdnD$par[1 : 4] < 0)

yhat <- estimateGrowth(y = TreeList$dbh, abdn = abdnD$par, xcomp = TreeList$ui3, xdr = TreeList$dbh / TreeList$dm) # # xdr = TreeList$rank, TreeList$dbh / TreeList$dm
varres <- var(TreeList$RGR - yhat, na.rm = TRUE)
(bias <- mean(yhat - TreeList$RGR, na.rm = TRUE)) 
bias / mean(TreeList$RGR, na.rm = TRUE)
(rmse <- sqrt(varres + bias^2))
rmse / mean(TreeList$RGR, na.rm = TRUE)
error.hats <- TreeList$RGR - yhat
(Eff <- 1 - (sum((TreeList$RGR - yhat)^2, na.rm = T) / sum((TreeList$RGR - mean(TreeList$RGR, na.rm = T))^2, na.rm = T))) # Efficiency

(se <- sqrt(sum(error.hats^2, na.rm = T)/(length(error.hats) - 4)))

# c1
myC1 <- estimateC1(abdnD$par[1 : 2], TreeList$ui3, incRegres$dr)
range(myC1)
# c2
myC2 <- estimateC2(abdnD$par[3 : 4], TreeList$dbh, incRegres$dr)
range(myC2)

par(mar = c(2, 5, 0.5, 0.5), mfrow = c(1, 1))
plot(yhat, TreeList$RGR, ylab = "", xlab = "", cex = .9, col = "black", pch = 16, axes = FALSE)
axis(1, cex.axis = 1.8)
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
abline(0, 1)
box(lwd = 2)
  
par(mar = c(2, 5.3, 0.5, 0.5), mfrow = c(1, 1))
plot(TreeList$dbh, error.hats, ylab = "", xlab = "", cex = .9, col = "black", pch = 16, axes = FALSE) #, ylim = c(0, 0.3), xlim = c(0, 0.3))
axis(1, cex.axis = 1.8)
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
abline(h = 0, col = "red")
box(lwd = 2)
  
par(mar = c(2, 2.5, 0.5, 0.5), mfrow = c(1, 1))
plot(TreeList$dbh, error.hats / sd(TreeList$RGR, na.rm = T), ylab = "", xlab = "", cex = .9, col = "black", pch = 16, axes = FALSE) #, ylim = c(0, 0.3), xlim = c(0, 0.3))
axis(1, cex.axis = 1.8)
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
abline(h = 0, col = "red")
box(lwd = 2)

qqnorm(error.hats, las = 1, cex.axis = 0.9)
qqline(error.hats)
library(car)
qqPlot(error.hats)
  
# Influence of competition
par(mar = c(2, 3.1, 0.5, 0.5), mfrow = c(1, 1))
curve(relativeGrowthFunction(x, estimateC1(myPar = abdnD$par[1 : 2], indepVar1 = 0.25, indepVar2 = 0), estimateC2(myPar = abdnD$par[3 : 4], indepVar1 = x, indepVar2 = 0), xcomp = 0.25), from = 0, 50, lwd = 2, lty = 1, col = "red", 
      axes = FALSE, ylab = "", xlab = "", xlim = c(0, 60), ylim = c(0, 1.1))
curve(relativeGrowthFunction(x, estimateC1(myPar = abdnD$par[1 : 2], indepVar1 = 0.50, indepVar2 = 0), estimateC2(myPar = abdnD$par[3 : 4], indepVar1 = x, indepVar2 = 0), xcomp = 0.50), from = 0, to = 50, lwd = 2, lty = 1, col = "blue", add = TRUE)
curve(relativeGrowthFunction(x, estimateC1(myPar = abdnD$par[1 : 2], indepVar1 = 0.70, indepVar2 = 0), estimateC2(myPar = abdnD$par[3 : 4], indepVar1 = x, indepVar2 = 0), xcomp = 0.70), from = 0, to = 50, lwd = 2, lty = 1, col = "black", add = TRUE)
legend("topright", legend = c("0.25", "0.50", "0.70"),
       col = c("red", "blue", "black"), lty = 1, lwd = 2, cex = 1.8, bty = "n")  
axis(1, lwd = 2, cex.axis = 1.8) 
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2)

# With d / dm
par(mar = c(2, 3.1, 0.5, 0.5), mfrow = c(1, 1))
curve(relativeGrowthFunction(x, estimateC1(myPar = abdnD$par[1 : 2], indepVar1 = 0.25, indepVar2 = 0), estimateC2(myPar = abdnD$par[3 : 4], indepVar1 = x, indepVar2 = 0), xcomp = 0.25), from = 0, 50, lwd = 2, lty = 1, col = "red", 
      axes = FALSE, ylab = "", xlab = "", xlim = c(0, 60), ylim = c(0, 1.1))
curve(relativeGrowthFunction(x, estimateC1(myPar = abdnD$par[1 : 2], indepVar1 = 0.75, indepVar2 = 0), estimateC2(myPar = abdnD$par[3 : 4], indepVar1 = x, indepVar2 = 0), xcomp = 0.75), from = 0, to = 50, lwd = 2, lty = 1, col = "blue", add = TRUE)
curve(relativeGrowthFunction(x, estimateC1(myPar = abdnD$par[1 : 2], indepVar1 = 1.00, indepVar2 = 0), estimateC2(myPar = abdnD$par[3 : 4], indepVar1 = x, indepVar2 = 0), xcomp = 1.50), from = 0, to = 50, lwd = 2, lty = 1, col = "black", add = TRUE)
legend("topright", legend = c("0.25", "0.75", "1.50"),
       col = c("red", "blue", "black"), lty = 1, lwd = 2, cex = 1.8, bty = "n")  
axis(1, lwd = 2, cex.axis = 1.8) 
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2)

# Influence of size
par(mar = c(2, 3.1, 0.5, 0.5), mfrow = c(1, 1))
curve(relativeGrowthFunction(20, estimateC1(myPar = abdnD$par[1 : 2], indepVar1 = x, indepVar2 = 0), estimateC2(myPar = abdnD$par[3 : 4], indepVar1 = 20, indepVar2 = 0), xcomp = x), from = 0, 1, lwd = 2, lty = 1, col = "red", 
      axes = FALSE, ylab = "", xlab = "", xlim = c(0, 1), ylim = c(0, 0.5))
curve(relativeGrowthFunction(50, estimateC1(myPar = abdnD$par[1 : 2], indepVar1 = x, indepVar2 = 0), estimateC2(myPar = abdnD$par[3 : 4], indepVar1 = 50, indepVar2 = 0), xcomp = x), from = 0, to = 1, lwd = 2, lty = 1, col = "blue", add = TRUE)
curve(relativeGrowthFunction(80, estimateC1(myPar = abdnD$par[1 : 2], indepVar1 = x, indepVar2 = 0), estimateC2(myPar = abdnD$par[3 : 4], indepVar1 = 80, indepVar2 = 0), xcomp = x), from = 0, to = 1, lwd = 2, lty = 1, col = "black", add = TRUE)
legend("topright", legend = c("20 cm", "50 cm", "80 cm"),
       col = c("red", "blue", "black"), lty = 1, lwd = 2, cex = 1.8, bty = "n")  
axis(1, lwd = 2, cex.axis = 1.8) 
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2)

# With d / dm
par(mar = c(2, 3.1, 0.5, 0.5), mfrow = c(1, 1))
curve(relativeGrowthFunction(20, estimateC1(myPar = abdnD$par[1 : 2], indepVar1 = x / 1.1, indepVar2 = 0), estimateC2(myPar = abdnD$par[3 : 4], indepVar1 = 20, indepVar2 = 0), xcomp = x), from = 0, 1.5, lwd = 2, lty = 1, col = "red", 
      axes = FALSE, ylab = "", xlab = "", xlim = c(0, 1.5), ylim = c(0, 0.5))
curve(relativeGrowthFunction(50, estimateC1(myPar = abdnD$par[1 : 2], indepVar1 = x / 1.1, indepVar2 = 0), estimateC2(myPar = abdnD$par[3 : 4], indepVar1 = 50, indepVar2 = 0), xcomp = x), from = 0, to = 1.5, lwd = 2, lty = 1, col = "blue", add = TRUE)
curve(relativeGrowthFunction(80, estimateC1(myPar = abdnD$par[1 : 2], indepVar1 = x / 1.1, indepVar2 = 0), estimateC2(myPar = abdnD$par[3 : 4], indepVar1 = 80, indepVar2 = 0), xcomp = x), from = 0, to = 1.5, lwd = 2, lty = 1, col = "black", add = TRUE)
legend("topright", legend = c("20 cm", "50 cm", "80 cm"),
       col = c("red", "blue", "black"), lty = 1, lwd = 2, cex = 1.8, bty = "n")  
axis(1, lwd = 2, cex.axis = 1.8) 
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2)

# c1
par(mar = c(2, 4.5, 0.5, 0.5), mfrow = c(1, 1))
curve(estimateC1(abdnD$par[1 : 2], x), from = min(TreeList$ui3, na.rm = T), to =  max(TreeList$ui3, na.rm = T), 
  lwd = 2, ylab = "", xlab = "", col = "black", axes = FALSE)
axis(1, lwd = 2, cex.axis = 1.8) 
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2)
  
# c2
par(mar = c(2, 4.5, 0.5, 0.5), mfrow = c(1, 1))
curve(estimateC2(abdnD$par[3 : 4], x), from = min(TreeList$dbh, na.rm = T), to =  max(TreeList$dbh, na.rm = T), 
    lwd = 2, ylab = "", xlab = "", col = "black", axes = FALSE)
axis(1, lwd = 2, cex.axis = 1.8) 
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2)
  



# Maximum likelihood
par0 <- c(par0, 5) # LS start parameters
abdn.ML <- optim(par0, loss.ML, dbh = TreeList$dbh, period = TreeList$period, RGR = TreeList$RGR, dend = TreeList$d.end, xcomp = TreeList$ui3, xdr = TreeList$dbh / TreeList$dm, # xdr = TreeList$rank, TreeList$dbh / TreeList$dm
                 fact = myFact, denom = myDenom, hessian = T, control = list(maxit = 900, reltol = .Machine$double.eps, temp = 80000, trace = 6, fnscale = -1, REPORT = 5000))

# Comparison of LS and ML
abdn.ML$par
abdnD$par
# Evaluate the L2 loss function at the ML location (ignoring variance)
lossD(abdn.ML$par[-length(abdn.ML$par)], dbh = TreeList$dbh, period = TreeList$period, RGR = TreeList$RGR, dend = TreeList$d.end, xcomp = TreeList$ui3, xdr = TreeList$dbh / TreeList$dm, varr = 1,
      fact = myFact, denom = myDenom) # TreeList$rank, TreeList$dbh / TreeList$dm
abdnD$value
  
yhat <- estimateGrowth(y = TreeList$dbh, abdn = abdn.ML$par, xcomp = TreeList$ui3, xdr = TreeList$dbh / TreeList$dm) # xdr = TreeList$rank, TreeList$dbh / TreeList$dm
varres <- var(TreeList$RGR - yhat, na.rm = TRUE)
(bias <- mean(yhat - TreeList$RGR, na.rm = TRUE)) 
bias / mean(TreeList$RGR, na.rm = TRUE)
(rmse <- sqrt(varres + bias^2))
rmse / mean(TreeList$RGR, na.rm = TRUE)
error.hats <- TreeList$RGR - yhat
(Eff <- 1 - (sum((TreeList$RGR - yhat)^2, na.rm = T) / sum((TreeList$RGR - mean(TreeList$RGR, na.rm = T))^2, na.rm = T))) # Efficiency
  
(se <- sqrt(sum(error.hats^2, na.rm = T)/(length(error.hats) - 4)))

# c1
myC1 <- estimateC1(abdn.ML$par[1 : 2], TreeList$ui3, incRegres$dr)
range(myC1)
# c2
myC2 <- estimateC2(abdn.ML$par[3 : 4], TreeList$dbh, incRegres$dr)
range(myC2)

par(mar = c(2, 5, 0.5, 0.5), mfrow = c(1, 1))
plot(yhat, TreeList$RGR, ylab = "", xlab = "", cex = .9, col = "black", pch = 16, axes = FALSE)
axis(1, cex.axis = 1.8)
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
abline(0, 1)
box(lwd = 2)
  
par(mar = c(2, 5.3, 0.5, 0.5), mfrow = c(1, 1))
plot(TreeList$dbh, error.hats, ylab = "", xlab = "", cex = .9, col = "black", pch = 16, axes = FALSE) #, ylim = c(0, 0.3), xlim = c(0, 0.3))
axis(1, cex.axis = 1.8)
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
abline(h = 0, col = "red")
box(lwd = 2)
  
par(mar = c(2, 2.5, 0.5, 0.5), mfrow = c(1, 1))
plot(TreeList$dbh, error.hats / sd(TreeList$RGR), ylab = "", xlab = "", cex = .9, col = "black", pch = 16, axes = FALSE) #, ylim = c(0, 0.3), xlim = c(0, 0.3))
axis(1, cex.axis = 1.8)
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
abline(h = 0, col = "red")
box(lwd = 2)

# c2
par(mar = c(2, 4.5, 0.5, 0.5), mfrow = c(1, 1))
curve(estimateC2(abdn.ML$par[3 : 4], x), from = min(TreeList$dbh, na.rm = T), to =  max(TreeList$dbh, na.rm = T), 
      lwd = 2, ylab = "", xlab = "", col = "black", axes = FALSE)
axis(1, lwd = 2, cex.axis = 1.8) 
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2)

