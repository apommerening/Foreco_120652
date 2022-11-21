# Spatially explicit RGR model. Umea, 22.10.2021. Updated on 17.11.2021.
rm(list = ls())									
options(digits = 6, width = 50)

# Packages
# library(spatstat)
library(MASS)
library(Rcpp)
# install.packages("truncnorm", dep = T)
library(truncnorm)
# install.packages("robustbase", dep = T)
library(robustbase) 

options(digits = 14, width = 100)

# Sys.setenv("LANGUAGE"="EN")

dataFile <- "/Users/arng0001/Dropbox/Gwaith/Cyhoeddi/RGRmodel/Data/"
codeFile <- "/Users/arng0001/Dropbox/Gwaith/Cyhoeddi/RGRmodel/Code/"
outputFile <- "/Users/arng0001/Dropbox/Gwaith/Cyhoeddi/RGRmodel/Results/"

sourceCpp(paste(codeFile, "findNeighbours.cpp", sep = ""))


calcBasalArea <- function(dbh) 
  return(pi * (dbh / 200)^2)


calcTrigDiff <- function(mark, neighbours, distance, k, alpha, dexp) {
  td <- c()
  for (i in 1 : length(mark)) {
    tsum <- 0
    distsum <- 0
    for (j in 1 : k) {
      w <- calcBasalArea(mark[neighbours[i, j]]) / distance[i, j]^dexp
      tsum <- tsum + mark[i]^(2 * alpha) / (mark[i]^(2 * alpha) + mark[neighbours[i, j]]^(2 * alpha)) * w
      distsum <- distsum + w
    }
    td[i] <- tsum / distsum 
  }
  return(td)
}

relativeGrowthFunction <- function(y, xk, xp, xcomp) {
  return(exp(-xk * y * (1 - exp(-xp * xcomp)))) 
}

# Settings
set.seed(round(runif(1, min = 1, max = 10000)))
replications <- 4499 # 2499 # Number of replications
kk <- 30 # Number of nearest neighbours
silent <- FALSE # Switch on/off messages.

# Data ellipses - slopes of semi-major axes
# install.packages("smatr", dep = T)
# install.packages("siar", dep = T)
source(paste(codeFile, "EllipseCode.R", sep = ""))

# Read and produce start configuration
timeSeries <- "giRegresGD.txt" 
TreeList  <- read.table(paste(dataFile, timeSeries, sep = ""), header = T)
no <- 410311 
TreeList  <- TreeList[TreeList$plotno == no, ]  
TreeList  <- TreeList[order(TreeList$Treeno, TreeList$year, decreasing = FALSE), ]
years <- table(TreeList$year)
myYears <- as.numeric(names(years))
simYears <- as.numeric(names(years[length(years)])) - as.numeric(names(years[1])) # Number of simulation years
xmax <- TreeList$xmax[1] # Observation window defined by xmax 
ymax <- TreeList$ymax[1] # Observation window defined by ymax
StartTreeList  <- TreeList[TreeList$year == as.numeric(names(years[1])), ]
StartTreeList <- StartTreeList[c("Treeno", "Species", "x", "y", "dbh")]
EndTreeList <- TreeList[TreeList$year == myYears[length(years)],  ]
range(EndTreeList$dbh)
endDensity <- 67 # 100 # 80
dmean <-  density(EndTreeList$dbh, n = 32, bw = 8, kernel = "gaussian", from = 0, to = endDensity) 
dmean$x
dmean$y
plot(dmean)
rm(TreeList, timeSeries, years, EndTreeList)

# Read model parameters
source(paste(codeFile, "ModelParametersGD", no, ".R", sep = "")) 
# Containers for results
dm <- meanSlope <- meanDyer <- cvr <- cvs <- outliers <- matrix(nrow = simYears + 1, ncol = replications) 
dend <- list()

for (j in 1 : replications) {
  TreeList <- StartTreeList
  for (i in 1 : (myYears[length(myYears)] - myYears[1] + 1)) {
    # Mortality
    p <- 1 - exp(-exp(parameters$m0 + parameters$m1 * 1 / TreeList$dbh + parameters$m2 * TreeList$prevRGR)) # Gompit
    dying <- runif(length(p), min = 0, max = 1) < p
    if(i > 1) # Assuming that mortality has already happened in year 1
      TreeList <- TreeList[!dying,] # Keep only alive ones
    rm(p, dying)
    # Growth
    dn <- findNeighbours(xmax, ymax, TreeList$x, TreeList$y, kk, 1) # 1 means periodic boundary conditions
    neighbours <- distances <- matrix(NA, nrow = length(TreeList$x), ncol = kk)
    for (k in 1 : kk) {
      distances[,k] <- as.vector(dn[[k]]) # 4
      neighbours[,k] <- 1 + as.vector(dn[[kk + k]]) # 4
    }
    TreeList$dom <- calcTrigDiff(TreeList$dbh, as.matrix(neighbours), as.matrix(distances), kk, parameters$alpha, parameters$delta) 
    c1 <- parameters$a0 * TreeList$dom^(-parameters$a1) # c1 <- parameters$a0 + parameters$a1 / TreeList$dom
    c2 <- parameters$b0 * TreeList$dbh^(-parameters$b1) # c2 <- parameters$b0 + parameters$b1 / TreeList$dbh
    TreeList$RGR <- relativeGrowthFunction(TreeList$dbh, c1, c2, TreeList$dom)  
    TreeList$RGR <- TreeList$RGR + rnorm(n = length(TreeList$dbh), mean = 0, sd = parameters$se)
    q3 <- quantile(TreeList$RGR, probs = 0.75, na.rm = T)
    q1 <- quantile(TreeList$RGR, probs = 0.25, na.rm = T)
    outlierboundaryplus <- q3 + 1.5 * IQR(TreeList$RGR, na.rm = T)
    outlierboundaryminus <- q1 - 1.5 * IQR(TreeList$RGR, na.rm = T) # 0
    outliers[i, j] <- (length(TreeList$RGR[TreeList$RGR > outlierboundaryplus]) + length(TreeList$RGR[TreeList$RGR < outlierboundaryminus])) / length(TreeList$RGR)
    rm(q1, q3, outlierboundaryplus, outlierboundaryminus)               
    rm(dn, neighbours, distances, c1, c2)
    # Remember this year's RGR for next year's mortality calculation
    TreeList$prevRGR <- TreeList$RGR
    # Calculate and save end-of-year results 
    dm[i, j] <- mean(TreeList$dbh) # Calculate and save mean stem diameter
    sma.lm1 <- sma(RGR ~ dbh, data = TreeList, method = "SMA")
    slope1 <- coef(sma.lm1)[[2]]
    fit <- fit.ellipse(TreeList$dbh, TreeList$RGR)
    slope2 <- fit$angle
    elli <- standard.ellipse(TreeList$dbh, TreeList$RGR, confs = NULL, steps = 2)
    slope3 <- elli$theta
    TreeList$rdbh <- TreeList$dbh / mean(TreeList$dbh, na.rm = T)
    sma.lm2 <- sma(RGR ~ rdbh, data = TreeList, method = "SMA")
    Dyer1 <- coef(sma.lm2)[[2]]
    sma.lm3 <- rlm(RGR ~ rdbh, data = TreeList, method = "MM")
    Dyer2 <- coef(sma.lm3)[[2]]
    sma.lm4 <- lmrob(RGR ~ rdbh, data = TreeList)
    Dyer3 <- coef(sma.lm4)[[2]]
    sma.lm5 <- lqs(RGR ~ rdbh, data = TreeList, method = "lts")
    Dyer4 <- coef(sma.lm5)[[2]]
    # Calculate and save mean slopes
    meanSlope[i, j] <- mean(c(slope1, slope2, slope3), na.rm = T) 
    meanDyer[i, j] <- mean(c(Dyer1, Dyer2, Dyer3, Dyer4), na.rm = T) 
    # Calculate and save RGR and size correlation coefficients
    cvr[i, j] <- sd(TreeList$RGR, na.rm = T) / mean(TreeList$RGR, na.rm = T)
    cvs[i, j] <- sd(TreeList$dbh, na.rm = T) / mean(TreeList$dbh, na.rm = T)
    rm(slope1, slope2, slope3, elli, sma.lm1, sma.lm2, sma.lm3, sma.lm4, sma.lm5, 
       Dyer1, Dyer2, Dyer3, Dyer4, fit) 
    # Save diameter distribution of end year
    if((myYears[1] + i - 1) == myYears[length(myYears)])
      dend[[j]] <- TreeList$dbh 
    if(!silent) 
      cat("Simulation: ", j, " Year: ", (myYears[1] + i - 1), " Number of trees: ", length(TreeList$dbh), "\n")
    # Growth update for next year
    TreeList$dbh <- TreeList$dbh * exp(TreeList$RGR) # Modifier
  }
  rm(TreeList)
}

# Save
saveRDS(dend, paste(outputFile, "ResultsModelNonZeroAsymptote/PinusRadiata/dend", no, ".rds", sep = ""))
saveRDS(dm, paste(outputFile, "ResultsModelNonZeroAsymptote/PinusRadiata/dm", no, ".rds", sep = ""))
saveRDS(meanSlope, paste(outputFile, "ResultsModelNonZeroAsymptote/PinusRadiata/meanSlope", no, ".rds", sep = "")) 
saveRDS(meanDyer, paste(outputFile, "ResultsModelNonZeroAsymptote/PinusRadiata/meanDyer", no, ".rds", sep = "")) 

# Read
dend <- readRDS(paste(outputFile, "ResultsModelNonZeroAsymptote/PinusRadiata/dend", no, ".rds", sep = ""))
dm <- readRDS(paste(outputFile, "ResultsModelNonZeroAsymptote/PinusRadiata/dm", no, ".rds", sep = ""))
meanSlope <- readRDS(paste(outputFile, "ResultsModelNonZeroAsymptote/PinusRadiata/meanSlope", no, ".rds", sep = ""))
meanDyer <- readRDS(paste(outputFile, "ResultsModelNonZeroAsymptote/PinusRadiata/meanDyer", no, ".rds", sep = ""))


# install.packages("remotes", dep = T)
# library(remotes)
# install_github('myllym/GET')
library(GET)

# End diameter distribution
myPoints <- dmean$x # From the start of the script
mean_d <- matrix(nrow = length(myPoints), ncol = replications) 
for (i in 1 : replications) {
  obs <- density(dend[[i]], n = length(myPoints), bw = 8, kernel = "gaussian", from = 0, to = endDensity) # Check "to" with above
  mean_d[,i] <- obs$y
  rm(obs)
}
any(is.infinite(dmean$y))
any(is.na(dmean$y))
any(is.na(dmean))
any(is.infinite(mean_d))
any(is.na(mean_d))
cset <- create_curve_set(list(r = myPoints, obs = dmean$y, sim_m = mean_d))
res <- global_envelope_test(cset)
# pdf(file = paste(outputFile, no, "dbhDistribution.pdf", sep = ""))
par(mar = c(2, 4.0, 0.5, 0.5), mfrow = c(1, 1))
plot(res$r, res$obs, type = "l", lwd = 2, ylab = "", xlab = "", main = "", axes = F, xlim = c(0, 90), ylim = c(0, 0.05))
polygon(c(res$r, rev(res$r)), c(res$lo, rev(res$hi)), col = "lightgray", border = "lightgray")
lines(res$r, res$obs, lwd = 2)
# abline(h = 0, lwd = 1, lty = 2)
axis(1, cex.axis = 1.8)
axis(2, las = 1, cex.axis = 1.8)
box(lwd = 2)
# dev.off()

# Define x range and scale for all following tests and graphs
range(parameters$dm)
range(rowMeans(dm), na.rm = T)
myX <- seq(7, 30, 0.5) # Adjust manually

# Slope parameter
smean <- ksmooth(x = parameters$dm, y = parameters$meanSlope, kernel = "normal", bandwidth = 4, x.points = myX)
mean_m <- r <- matrix(nrow = length(myX), ncol = replications) 
for (i in 1 : replications) {
  dummy <- ksmooth(x = dm[,i], y = meanSlope[,i], kernel = "normal", bandwidth = 4, x.points = myX)
  r[,i] <- dummy$x
  mean_m[,i] <- dummy$y
  rm(dummy)
}
rowMeans(r)
any(is.infinite(mean_m))
any(is.na(mean_m))
cset <- create_curve_set(list(r = myX, obs = smean$y, sim_m = mean_m))
res <- global_envelope_test(cset)
options(scipen = 1)
# pdf(file = paste(outputFile, no, "Slope.pdf", sep = ""))
par(mar = c(2, 4.5, 0.5, 0.5), mfrow = c(1, 1))
plot(res$r, res$obs, type = "l", lwd = 2, ylab = "", xlab = "", main = "", axes = F, xlim = c(0, 40), ylim = c(-0.12, 0.03))
polygon(c(res$r, rev(res$r)), c(res$lo, rev(res$hi)), col = "lightgray", border = "lightgray")
lines(res$r, res$obs, lwd = 2)
abline(h = 0, lwd = 1, lty = 2)
axis(1, cex.axis = 1.8)
axis(2, las = 1, cex.axis = 1.8)
box(lwd = 2)
# dev.off()

# Dyer slopes with d / dm
smean <- ksmooth(x = parameters$dm, y = parameters$meanDyer, kernel = "normal", bandwidth = 4, x.points = myX)
mean_m <- matrix(nrow = length(myX), ncol = replications) 
for (i in 1 : replications) {
  dummy <- ksmooth(x = dm[,i], y = meanDyer[,i], kernel = "normal", bandwidth = 4, x.points = myX)
  mean_m[,i] <- dummy$y
  rm(dummy)
}
rowMeans(r)
any(is.infinite(mean_m))
any(is.na(mean_m))
cset <- create_curve_set(list(r = myX, obs = smean$y, sim_m = mean_m))
res <- global_envelope_test(cset)
options(scipen = 1)
# pdf(file = paste(outputFile, no, "Dyer.pdf", sep = ""))
par(mar = c(2, 3.8, 0.5, 0.5), mfrow = c(1, 1))
plot(res$r, res$obs, type = "l", lwd = 2, ylab = "", xlab = "", main = "", axes = F, xlim = c(0, 40), ylim = c(-0.5, 0.3))
polygon(c(res$r, rev(res$r)), c(res$lo, rev(res$hi)), col = "lightgray", border = "lightgray")
lines(res$r, res$obs, lwd = 2)
abline(h = 0, lwd = 1, lty = 2)
axis(1, cex.axis = 1.8)
axis(2, las = 1, cex.axis = 1.8)
box(lwd = 2)
# dev.off()

# Coefficient of variation of RGR
smean <- ksmooth(x = parameters$dm, y = parameters$cvr, kernel = "normal", bandwidth = 4, x.points = myX)
mean_m <- matrix(nrow = length(myX), ncol = replications) 
for (i in 1 : replications) {
  dummy <- ksmooth(x = dm[,i], y = cvr[,i], kernel = "normal", bandwidth = 4, x.points = myX)
  mean_m[,i] <- dummy$y
  rm(dummy)
}
any(is.infinite(mean_m))
any(is.na(mean_m))
cset <- create_curve_set(list(r = myX, obs = smean$y, sim_m = mean_m))
res <- global_envelope_test(cset)
# pdf(file = paste(outputFile, no, "RGRCv.pdf", sep = ""))
par(mar = c(2, 3.3, 0.5, 0.5), mfrow = c(1, 1))
plot(res$r, res$obs, type = "l", lwd = 2, ylab = "", xlab = "", main = "", axes = F, xlim = c(0, 40), ylim = c(0, 2))
polygon(c(res$r, rev(res$r)), c(res$lo, rev(res$hi)), col = "lightgray", border = "lightgray")
lines(res$r, res$obs, lwd = 2)
axis(1, cex.axis = 1.8)
axis(2, las = 1, cex.axis = 1.8)
box(lwd = 2)
# dev.off()

# Coefficient of variation of size
smean <- ksmooth(x = parameters$dm, y = parameters$cvs, kernel = "normal", bandwidth = 4, x.points = myX)
mean_m <- matrix(nrow = length(myX), ncol = replications) 
for (i in 1 : replications) {
  dummy <- ksmooth(x = dm[,i], y = cvs[,i], kernel = "normal", bandwidth = 4, x.points = myX)
  mean_m[,i] <- dummy$y
  rm(dummy)
}
any(is.infinite(mean_m))
any(is.na(mean_m))
cset <- create_curve_set(list(r = myX, obs = smean$y, sim_m = mean_m))
res <- global_envelope_test(cset)
# pdf(file = paste(outputFile, no, "SizeCv.pdf", sep = ""))
par(mar = c(2, 3.2, 0.5, 0.5), mfrow = c(1, 1))
plot(res$r, res$obs, type = "l", lwd = 2, ylab = "", xlab = "", main = "", axes = F, xlim = c(0, 40), ylim = c(0, 0.5))
polygon(c(res$r, rev(res$r)), c(res$lo, rev(res$hi)), col = "lightgray", border = "lightgray")
lines(res$r, res$obs, lwd = 2)
axis(1, cex.axis = 1.8)
axis(2, las = 1, cex.axis = 1.8)
box(lwd = 2)
# dev.off()

# Ratio of coefficients of correlation
cvratio <- matrix(nrow = dim(cvr)[1], ncol = replications) 
for (i in 1 : dim(cvr)[1]) 
  for (j in 1 : replications) 
    cvratio[i, j] <- cvr[i, j] / cvs[i, j]
smean <- ksmooth(x = parameters$dm, y = parameters$cvr / parameters$cvs, kernel = "normal", bandwidth = 4, x.points = myX)
mean_m <- matrix(nrow = length(myX), ncol = replications) 
for (i in 1 : replications) {
  dummy <- ksmooth(x = dm[,i], y = cvratio[,i], kernel = "normal", bandwidth = 4, x.points = myX)
  mean_m[,i] <- dummy$y
  rm(dummy)
}
any(is.infinite(mean_m))
any(is.na(mean_m))
cset <- create_curve_set(list(r = myX, obs = smean$y, sim_m = mean_m))
res <- global_envelope_test(cset)
# pdf(file = paste(outputFile, no, "CvRatio.pdf", sep = ""))
par(mar = c(2, 2.8, 0.5, 0.5), mfrow = c(1, 1))
plot(res$r, res$obs, type = "l", lwd = 2, ylab = "", xlab = "", main = "", axes = F, xlim = c(00, 40), ylim = c(0, 12))
polygon(c(res$r, rev(res$r)), c(res$lo, rev(res$hi)), col = "lightgray", border = "lightgray")
lines(res$r, res$obs, lwd = 2)
axis(1, cex.axis = 1.8)
axis(2, las = 1, cex.axis = 1.8)
box(lwd = 2)
# dev.off()

# Checking outlier devlopment
rowMeans(outliers)
mean(rowMeans(outliers))
