rm(list = ls())
options(digits = 6, width = 50)

# File path
sourcePath <- dataPath <- "C:/Users/arng0001/Dropbox/Gwaith/Cyhoeddi/RGRmodel/Data/"
destPath <- "C:/Users/arng0001/Dropbox/Gwaith/Cyhoeddi/RGRmodel/Code/"
outputFile <- "C:/Users/arng0001/Dropbox/Gwaith/Cyhoeddi/RGRmodel/Results/"

# install.packages("spatstat", dep = T)
library(spatstat)
library(Rcpp) 

# Sys.setenv("LANGUAGE"="EN")

sourceCpp(paste(destPath, "findNeighbours.cpp", sep = ""))
source(paste(destPath, "GrowthFunctions.R", sep = ""))

fileName <- c("RadiataPineData.txt")
incRegres <- NA
myData  <- read.table(paste(sourcePath, fileName, sep = ""), header = T)
names(myData)
myData$plot_id <- sprintf("%02d", as.numeric(myData$plot_id))
myData$id <- paste(myData$trial_id, myData$block_id, myData$plot_id, sep = "")
table(myData$id)
PlotNames <- as.numeric(unique(myData$id))
myData$plotno <- as.numeric(myData$id)
table(myData$plotno)
TreeList  <- myData
rm(myData)
# i <- 1
# j <- 1
myAngle <- c(0, 0, 0, 0.02, 0.08, 0.13, 0.12, 0.12, 0.12, 0, 0, 0.02)
firstData <- T
# Calculate mean annual growth rates
for (i in 1 : length(PlotNames)) {
  myData  <- TreeList[TreeList$id == PlotNames[i], ]
  xmax <- ceiling(max(myData$X))
  ymax <- ceiling(max(myData$Y))
  xmin <- floor(min(myData$X))
  ymin <- floor(min(myData$Y))
  testp <- ppp(x = myData$X, y = myData$Y, window = owin(c(xmin, xmax), c(ymin, ymax)))
  if(myAngle[i] > 0) {
    testp <- rotate(testp, angle = myAngle[i]) # 0.75
    myData$X <- testp$x
    myData$Y <- testp$y
  }
  rm(testp)
  myData$X <- myData$X - min(myData$X) + 0.1
  myData$Y <- myData$Y - min(myData$Y) + 0.1
  rm(xmin, ymin)
  xmax <- ceiling(max(myData$X))
  ymax <- ceiling(max(myData$Y))
  testp <- ppp(x = myData$X, y = myData$Y, window = owin(c(0, xmax), c(0, ymax)))
  surveyYears <- as.numeric(names(table(substr(myData$Date, 1, 4))))
  myData$year <- as.numeric(substr(myData$Date, 1, 4))
  myData$x <- myData$X
  myData$y <- myData$Y
  myData$Treeno <- myData$tree_id
  myData$Species <- 1
  myData$xmax <- xmax
  myData$ymax <- ymax
  myData <- myData[c("plotno", "Treeno", "Species", "dbh", "x", "y", "year", "xmax", "ymax")]
  firstSurvey <- T
  years <- table(myData$year)
  myYears <- as.numeric(names(years))
  lostTrees <- c()
  uniqueIDs <- unique(myData$Treeno)
  rm(years, testp)
  for (k in 1 : length(uniqueIDs)) { # Check diameter trajectories
    if(k == 1) {
      par(mar = c(2.0, 2.7, 0.5, 0.8))
      plot(myYears, myData$dbh[myData$Treeno == uniqueIDs[k]], type = "l", axes = FALSE, ylab = "", xlab = "", main = "", 
         cex = 1.0, lwd = 1, col = "black", ylim = c(0, 60), xlim = c(1990, 2015))
      axis(side = 1, lwd = 2, las = 1, cex.axis = 1.7)
      axis(side = 2, lwd = 2, las = 1, cex.axis = 1.7)
      box(lwd = 2)
    }
    else {
      if(length(myYears) == length(myData$dbh[myData$Treeno == uniqueIDs[k]]))
        lines(myYears, myData$dbh[myData$Treeno == uniqueIDs[k]], lwd = 1, col = "black")
      else {
        lostTrees <- c(lostTrees, uniqueIDs[k])
        lines(myYears[1 : length(myData$dbh[myData$Treeno == uniqueIDs[k]])], myData$dbh[myData$Treeno == uniqueIDs[k]], lwd = 1, col = "red")
      }
    }
  }
  # dev.off()
  RGR <- matrix(data = NA, nrow = length(myYears) - 1, ncol = length(uniqueIDs))
  for (k in 1 : length(uniqueIDs)) { # Check RGR trajectories
    dummy <- diff(log(myData$dbh[myData$Treeno == uniqueIDs[k]]), lag = 1, differences = 1) / diff(myData$year[myData$Treeno == uniqueIDs[k]], lag = 1, differences = 1)
    RGR[, k] <- c(dummy, rep(NA, length(myYears) - 1 - length(dummy)))
    rm(dummy)
  }
  mRGR <- apply(RGR, 1, median, na.rm=T)
  rm(RGR)
  for (k in 1 : length(uniqueIDs)) { # Check RGR trajectories
    if(k == 1) {
      par(mar = c(2.0, 3.2, 0.5, 0.8))
      plot(myYears[-length(myYears)], diff(log(myData$dbh[myData$Treeno == uniqueIDs[k]]), lag = 1, differences = 1) / diff(myData$year[myData$Treeno == uniqueIDs[k]], lag = 1, differences = 1), type = "l", axes = FALSE, ylab = "", xlab = "", main = "", 
           cex = 1.0, lwd = 1, col = "black", ylim = c(-0.15, 2.0), xlim = c(1990, 2015))
      abline(h = 0, lwd = 1, lty = 2)
      axis(side = 1, lwd = 2, las = 1, cex.axis = 1.7)
      axis(side = 2, lwd = 2, las = 1, cex.axis = 1.7)
      box(lwd = 2)
    }
    else {
      if(length(myYears) == length(myData$dbh[myData$Treeno == uniqueIDs[k]]))
        lines(myYears[-length(myYears)], diff(log(myData$dbh[myData$Treeno == uniqueIDs[k]]), lag = 1, differences = 1) / diff(myData$year[myData$Treeno == uniqueIDs[k]], lag = 1, differences = 1), lwd = 1, col = "black")
      else if(length(myData$dbh[myData$Treeno == uniqueIDs[k]]) > 1) {
        lostTrees <- c(lostTrees, uniqueIDs[k])
        lines(myYears[1 : (length(myData$dbh[myData$Treeno == uniqueIDs[k]]) - 1)], diff(log(myData$dbh[myData$Treeno == uniqueIDs[k]]), lag = 1, differences = 1) / diff(myData$year[myData$Treeno == uniqueIDs[k]], lag = 1, differences = 1), lwd = 1, col = "red")
      }
    }
  }
  for (k in 1 : length(uniqueIDs)) { # Check RGR trajectories
    if(k == 1) {
      par(mar = c(2.0, 3.5, 0.5, 0.8))
      plot(myYears[-length(myYears)], (diff(log(myData$dbh[myData$Treeno == uniqueIDs[k]]), lag = 1, differences = 1) / diff(myData$year[myData$Treeno == uniqueIDs[k]], lag = 1, differences = 1)) - mRGR, type = "l", axes = FALSE, ylab = "", xlab = "", main = "", 
           cex = 1.0, lwd = 1, col = "black", ylim = c(-0.5, 1.3), xlim = c(1990, 2015))
      abline(h = 0, lwd = 1, lty = 2)
      axis(side = 1, lwd = 2, las = 1, cex.axis = 1.7)
      axis(side = 2, lwd = 2, las = 1, cex.axis = 1.7)
      box(lwd = 2)
    }
    else {
      if(length(myYears) == length(myData$dbh[myData$Treeno == uniqueIDs[k]]))
        lines(myYears[-length(myYears)], (diff(log(myData$dbh[myData$Treeno == uniqueIDs[k]]), lag = 1, differences = 1) / diff(myData$year[myData$Treeno == uniqueIDs[k]], lag = 1, differences = 1)) - mRGR, lwd = 1, col = "black")
      else if(length(myData$dbh[myData$Treeno == uniqueIDs[k]]) > 1) {
        lostTrees <- c(lostTrees, uniqueIDs[k])
        lines(myYears[1 : (length(myData$dbh[myData$Treeno == uniqueIDs[k]]) - 1)], (diff(log(myData$dbh[myData$Treeno == uniqueIDs[k]]), lag = 1, differences = 1) / diff(myData$year[myData$Treeno == uniqueIDs[k]], lag = 1, differences = 1)) - mRGR[1 : (length(myData$dbh[myData$Treeno == uniqueIDs[k]]) - 1)], lwd = 1, col = "red")
      }
    }
  }
  rm(uniqueIDs, mRGR)
  for (j in 1 : (length(surveyYears) - 1)) {
    dataOneYear <- myData[myData$year == surveyYears[j], ]
    dataOneYear <- dataOneYear[order(dataOneYear$Treeno, decreasing = FALSE), ]
    if(firstSurvey == F)
        dataOneYear <- merge(dataOneYear, dataOtherYear.s, all.x = TRUE, sort = TRUE)
    else {
      dataOneYear$prevRGR <- NA
    }
      firstSurvey <- F
      dataOtherYear <- myData[myData$year == surveyYears[j + 1], ]
      dataOtherYear <- dataOtherYear[order(dataOtherYear$Treeno, decreasing = FALSE), ]
      dataOtherYear$dbhOther <- dataOtherYear$dbh
      dataOtherYear.s <- dataOtherYear[c("Treeno", "dbhOther")]
      dataOneYear <- merge(dataOneYear, dataOtherYear.s, all.x = TRUE, sort = TRUE)  
      dataOneYear$dead <- 1 # glm() seems to go for mort = 1
      dataOneYear$dead[!is.na(dataOneYear$dbhOther)] <- 0
      dataOneYear$dm <- median(dataOneYear$dbh, na.rm = T) # mean(dataOneYear$dbh, na.rm = T)
      dataOneYear$rank <- rank(dataOneYear$dbh, ties.method = "average") / length(dataOneYear$dbh)
      dataOneYear$period <- surveyYears[j + 1] - surveyYears[j]
      dataOneYear$AGR <- (dataOneYear$dbhOther - dataOneYear$dbh) / (surveyYears[j + 1] - surveyYears[j])
      dataOneYear$RGR <- (log(dataOneYear$dbhOther) - log(dataOneYear$dbh)) / (surveyYears[j + 1] - surveyYears[j])
      dataOneYear$d.end <- dataOneYear$dbhOther
      dataOneYear$AGR.g <- (calcBasalArea(dataOneYear$dbhOther) - calcBasalArea(dataOneYear$dbh)) / (surveyYears[j + 1] - surveyYears[j])
      dataOneYear$RGR.g <- log(calcBasalArea(dataOneYear$dbhOther) / calcBasalArea(dataOneYear$dbh)) / (surveyYears[j + 1] - surveyYears[j])
      dataOtherYear.s <- dataOneYear[c("Treeno", "RGR")]
      dataOtherYear.s$prevRGR <- dataOtherYear.s$RGR
      dataOtherYear.s <- dataOtherYear.s[c("Treeno", "prevRGR")]
      if(firstData == TRUE)
        incRegres <- dataOneYear[c("plotno", "Treeno", "Species", "x", "y", "dbh", "d.end", "dm", "rank", "AGR", "RGR", "prevRGR", "AGR.g", "RGR.g", "dead", "year", "period", "xmax", "ymax")] 
      else {
        dummy <- dataOneYear[c("plotno", "Treeno", "Species", "x", "y", "dbh", "d.end", "dm", "rank", "AGR", "RGR", "prevRGR", "AGR.g", "RGR.g", "dead", "year", "period", "xmax", "ymax")] 
        incRegres <- rbind(incRegres, dummy)
      }
      firstData <- FALSE
      rm(dataOneYear, dataOtherYear)
    }
}
rm(myData)

head(incRegres)
tail(incRegres)
table(incRegres$plotno)

write.table(incRegres, file = paste(dataPath, "giRegresGD.txt", sep = ""))

no <- 410311 
incRegres  <- incRegres[incRegres$plotno == no, ] 
source(paste(destPath, "ModelParametersGD", no, ".R", sep = "")) 
alpha <- parameters$alpha
exponent <- parameters$delta
# Dominance
kk <- 30
incRegres$ui3 <- incRegres$dr <- matrix(NA, nrow = length(incRegres$Treeno), ncol = length(table(incRegres$year)))
# i <- 1
for (i in 1 : length(table(incRegres$year))) {
  incRegres1 <- incRegres[incRegres$year == as.numeric(names(table(incRegres$year)[i])),]
  dr <- incRegres1$dbh / mean(incRegres1$dbh, na.rm = T) # median(incRegres1$dbh, na.rm = T)
  dn <- findNeighbours(incRegres1$xmax[1], incRegres1$ymax[1], incRegres1$x, incRegres1$y, kk, 1) # 1 means periodic boundary conditions
  neighbours <- distances <- matrix(NA, nrow = length(incRegres1$x), ncol = kk)
  for (j in 1 : kk) {
    distances[,j] <- as.vector(dn[[j]]) # 4
    neighbours[,j] <- 1 + as.vector(dn[[kk + j]]) # 4
  }
  ui3 <- calcTrigDiff(incRegres1$dbh, as.matrix(neighbours), as.matrix(distances), kk, alpha, exponent) # alpha[z] # alpha = 0.5 so that exponent = 1
  df2 <- data.frame(Treeno = incRegres1$Treeno, year = incRegres1$year, ui3)
  test <- merge(incRegres, df2, by = c("Treeno", "year"), all.x = T, all.y = F)
  incRegres$ui3[,i] <- test$ui3.y
  rm(test, df2)
  df2 <- data.frame(Treeno = incRegres1$Treeno, year = incRegres1$year, dr)
  test <- merge(incRegres, df2, by = c("Treeno", "year"), all.x = T, all.y = F)
  incRegres$dr[,i] <- test$dr.y
  rm(incRegres1, neighbours, distances, ui3, df2, dr, dn, test)
}
ui3 <- rowSums(incRegres$ui3, na.rm = T)  
incRegres$ui3
ui3[1 : 10]
incRegres$ui3 <- ui3
rm(ui3)
range(incRegres$ui3)
dr <- rowSums(incRegres$dr, na.rm = T)  
incRegres$dr <- dr
rm(dr)
range(incRegres$dr)


incRegres <- incRegres[incRegres$dbh > 0, ] 
incRegres$RGR[incRegres$RGR < 0] <- 0

head(incRegres)
incRegres[100 : 200, c(7, 9, 12)]
any(is.na(incRegres$d.end))
which(is.na(incRegres$d.end))
table(incRegres$year)
tapply(incRegres$Treeno, incRegres$year, length)
length(incRegres$dead[incRegres$dead == 1]) / length(incRegres$dead)
length(incRegres$dead[incRegres$dead == 1])
hist(incRegres$prevRGR)



xres <- glm(dead ~ dbh + ui3, family = binomial, data = incRegres, offset = log(period))
(r2 <- (xres$null.deviance - deviance(xres)) / xres$null.deviance)
dp <- sum(residuals(xres, type = "pearson")^2) / xres$df.res
summary(xres)

xres <- glm(dead ~ dbh + dr, family = binomial, data = incRegres, offset = log(period))
(r2 <- (xres$null.deviance - deviance(xres)) / xres$null.deviance)
dp <- sum(residuals(xres, type = "pearson")^2) / xres$df.res
summary(xres)

# This is the one
incRegres$rdbh <- 1 / incRegres$dbh
options(digits = 4)
xres <- glm(dead ~ rdbh + prevRGR, family = binomial, data = incRegres, offset = log(period))
(r2 <- (xres$null.deviance - deviance(xres)) / xres$null.deviance)
dp <- sum(residuals(xres, type = "pearson")^2) / xres$df.res
summary(xres)
options(digits = 22)
summary(xres)$coefficients[1 : 3]

incRegres$rdbh <- 1 / incRegres$dbh
gompit.m <- glm(dead ~ rdbh + prevRGR,
                  offset = log(period), 
                  data = incRegres, family = binomial(link = "cloglog"))
summary(gompit.m)
(r2 <- (gompit.m$null.deviance - deviance(gompit.m)) / gompit.m$null.deviance)
dp <- sum(residuals(gompit.m, type = "pearson")^2) / gompit.m$df.res
options(digits = 22)
summary(gompit.m)$coefficients[1 : 3]

# Check
obs <- hist(incRegres$dbh[incRegres$dead == 1], breaks = seq(0, 40, by = 5), include.lowest = T, right = F, plot = F)
obs$counts <- obs$counts / length(incRegres$dbh) # length(incRegres$dbh[incRegres$dead == 1])
par(mar = c(2, 4.8, 0.5, 0.5))
plot(obs, las = 1, cex.axis = 1.7, ylim = c(0, 0.02), main = "", xlab = "", ylab = "")
box(lwd = 2)

# Set up cut-off values 
breaks <- seq(0, 40, by = 5) # seq(5, 65, by = 5)
# Specify interval/bin labels
tags <- seq(1, length(breaks) - 1, 1) # c("[0-2)","[2-4)", "[4-6)", "[6-8)", "[8-10)", "[10-12)","[12-14)", "[14-16)","[16-18)", "[18-20)")
length(breaks)
length(tags)
# Bucketing values into bins
group_tags <- cut(incRegres$dbh, breaks = breaks, include.lowest = T, right = F, labels = tags)
any(is.na(group_tags))

# Checking whether tags are correct
deadOnes <- incRegres$dbh[incRegres$dead == 1]
group_tags <- cut(deadOnes, breaks = breaks, include.lowest = T, right = F, labels = tags)
for (i in 1 : length(tags))
  print(length(deadOnes[group_tags == tags[i]]) / length(incRegres$dbh))
obs$counts

group_tags <- cut(incRegres$dbh, breaks = breaks, include.lowest = T, right = F, labels = tags)
dbh <- tapply(incRegres$dbh, group_tags, mean, na.rm = T)
prevRGR <- tapply(incRegres$prevRGR, group_tags, mean, na.rm = T)
# Logit
ei <- 1 - (1 / (1 + exp(summary(xres)$coefficients[1] + summary(xres)$coefficients[2] * 1 / dbh + summary(xres)$coefficients[3] * prevRGR))) # + summary(xres)$coefficients[4] * dbh^2)))
# Cloglog
x <- summary(gompit.m)$coefficients[1] + summary(gompit.m)$coefficients[2] * 1 / dbh + summary(gompit.m)$coefficients[3] * prevRGR
dt <- 1
ei <- 1 - exp(-exp(x) * dt)
ei[is.na(ei)] <- 0
length(obs$counts)

par(mar = c(2, 4.8, 0.5, 0.5))
barplot(rbind(obs$counts, ei[1 : length(obs$counts)]), names = round(dbh[1 : length(obs$counts)], digits = 0), 
        beside = TRUE, ylab = "", cex.names = 1.7, col = c("white", "black"), main = "", las = 1, cex.axis = 1.7,
        ylim = c(0, 0.01))
box(lwd = 2)
# Genetic distance (Gregorius, 1974; Westphal et al., 2006)
0.5 * sum(abs(obs$counts - ei[1 : length(obs$counts)])) # With dbh^2: 0.037358658283735316, without dbh^2: 0.033904079607174981 










xres <- glm(dead ~ rdbh + prevRGR + dbh2, family = poisson, data = incRegres, offset = log(period))
(r2 <- (xres$null.deviance - deviance(xres)) / xres$null.deviance)
dp <- sum(residuals(xres, type = "pearson")^2) / xres$df.res
summary(xres)
options(digits = 22)
summary(xres)$coefficients[1 : 4]

dbh <- c(10, 20, 30, 40, 50)
1 / dbh

xres <- glm(dead ~ rdbh + dr, family = binomial, data = incRegres, offset = log(period))
(r2 <- (xres$null.deviance - deviance(xres)) / xres$null.deviance)
dp <- sum(residuals(xres, type = "pearson")^2) / xres$df.res
summary(xres)

par(mar = c(2, 5, 0.5, 0.5))
plot(log(incRegres$dbh[incRegres$RGR > 0]),incRegres$RGR[incRegres$RGR > 0], pch = 16, xlab = "", ylab = "", axes = FALSE)	# xlim = c(0, 100), ylim = c(0, 35),
axis(1, lwd = 2, cex.axis = 1.8)
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2)

lmr <- lm(RGR ~ log(incRegres$dbh), data = incRegres)
summary(lmr)
curve(lmr$coefficients[1] + x * lmr$coefficients[2], from = 2,  to = 4, lwd = 2, lty = 1, col = "red", add = TRUE) 
