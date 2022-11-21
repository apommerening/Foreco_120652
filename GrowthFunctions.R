calcBasalArea <- function(dbh) 
  return(pi * (dbh / 200)^2)


calcTrigDiff <- function(mark, neighbours, distance, k, alpha, dexp) {
  td <- c()
  for (i in 1 : length(mark)) {
    tsum <- 0
    distsum <- 0
    for (j in 1 : k) {
      w <- calcBasalArea(mark[neighbours[i, j]]) / distance[i, j]^dexp
      tsum <- tsum + mark[i]^(2 * alpha) / (mark[i]^(2 * alpha) + mark[neighbours[i, j]]^(2 * alpha)) * w #  exp(-w) # w
      distsum <- distsum + w
    }
    td[i] <- tsum / distsum 
  }
  return(td)
}

relativeGrowthFunction <- function(y, xk, xp, xcomp) {
  return(exp(-xk * y * (1 - exp(-xp * xcomp)))) 
}

estimateStateSpace <- function(dbh, RGR, dend, period, xk, xp, xcomp) {	# , b			
  pd <- c()
  for(i in 1 : length(dbh)) { # length of observations
    M <- 1
    yest <- dbh[i]
    for(j in 1 : period[i]) { # period
      p <- relativeGrowthFunction(yest, xk[i], xp[i], xcomp[i]) 
      yest <- yest * exp(p)
      M <- M * exp(p)
    }
    pd[i] <- log(M^(1 / period[i]))
  }
  return(pd)
}

estimateC1 <- function(myPar, indepVar1, indepVar2) {
  return(myPar[1] * indepVar1^{-myPar[2]})
}

estimateC2 <- function(myPar, indepVar1, indepVar2) {
  return(myPar[1] * indepVar1^{-myPar[2]})
}

lossD <- function(abdn, dbh, dend, period, RGR, xcomp, xdr, varr, fact, denom) {
  c1 <- estimateC1(abdn[1 : 2], xcomp, xdr) # abdn[1] + abdn[2] * xcomp^{-abdn[3]} # abdn[1] + abdn[2] / xcomp # abdn[1] + abdn[2] * xcomp # abdn[1] * xcomp^{-abdn[2]} # xcomp
  c2 <- estimateC2(abdn[3 : 4], dbh, xdr) # abdn[4] + abdn[5] * dbh^{-abdn[6]} # abdn[3] + abdn[4] / dbh # abdn[3] * dbh^{-abdn[4]} # abdn[3] + abdn[4] * dbh # abdn[3] * dbh^{-abdn[4]} # dbh
  di <- estimateStateSpace(dbh, RGR, dend, period, c1, c2, xcomp) # xcomp, xdr
  dev <- 1 / varr * (RGR - di)^2
  ret <- sum(dev, na.rm = T) / sum(1 / varr, na.rm = T)
  if(denom == 0)
    ret <- ret + fact * any(abdn[1 : 4] < 0) else ret <- ret + fact * any(abdn[1 : 4] < 0) + sum(abs(abdn[1 : 4]) / denom)
  return(ret)
}  

loss.ML <- function(abdn, dbh, dend, period, RGR, xcomp, xdr, fact, denom) {
  c1 <- estimateC1(abdn[1 : 2], xcomp, xdr) # abdn[1] + abdn[2] * xcomp^{-abdn[3]} # abdn[1] + abdn[2] / xcomp # abdn[1] + abdn[2] * xcomp # abdn[1] * xcomp^{-abdn[2]} # xcomp
  c2 <- estimateC2(abdn[3 : 4], dbh, xdr) # abdn[4] + abdn[5] * dbh^{-abdn[6]} # abdn[3] + abdn[4] / dbh # abdn[3] * dbh^{-abdn[4]} # abdn[3] + abdn[4] * dbh # abdn[3] * dbh^{-abdn[4]} # dbh
  di <- estimateStateSpace(dbh, RGR, dend, period, c1, c2, xcomp) # xcomp, xdr 
  dev <- dnorm(RGR, mean = di, sd = exp(abdn[length(abdn)]), log = TRUE)
  if(denom == 0)
    ret <- sum(dev, na.rm = TRUE) - fact * any(abdn[1 : 4] < 0) else ret <- sum(dev, na.rm = TRUE) - 2 * sum(abdn[1 : 4] / denom) - fact * any(abdn[1 : 4] < 0)
  return(ret)
}

estimateGrowth <- function(y, abdn, xcomp, xdr) {
  c1 <- estimateC1(abdn[1 : 2], xcomp, xdr) # abdn[1] + abdn[2] * xcomp^{-abdn[3]} # abdn[1] + abdn[2] / xcomp # abdn[1] + abdn[2] * xcomp # abdn[1] * xcomp^{-abdn[2]} # xcomp
  c2 <- estimateC2(abdn[3 : 4], y, xdr) # abdn[4] + abdn[5] * y^{-abdn[6]} # abdn[3] + abdn[4] / y # abdn[3] * y^{-abdn[4]} # abdn[3] + abdn[4] * y # abdn[3] * y^{-abdn[4]} # y
  return(relativeGrowthFunction(y, c1, c2, xcomp)) # xcomp, xdr
}
