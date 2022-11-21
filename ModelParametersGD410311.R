# Model parameters and results GD 410311 data, kk = 30
#----------------------
meanSlope <- c(-0.0563344083131597, -0.0150156162222864, 0.0084634313659726, 0.0089951814592860, 0.0041660549741623, 0.0033906605488606,
               0.0045646528546510, 0.0034539626027548, 0.0026180952022818, 0.0032113097646365, 0.0036023320009480, 0.0018821886414112,
               0.0022092206400770, 0.0016787635365009, 0.0011892076066931, 0.0012089918762369) # Mean observed slope for each survey year
meanDyer <- c(-0.367558619551449, -0.099481819225503, 0.109508318834669, 0.140709745380605, 0.072752155557361, 0.066474504686719,
              0.100757521879589, 0.086980568439576, 0.091286664994956, 0.068316947650614, 0.088720261743394, 0.059908804451310,
              0.072602134359736, 0.049855800156038, 0.040946023188257, 0.032279900543308) # Mean observed slope of RGR - d/dm for each survey year
cvr <- c(0.15990793178492, 0.16914956034061, 0.26031271068214, 0.27341388442276, 0.34543116676759, 0.39474269102557, 0.49831976796544,
         0.69605472982769, 0.54571896152984, 0.77263145670552, 0.84606691729913, 1.15216698731959, 1.14546561979825, 0.81135629410899,
         0.84285636191031, 1.02495508980313) # RGR coefficient of variation
cvs <- c(0.15563177622112, 0.12785917087270, 0.12686459326890, 0.13602843708770, 0.15061790523917, 0.15812474112367, 0.16516046125379, 
         0.17726366726070, 0.18211572198595, 0.19359246361793, 0.20225394742283, 0.20654581854713, 0.21453844208032, 0.23047160716085,
         0.21899098989998, 0.23369315572114) # Size coefficient of variation
dm <- c(6.8845679012346, 11.5012345679012, 14.9000000000000, 17.5728395061728, 19.8358024691358, 21.2753086419753, 22.3648148148148, 
        23.4648148148148, 24.2981366459627, 25.5118750000000, 26.2798742138365, 27.1647435897436, 27.8535947712418, 28.3907894736842,
        29.0417218543046, 30.5943661971831) # Mean arithmetic diameter for each survey year
corrGroDomDyer <- 0.5258487471882 # Correlation between growth dominance and the Dyer slope
alpha <- 1 # A/symmetry exponent in hyperbolic tangent index; default
delta <- 2 # Distance exponent; default
c1 <- 0.13013786382271 # c1 parameter of median population trend
c2 <- 238.53880973128969 # c2 parameter of median population trend
a0 <- 0.050097119486697  # First parameter for estimating c1
a1 <- 1.343780998449513  # Second parameter for estimating c1
b0 <- 0.212724637573473  # First parameter for estimating c2
b1 <- -1.333893626485041 # Second parameter for estimating c2
rangeC1 <- c(0.077847414335601, 0.837113103730480) # Minimum/maximum value of c1 estimated from a0 and a1
rangeC2 <- c(1.3068781728896, 34.5266174796971) # Minimum/maximum value of c2 estimated from b0 and b1
bias <- 0.0010450603943133 # Bias
rbias <- 0.011317854993976 # Relative bias
rmse <- 0.032486778788501 # RMSE
rrmse <- 0.35182718008486 # Relative RMSE
eff <- 0.94051368355886 # Efficiency
se <- 0.03250624060173 # Standard error of residuals of RGR estimation
m0 <- -6.079518980422371 # Mortality model, intercept
m1 <- 48.871136477431108 # Mortality model, 1 / dbh
m2 <- -30.091105209488887 # Mortality model, prevRGR
parameters <- list(meanSlope = meanSlope, meanDyer = meanDyer, cvr = cvr, cvs = cvs, dm = dm, corrGroDomDyer = corrGroDomDyer,
                   alpha = alpha, delta = delta, c1 = c1, c2 = c2, a0 = a0, a1 = a1, b0 = b0, b1 = b1, rangeC1 = rangeC1, 
                   rangeC2 = rangeC2, bias = bias, rbias = rbias, rmse = rmse, rrmse = rrmse, eff = eff, 
                   se = se, m0 = m0, m1 = m1, m2 = m2)

