################################################################################
# Reproduction of Gavin Simpsons work on additive models from
# http://www.fromthebottomoftheheap.net/ 
# With a view to creating a process for analysis of Landsat veg time series
#
# Bart Huntley

rm(list=ls())

## this section "Additive modelling and the HadCRUT3v global mean temperature 
## series" from:
## http://www.fromthebottomoftheheap.net/2011/06/12/additive-modelling-and-the-hadcrut3v-global-mean-temperature-series/

## This technique fits a local, and not global, model to the entire time series 
# and see how that informs our knowledge of trends in the recent period. NOTE
## this method is ONLY MODELLING ANNUAL DATA

## load the packages and code we need
require(mgcv)
require(nlme)
## load custom functions
tmp <- tempfile()
download.file("https://github.com/gavinsimpson/random_code/raw/master/derivFun.R",
              tmp)
source(tmp)
tmp <- tempfile()
download.file("https://github.com/gavinsimpson/random_code/raw/master/tsDiagGamm.R",
              tmp)
source(tmp)

## Global temperatures
URL <- url("http://www.cru.uea.ac.uk/cru/data/temperature/HadCRUT3v-gl.dat")
gtemp <- read.table(URL, fill = TRUE)
## Don't need the even rows
gtemp <- gtemp[-seq(2, nrow(gtemp), by = 2), ]
## set the Year as rownames
rownames(gtemp) <- gtemp[,1]
## Add colnames
colnames(gtemp) <- c("Year", month.abb, "Annual")
## Data for 2011 incomplete so work only with 1850-2010 data series
gtemp <- gtemp[-nrow(gtemp), ]
## Plot the data
ylab <- expression(Temperature~Anomaly~(1961-1990)~degree*C)
plot(Annual ~ Year, data = gtemp, type = "o", ylab = ylab)

## Looking at the plot, we can see that the level of the global annual mean 
## temperature record has varied substantially over the 160 years of observations. 
## To fit a global, linear trend to the entire data would make little sense — 
## clearly such a model would not provide a good fit to the data, failing to 
## describe the relationship in temperature over time.

# The additive model (without any correlation structure at this stage) is 
# fitted and summarised as follows
m1 <- gamm(Annual ~ s(Year, k = 20), data = gtemp)
summary(m1$gam)

# This smoother explains 89% of the variance in the data, uses 12 .52
# degrees of freedom and is statistically significant, in the sense that the 
# fitted smoother is different from a null model. We should not take this 
# p-value at face value however, as these data are a time series and the 
# standard errors on the fitted smoother are likely to be overly narrow. The ACF 
# and partial ACF can be used to determine what types of time series model might 
# be required for the residuals.

## look at autocorrelation in residuals:
acf(resid(m1$lme, type = "normalized"))
## ...wait... look at the plot, only then do...
pacf(resid(m1$lme, type = "normalized"))
## seems like like sharp cut-off in ACF and PACF - AR terms probably best

## ...so fit the AR1
m2 <- gamm(Annual ~ s(Year, k = 30), data = gtemp,
           correlation = corARMA(form = ~ Year, p = 1))
## ...and fit the AR2
m3 <- gamm(Annual ~ s(Year, k = 30), data = gtemp,
           correlation = corARMA(form = ~ Year, p = 2))

# NOTE TO SELF - these gam models contain a list of 2 items, 1 an object of class
# lme and the other of class gam. Therefore We can use the anova() method for 
# "lme" objects to assess whether the models with the correlation structures fit 
# the data better than the original model.
anova(m1$lme, m2$lme, m3$lme)
#Model df       AIC       BIC   logLik   Test   L.Ratio p-value
#m1$lme     1  4 -280.0600 -267.6605 144.0300                         
#m2$lme     2  5 -305.7886 -290.2893 157.8943 1 vs 2 27.728686  <.0001
#m3$lme     3  6 -304.7318 -286.1326 158.3659 2 vs 3  0.943197  0.3315
# AR(2) is best - m2

# plot
plot(m2$gam, residuals = TRUE, pch = 19, cex = 0.75)


# Some diagnostic plots can be be produced using tsDiagGamm() function 
with(gtemp, tsDiagGamm(m2, timevar = Year, observed = Annual))

plot(Annual ~ Year, data = gtemp, type = "p", ylab = ylab)
pdat <- with(gtemp,
             data.frame(Year = seq(min(Year), max(Year),
                                   length = 200)))
p1 <- predict(m1$gam, newdata = pdat)
p2 <- predict(m2$gam, newdata = pdat)
lines(p1 ~ Year, data = pdat, col = "red")
lines(p2 ~ Year, data = pdat, col = "blue")
legend("topleft",
       legend = c("Uncorrelated Errors","AR(1) Errors"),
       bty = "n", col = c("red","blue"), lty = 1)
# This is important as the Deriv function from source doesn't work
Deriv <- function(mod, n = 200, eps = 1e-7, newdata) {
        #if(isTRUE(all.equal(class(mod), "list")))
        mod <- mod$gam
        m.terms <- attr(terms(mod), "term.labels")
        if(missing(newdata)) {
                newD <- sapply(model.frame(mod)[, m.terms, drop = FALSE],
                               function(x) seq(min(x), max(x), length = n))
                names(newD) <- m.terms
        } else {
                newD <- newdata
        }
        X0 <- predict(mod, data.frame(newD), type = "lpmatrix")
        newD <- newD + eps
        X1 <- predict(mod, data.frame(newD), type = "lpmatrix")
        Xp <- (X1 - X0) / eps
        Xp.r <- NROW(Xp)
        Xp.c <- NCOL(Xp)
        ## dims of bs
        bs.dims <- sapply(mod$smooth, "[[", "bs.dim") - 1
        # number of smooth terms
        t.labs <- attr(mod$terms, "term.labels")
        nt <- length(t.labs)
        ## list to hold the derivatives
        lD <- vector(mode = "list", length = nt)
        names(lD) <- t.labs
        for(i in seq_len(nt)) {
                Xi <- Xp * 0
                want <- grep(t.labs[i], colnames(X1))
                Xi[, want] <- Xp[, want]
                df <- Xi %*% coef(mod)
                df.sd <- rowSums(Xi %*% mod$Vp * Xi)^.5
                lD[[i]] <- list(deriv = df, se.deriv = df.sd)
                ## Xi <- Xp * 0 ##matrix(0, nrow = Xp.r, ncol = Xp.c)
                ## J <- bs.dims[i]
                ## Xi[,(i-1) * J + 1:J + 1] <- Xp[,(i-1) * J + 1:J +1]
                ## df <- Xi %*% coef(mod)
                ## df.sd <- rowSums(Xi %*% mod$Vp * Xi)^.5
                ## lD[[i]] <- list(deriv = df, se.deriv = df.sd)
        }
        class(lD) <- "Deriv"
        lD$gamModel <- mod
        lD$eps <- eps
        lD$eval <- newD - eps
        return(lD)
}

# plot that draws a time series of first derivatives with a confidence interval
# Periods where zero is not included in confidence interval can be coloured to 
# show important periods of change (red for decreasing, and blue for increasing).
# sizer is colouring on/off and alpha is coverage for confidence intervals
m2.d <- Deriv(m2, n = 200)
plot(m2.d, sizer = TRUE, alpha = 0.01)

# this plot more meaningful as it is over the original series
plot(Annual ~ Year, data = gtemp, type = "p", ylab = ylab)
lines(p2 ~ Year, data = pdat)
CI <- confint(m2.d, alpha = 0.01)
S <- signifD(p2, m2.d$Year$deriv, CI$Year$upper, CI$Year$lower,
             eval = 0)
lines(S$incr ~ Year, data = pdat, lwd = 3, col = "blue")
lines(S$decr ~ Year, data = pdat, lwd = 3, col = "red")
# The derivatives suggest two periods of significant increase in temperature 
# (at the 99% level)

# Another observation worth making is that the fitted spline is based on the ML 
# estimates of the coefficients that describe the spline. Each of these 
# coefficients is subject to uncertainty. The set of coefficients and their 
# standard errors form a multivariate normal distribution, from which we can 
# sample new values of the coefficients that are consistent with the fitted 
# model but will describe slightly different splines through the data and 
# consequently, slightly different trends.

library(MASS)
## simulate from posterior distribution of beta
Rbeta <- mvrnorm(n = 1000, coef(m2$gam), vcov(m2$gam))
Xp <- predict(m2$gam, newdata = pdat, type = "lpmatrix")
sim1 <- Xp %*% t(Rbeta)

## plot the observation and 25 of the 1000 trends
set.seed(321)
want <- sample(1000, 25)
ylim <- range(sim1[,want], gtemp$Annual)
plot(Annual ~ Year, data = gtemp, ylim = ylim, ylab = ylab)
matlines(pdat$Year, sim1[,want], col = "black", lty = 1, pch = NA)

# In summary , by using a model that is fitted to the entire period of data but 
# which can adapt to local features of the time series provides a powerful means 
# of estimating trends in temperature data. The thin-plate spline that describes 
# the fitted trend is defined by a set of coefficients that we can use to explore 
# the uncertainty in the model via simulation. Because the model can be 
# expressed as a linear mixed model we can exploit the lme() function to fit 
# correlation structures in the model residuals to account for the autocorrelation 
# in the data.