################################################################################
# Reproduction of Gavin Simpsons work on additive models from
# http://www.fromthebottomoftheheap.net/ 
# With a view to creating a process for analysis of Landsat veg time series
#
# Bart Huntley

rm(list=ls())

## this section "Modelling seasonal data with GAMs" from
## http://www.fromthebottomoftheheap.net/2014/05/09/modelling-seasonal-data-with-gam/

CET <- url("http://www.metoffice.gov.uk/hadobs/hadcet/cetml1659on.dat")
writeLines(readLines(CET, n = 10))

cet <- read.table(CET, sep = "", skip = 6, header = TRUE, 
                  fill = TRUE, na.string = c(-99.99, -99.9))
names(cet) <- c(month.abb, "Annual")
## remove last row of incomplete data
cet <- cet[-nrow(cet), ]
## get rid of the annual too - store for plotting
rn <- as.numeric(rownames(cet))
Years <- rn[1]:rn[length(rn)]
annCET <- data.frame(Temperature = cet[, ncol(cet)],
                     Year = Years) #annual
cet <- cet[, -ncol(cet)] #monthly


## Stack the data (long format required) - include new time variables
cet <- stack(cet)[,2:1]
names(cet) <- c("Month","Temperature")
## add in Year and nMonth for numeric month and a proper Date class
cet <- transform(cet, Year = (Year <- rep(Years, times = 12)),
                 nMonth = rep(1:12, each = length(Years)),
                 Date = as.Date(paste(Year, Month, "15", sep = "-"),
                                format = "%Y-%b-%d"))
## sort into temporal order
cet <- cet[with(cet, order(Date)), ] ## Really Important for plotting!
## Add in a Time variable
cet <- transform(cet, Time = as.numeric(Date) / 1000) #keep number small /1000
## check temporal order
head(cet)

## plot data
ylab <- expression(Temperature ~ (degree*C))
# annual
plot(Temperature ~ Year, data = annCET, type = "l",
     ylab = ylab, main = "CET")

## modelling time - fit naive model that assumes all observations are 
## independent. Other models will be tested against this
library(mgcv)
m <- gamm(Temperature ~ s(nMonth, bs = "cc", k = 12) + s(Time),
          data = cet)
# k set to max number of obs of months in year
# bs is basis type for spline setting - "cc" is cubic spline which will ensure
# no discontinuity between Dec and Jan

summary(m$gam)

## plot the crappy model - trend over time is too wiggly (overfitted)
layout(matrix(1:2, ncol = 2))
plot(m$gam, scale = 0) # remove scale arg to see same scale in each plot
layout(1)

## look at the residuals - obvious autocorrelation that model hasn't dealt with
## suggest an AR(p) model required
layout(matrix(1:2, ncol = 2))
acf(resid(m$lme), lag.max = 36, main = "ACF")
pacf(resid(m$lme), lag.max = 36, main = "pACF")
layout(1)

## fit 3 low order AR models
ctrl <- list(niterEM = 0, msVerbose = TRUE, optimMethod="L-BFGS-B")
# AR(1)
m1 <- gamm(Temperature ~ s(nMonth, bs = "cc", k = 12) + s(Time, k = 20),
           data = cet, correlation = corARMA(form = ~ 1|Year, p = 1),
           control = ctrl)
# AR(2)
m2 <- gamm(Temperature ~ s(nMonth, bs = "cc", k = 12) + s(Time, k = 20),
           data = cet, correlation = corARMA(form = ~ 1|Year, p = 2),
           control = ctrl)
# AR(3)               
m3 <- gamm(Temperature ~ s(nMonth, bs = "cc", k = 12) + s(Time, k = 20),
           data = cet, correlation = corARMA(form = ~ 1|Year, p = 3),
           control = ctrl)

## NOTES
# corARMA(form = ~ 1|Year, p = x) means fit an ARMA process to the residuals,
# where p indicates the order for the AR part of the ARMA model, and 
# form = ~ 1|Year means that the ARMA is nested within each year.  This 
# is potentially risky as we don’t consider residual variation from year to year.

## check models against the naive model
anova(m$lme, m1$lme, m2$lme, m3$lme)

#       Model df      AIC      BIC    logLik   Test   L.Ratio p-value
# m$lme      1  5 14978.99 15010.78 -7484.492                         
# m1$lme     2  6 14693.56 14731.72 -7340.780 1 vs 2 287.42398  <.0001
# m2$lme     3  7 14668.43 14712.95 -7327.217 2 vs 3  27.12703  <.0001
# m3$lme     4  8 14667.91 14718.79 -7325.954 3 vs 4   2.52623   0.112

# above reads AR(1) better than naive (1 vs 2) and even better is AR(2) as 
# further improvements 2 vs 3. Remember 1 is naive, 2 is AR(1) etc

## plot the AR(2) model
layout(matrix(1:2, ncol = 2))
plot(m2$gam, scale =0)
layout(1)

# Looking now at the normalized residuals (which take into account the 
# covariance matrix of the residuals). No significant autocorrelation
layout(matrix(1:2, ncol = 2))
res <- resid(m2$lme, type = "normalized")
acf(res, lag.max = 36, main = "ACF - AR(2) errors")
pacf(res, lag.max = 36, main = "pACF - AR(2) errors")
layout(1)


## extract model terms and plot
## we will use model to "predict" given 200 evenly spaced points along the 
## time variable
want <- seq(1, nrow(cet), length.out = 200) #index for 200
pdat <- with(cet,
             data.frame(Time = Time[want], Date = Date[want],
                        nMonth = nMonth[want]))

## predict trend contributions
p  <- predict(m$gam,  newdata = pdat, type = "terms", se.fit = TRUE)
p1 <- predict(m1$gam, newdata = pdat, type = "terms", se.fit = TRUE)
p2 <- predict(m2$gam, newdata = pdat, type = "terms", se.fit = TRUE)
p3 <- predict(m3$gam, newdata = pdat, type = "terms", se.fit = TRUE)

## combine with the predictions data, including fitted and SEs
pdat <- transform(pdat,
                  p  = p$fit[,2],  se  = p$se.fit[,2],
                  p1 = p1$fit[,2], se1 = p1$se.fit[,2],
                  p2 = p2$fit[,2], se2 = p2$se.fit[,2],
                  p3 = p3$fit[,2], se3 = p3$se.fit[,2])

## plot
op <- par(mar = c(5,4,2,2) + 0.1)
ylim <- with(pdat, range(p, p1, p2, p3))
ylim[1] <- floor(ylim[1])
ylim[2] <- ceiling(ylim[2])
ylab <- expression(Temperature ~ (degree*C ~ centred)) 
plot(Temperature - mean(Temperature) ~ Date, data = cet, type = "n", 
     ylab = ylab, ylim = ylim)
lines(p  ~ Date, data = pdat, col = "black")
lines(p1 ~ Date, data = pdat, col = "red")
lines(p2 ~ Date, data = pdat, col = "blue")
lines(p3 ~ Date, data = pdat, col = "forestgreen", lwd = 1)
legend("topleft", legend = c("Uncorrelated Errors", paste0("AR(", 1:3, ") Errors")),
       bty = "n", col = c("black","red","blue","forestgreen"),
       lty = 1, lwd = c(1,1,1))
par(op)

## The purpose of the above is to really get used to working with GAM's. 
## It illustrates nicely the effect of an ARMA model on un-correlated residuals 
## and the reduction of overfitting (wiggliness). What we really need to now look 
## at is changes in the seasonal temperature with the trend.