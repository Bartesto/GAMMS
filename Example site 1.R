rm(list=ls())

library(mgcv)
library(nlme)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(gridExtra)
# Functions
plotForecastErrors <- function(forecasterrors)
{
        # make a histogram of the forecast errors:
        mybinsize <- IQR(forecasterrors)/4
        mysd   <- sd(forecasterrors)
        mymin  <- min(forecasterrors) - mysd*5
        mymax  <- max(forecasterrors) + mysd*3
        # generate normally distributed data with mean 0 and standard deviation mysd
        mynorm <- rnorm(10000, mean=0, sd=mysd)
        mymin2 <- min(mynorm)
        mymax2 <- max(mynorm)
        if (mymin2 < mymin) { mymin <- mymin2 }
        if (mymax2 > mymax) { mymax <- mymax2 }
        # make a red histogram of the forecast errors, with the normally distributed data overlaid:
        mybins <- seq(mymin, mymax, mybinsize)
        hist(forecasterrors, col="red", freq=FALSE, breaks=mybins)
        # freq=FALSE ensures the area under the histogram = 1
        # generate normally distributed data with mean 0 and standard deviation mysd
        myhist <- hist(mynorm, plot=FALSE, breaks=mybins)
        # plot the normal curve as a blue line on top of the histogram of forecast errors:
        points(myhist$mids, myhist$density, type="l", col="blue", lwd=2)
}

setwd("Z:\\DOCUMENTATION\\BART\\R\\R_DEV\\GAMMS")        
# Deal with data
dveg <- read.csv("GAMM_month_interp_mtsd.csv", header = TRUE, stringsAsFactors = FALSE)
dveg <- dveg[,-1] # remove 1st column of row numbers

# Create an individual data frame for each site
for(i in 1:length(names(dveg))-1){
        d <- dveg[, c(1, 1 + i)]
        d$date <- as.Date(d$date, "%Y-%m-%d")
        name <- names(d)[2]
        d$site <- name
        names(d)[2] <- "veg"
        assign(paste0("d_", name), d)
}

df <- d_dhi_01
df <- df %>%
        mutate(nMonth = month(date), Month = month.abb[nMonth], 
               Year = year(date), Time = as.numeric(date) / 100)
df <- df[with(df, order(date)), ]

# Initial plots
p1 <- ggplot(df, aes(x = date, y = veg)) +
        geom_point() +
        geom_line() +
        theme_bw() +
        ggtitle("Site dhi_01") +
        ylab("vegetation %") +
        ylim(c(0, 80))
p1
ylim <- with(df, range(veg))
p2 <- ggplot(df, aes(x = date, y = veg)) +
        geom_point() +
        geom_line() +
        theme_bw() +
        ggtitle("Site dhi_01") +
        ylab("vegetation %") +
        ylim(ylim)
p2

# Naive model
m <- gamm(veg ~ s(nMonth, bs = "cc", k = 12) + s(Time),
          data = df)
summary(m$gam)

layout(matrix(1:4, ncol = 2))
plot(m$gam, scale = 0)
acf(resid(m$lme), lag.max = 36, main = "ACF")
pacf(resid(m$lme), lag.max = 36, main = "pACF")
layout(1)

# accounting for correlations
ctrl <- list(niterEM = 0, msVerbose = FALSE, optimMethod="L-BFGS-B")
# AR(1)
m1 <- gamm(veg ~ s(nMonth, bs = "cc", k = 12) + s(Time, k = 20),
           data = df, correlation = corARMA(form = ~ 1|Year, p = 1),
           control = ctrl)

# AR(2)
m2 <- gamm(veg ~ s(nMonth, bs = "cc", k = 12) + s(Time, k = 20),
           data = df, correlation = corARMA(form = ~ 1|Year, p = 2),
           control = ctrl)

# AR(3)
m3 <- gamm(veg ~ s(nMonth, bs = "cc", k = 12) + s(Time, k = 20),
           data = df, correlation = corARMA(form = ~ 1|Year, p = 3),
           control = ctrl)

anova(m$lme, m1$lme, m2$lme, m3$lme)

layout(matrix(1:4, ncol = 2))
plot(m1$gam, scale = 0)
res <- resid(m1$lme, type = "normalized")
acf(res, lag.max = 36, main = "ACF - AR(1) errors")
pacf(res, lag.max = 36, main = "pACF- AR(1) errors")
layout(1)

plotForecastErrors(resid(m1$lme))

# Model looks good move on and forecast for existing time period
want <- seq(1, nrow(df), length.out = 200) # evenly spaced sequence of points
pdat <- with(df,
             data.frame(Time = Time[want], Date = date[want],
                        nMonth = nMonth[want])) # new predictions data frame
p1 <- predict(m1$gam, newdata = pdat, type = "terms", se.fit = TRUE)
## combine with the predictions data, including fitted and SEs
pdat <- transform(pdat, p1 = p1$fit[,2], se1 = p1$se.fit[,2])

ggplot(df, aes(x = date, y = veg - mean(veg))) +
        geom_point() +
        geom_line() +
        theme_bw() +
        geom_line(data = pdat, aes(x = Date, y = p1), stat = "identity", 
                  colour = "blue")
summary(m1$gam)

# statistically significant changes in the time series
df.res <- df.residual(m1$gam)
crit.t <- qt(0.025, df.res, lower.tail = FALSE)
pdat <- transform(pdat,
                  upper = p1 + (crit.t * se1),
                  lower = p1 - (crit.t * se1))

tmpf <- tempfile()
download.file("https://gist.github.com/gavinsimpson/e73f011fdaaab4bb5a30/raw/82118ee30c9ef1254795d2ec6d356a664cc138ab/Deriv.R",
              tmpf)
source(tmpf)


Term <- "Time"
m1.d <- Deriv(m1)
m1.dci <- confint(m1.d, term = Term)
m1.dsig <- signifD(pdat$p1, d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper, m1.dci[[Term]]$lower)

ggplot(df, aes(x = date, y = veg - mean(veg))) +
        geom_point(colour = "light grey") +
        geom_line(colour = "light grey") +
        theme_bw() +
        geom_line(data = pdat, aes(x = Date, y = p1), stat = "identity", 
                  colour = "black") +
        geom_line(data = pdat, aes(x = Date, y = upper), stat = "identity", 
                  colour = "black", linetype = 2) +
        geom_line(data = pdat, aes(x = Date, y = lower), stat = "identity", 
                  colour = "black", linetype = 2) +
        geom_line(data = pdat, aes(x = Date, y = unlist(m1.dsig$incr)), stat = "identity", 
                  colour = "blue", size = 1.5) +
        geom_line(data = pdat, aes(x = Date, y = unlist(m1.dsig$decr)), stat = "identity", 
                  colour = "red", size = 1.5) 
ggsave("site01_sig_chng_ts.jpeg")

## percentage veg adaption
y_mod <- pdat$p1 + mean(df$veg)
y_up <- pdat$upper + mean(df$veg)
y_low <- pdat$lower + mean(df$veg)


ggplot(df, aes(x = date, y = veg)) +
        geom_point(colour = "dark grey") +
        geom_line(colour = "dark grey") +
        theme_bw() +
        geom_line(data = pdat, aes(x = Date, y = y_mod), stat = "identity", 
                  colour = "black") +
        geom_line(data = pdat, aes(x = Date, y = y_up), stat = "identity", 
                  colour = "black", linetype = 2) +
        geom_line(data = pdat, aes(x = Date, y = y_low), stat = "identity", 
                  colour = "black", linetype = 2) +
        geom_line(data = pdat, aes(x = Date, y = unlist(m1.dsig$incr) + mean(df$veg)), stat = "identity", 
                  colour = "blue", size = 1.5) +
        geom_line(data = pdat, aes(x = Date, y = unlist(m1.dsig$decr) + mean(df$veg)), stat = "identity", 
                  colour = "red", size = 1.5) +
        ylim(0, 80)

ggsave("site01_sig_chng_ts_origscale.jpeg")

# Spline interactions

knots <- list(nMonth = c(0.5, seq(1, 12, length = 10), 12.5))
nm <- gamm(veg ~ te(Year, nMonth, bs = c("cr","cc"), k = c(10,12)),
           data = df, method = "REML", knots = knots)
plot(acf(resid(nm$lme, type = "normalized")))

ctrl <- list(niterEM = 0, optimMethod="L-BFGS-B", maxIter = 100, msMaxIter = 100)

for (i in 1:4) {
        m <- gamm(veg ~ te(Year, nMonth, bs = c("cr","cc"), k = c(10,12)),
                   data = df, method = "REML", control = ctrl, knots = knots,
                   correlation = corARMA(form = ~ 1 | Year, p = i))
        assign(paste0("nm", i), m) 
}

anova(nm$lme, nm1$lme, nm2$lme, nm3$lme, nm4$lme)

plot(acf(resid(nm1$lme, type = "normalized")))
summary(nm1$gam)

# monthly
npdat <- with(df,
             data.frame(Year = rep(c(1990, 2014), each = 100),
                        nMonth = rep(seq(1, 12, length = 100), times = 2)))

pred <- predict(nm1$gam, newdata = npdat, se.fit = TRUE)
crit <- qt(0.975, df = df.residual(nm1$gam)) # ~95% interval critical t
npdat <- transform(npdat, fitted = pred$fit, se = pred$se.fit, fYear = as.factor(Year))
npdat <- transform(npdat,
                  upper = fitted + (crit * se),
                  lower = fitted - (crit * se))

p3 <- ggplot(npdat, aes(x = nMonth, y = fitted, group = fYear)) +
        geom_ribbon(mapping = aes(ymin = lower, ymax = upper,
                                  fill = fYear), alpha = 0.2) + # confidence band
        geom_line(aes(colour = fYear)) +    # predicted temperatures
        theme_bw() +                        # minimal theme
        theme(legend.position = "top") +    # push legend to the top
        labs(y = "vegetation cover %", x = NULL) +
        scale_fill_discrete(name = "Year") + # correct legend name
        scale_colour_discrete(name = "Year") +
        scale_x_continuous(breaks = 1:12,   # tweak where the x-axis ticks are
                           labels = month.abb, # & with what labels
                           minor_breaks = NULL)
p3
ggsave("site01_1990-2014.jpeg")

