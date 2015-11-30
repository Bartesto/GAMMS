################################################################################
# Reproduction of Gavin Simpsons work on additive models from
# http://www.fromthebottomoftheheap.net/ 
# With a view to creating a process for analysis of Landsat veg time series
#
# Bart Huntley

rm(list=ls())

## this section "Climate change and spline interactions" from:
## http://www.fromthebottomoftheheap.net/2015/11/21/climate-change-and-spline-interactions/#fn1

source(con <- url("http://bit.ly/loadCET", method = "libcurl"))
close(con)
cet <- loadCET()
library("mgcv")
library("ggplot2")

knots <- list(nMonth = c(0.5, seq(1, 12, length = 10), 12.5))

# naive model
m0 <- gamm(Temperature ~ te(Year, nMonth, bs = c("cr","cc"), k = c(10,12)),
           data = cet, method = "REML", knots = knots)

plot(acf(resid(m0$lme, type = "normalized")))

# need AR(p) model so will fit them out to lag 8 with loop
ctrl <- list(niterEM = 0, optimMethod="L-BFGS-B", maxIter = 100, msMaxIter = 100)
for (i in 1:8) {
        m <- gamm(Temperature ~ te(Year, nMonth, bs = c("cr","cc"), k = c(10,12)),
                  data = cet, method = "REML", control = ctrl, knots = knots,
                  correlation = corARMA(form = ~ 1 | Year, p = i))
        assign(paste0("m", i), m) 
}

# A generalised likelihood ratio test can be used to test for which correlation 
# structure fits best
anova(m0$lme, m1$lme, m2$lme, m3$lme, m4$lme, m5$lme, m6$lme, m7$lme, m8$lme)

# Model df      AIC      BIC    logLik   Test   L.Ratio p-value
# m1$lme     1  6 14849.98 14888.13 -7418.988                         
# m2$lme     2  7 14836.78 14881.29 -7411.389 1 vs 2 15.197206  0.0001
# m3$lme     3  8 14810.73 14861.60 -7397.365 2 vs 3 28.047345  <.0001
# m4$lme     4  9 14784.63 14841.86 -7383.314 3 vs 4 28.101617  <.0001
# m5$lme     5 10 14778.35 14841.95 -7379.177 4 vs 5  8.275739  0.0040
# m6$lme     6 11 14776.49 14846.44 -7377.244 5 vs 6  3.865917  0.0493
# m7$lme     7 12 14762.45 14838.77 -7369.227 6 vs 7 16.032362  0.0001
# m8$lme     8 13 14764.33 14847.01 -7369.167 7 vs 8  0.119911  0.7291

# AR(7) is the winner (m7)

plot(acf(resid(m7$lme, type = "normalized")))

m <- m7

plot(m$gam, pers = TRUE)


# predict monthly temp for 1914 and 2014
pdat <- with(cet,
             data.frame(Year = rep(c(1914, 2014), each = 100),
                        nMonth = rep(seq(1, 12, length = 100), times = 2)))

pred <- predict(m$gam, newdata = pdat, se.fit = TRUE)
crit <- qt(0.975, df = df.residual(m$gam)) # ~95% interval critical t
pdat <- transform(pdat, fitted = pred$fit, se = pred$se.fit, fYear = as.factor(Year))
# adds variables containing the upper and lower pointwise confidence bounds (95%)
pdat <- transform(pdat,
                  upper = fitted + (crit * se),
                  lower = fitted - (crit * se))

p1 <- ggplot(pdat, aes(x = nMonth, y = fitted, group = fYear)) +
        geom_ribbon(mapping = aes(ymin = lower, ymax = upper,
                                  fill = fYear), alpha = 0.2) + # confidence band
        geom_line(aes(colour = fYear)) +    # predicted temperatures
        theme_bw() +                        # minimal theme
        theme(legend.position = "top") +    # push legend to the top
        labs(y = expression(Temperature ~ (degree*C)), x = NULL) +
        scale_fill_discrete(name = "Year") + # correct legend name
        scale_colour_discrete(name = "Year") +
        scale_x_continuous(breaks = 1:12,   # tweak where the x-axis ticks are
                           labels = month.abb, # & with what labels
                           minor_breaks = NULL)
p1


#Predict trends for each month, 1914–2014
pdat2 <- with(cet,
              data.frame(Year = rep(1914:2014, each = 12),
                         nMonth = rep(1:12, times = 101)))

pred2 <- predict(m$gam, newdata = pdat2, se.fit = TRUE)
## add predictions & SEs to the new data ready for plotting
pdat2 <- transform(pdat2,
                   fitted = pred2$fit,  # predicted values
                   se = pred2$se.fit,   # standard errors
                   fMonth = factor(month.abb[nMonth], # month as a factor
                                   levels = month.abb))
pdat2 <- transform(pdat2,
                   upper = fitted + (crit * se), # upper and...
                   lower = fitted - (crit * se)) # lower confidence bounds

p2 <- ggplot(pdat2, aes(x = Year, y = fitted, group = fMonth)) +
        geom_line(aes(colour = fMonth)) +   # draw trend lines
        theme_bw() +                        # minimal theme
        theme(legend.position = "none") +   # no legend
        labs(y = expression(Temperature ~ (degree*C)), x = NULL) +
        facet_wrap(~ fMonth, ncol = 6) +    # facet on month
        scale_y_continuous(breaks = seq(4, 17, by = 1),
                           minor_breaks = NULL) # nicer ticks
p2

# Predict trends for each month, 1914–2014, by quarter
pdat2$Quarter <- NA
pdat2$Quarter[pdat2$nMonth %in% c(12,1,2)] <- "Winter"
pdat2$Quarter[pdat2$nMonth %in% 3:5] <- "Spring"
pdat2$Quarter[pdat2$nMonth %in% 6:8] <- "Summer"
pdat2$Quarter[pdat2$nMonth %in% 9:11] <- "Autumn"
pdat2 <- transform(pdat2,
                   Quarter = factor(Quarter,
                                    levels = c("Spring","Summer","Autumn","Winter")))


p3 <- ggplot(pdat2, aes(x = Year, y = fitted, group = fMonth)) +
        geom_line(aes(colour = fMonth)) +   # draw trend lines
        theme_bw() +                        # minimal theme
        theme(legend.position = "top") +    # legend on top
        scale_fill_discrete(name = "Month") + # nicer legend title
        scale_colour_discrete(name = "Month") +
        labs(y = expression(Temperature ~ (degree*C)), x = NULL) +
        facet_grid(Quarter ~ ., scales = "free_y") # facet by Quarter
p3





























































