### Functions
# Code by Shannon Hennessey - s.hennessey13@gmail.com
# Data from the Hawaii Monitoring and Reporting Collaborative (HIMARC)

### Function to extract fitted GAMM diagnostics
# model = fitted model object
# resp.name = indicator
# dri.name = driver
# out.name = name of model
GAMM.summary <- function(model, resp.name, dri.name, out.name) {
  summary.gam <- as.data.frame(cbind(out.name,                     # Model name
                                     resp.name,                    # Response variable
                                     dri.name,                     # Pressure variable
                                     summary(model)$sp.criterion,  # GCV score
                                     summary(model)$dev.expl,      # deviance explained
                                     summary(model)$edf[1],        # driver estimated degrees of freedom (can select other elements of list too for additional columns if have random effects)
                                     summary(model)$s.pv[1]))      # driver p-value (can select other elements of list too for additional columns if have random effects)
}


### Function to plot fitted GAMMs, the bootstrapped confidence interval, and the threshold range(s) and point estimate(s)
# ci = dataframe with the driver values, predicted response, and upper and lower bounds of the CI and 2nd derivative CI
# figure.dir = directory to write figure file to
# myresponse = numeric sequence of indicator values
# mydriver = numeric sequence of driver values
# resp.name = name of the indicator for plotting (character)
# dri.name = name of the driver for plotting (character)
# thresh.min = the minimum driver value of the lowest threshold range
# thresh.max = the maximum driver value of the highest threshold range
# thresh1 = the driver value of the first threshold 'best estimate'
# thresh2 = the driver value of the second threshold 'best estimate', if applicable
gamm.threshold.fig <- function(ci, figure.dir, myresponse, mydriver, resp.name, dri.name,
                               thresh.min, thresh.max, thresh1, thresh2){
  # set file directory
  png(file = figure.dir, width = 170, height = 160, units = "mm", res = 600)
  
  plot(myresponse ~ mydriver,
       ylim = c(min(myresponse), max(myresponse)), 
       type = "n",
       ylab = "",
       xlab = "",
       axes = F)
  axis(1, tck = -0.02, labels = TRUE, cex = 0.25) # 250 for most, 75 for CAME and browsers
  gam.ax <- pretty(c(min(myresponse), max(myresponse), 3))
  axis(2, gam.ax, cex = 0.25)
  mtext(resp.name, side = 2, line = 2, cex = 1.1)
  mtext(dri.name, side = 1, line = 2, cex = 1.1)
  
  # add original data points
  points(myresponse ~ mydriver, pch = 16, cex = 0.55, lwd = 1, col = gray(0.7))

  # shade CIs
  polygon(c(ci$driver, rev(ci$driver)), 
          c(ci[,"upper"], rev(ci[,"lower"])),
          col = adjustcolor(gray(0.8), alpha.f = 0.65), 
          border = NA)
  
  # GAMM predicted response line
  lines(ci$response ~ ci$driver, lty = 1,  lwd = 2)
  
  # restrict to just within final threshold interval
  ci.final <- ci[which(ci$driver < thresh.max & ci$driver > thresh.min),] 
  
  # threshold arrow 
  arrows(x0 = thresh1, y0 = ci.final$response[which(round(ci.final$driver, digits = 4) == thresh1)], 
         x1 = thresh1, y1 = min(myresponse) + 0.01, #y1 = -5 for gamma, -0.1 binom
         length = 0.09, angle = 35, code = 2, col = "red", lty = 1, lwd = 2.5) # lwd?
  
  try(arrows(x0 = thresh2, y0 = ci.final$response[which(round(ci.final$driver, digits = 4) == thresh2)], 
             x1 = thresh2, y1 = min(myresponse) + 0.01, length = 0.09, angle = 35,
             code = 2, col = "red", lty = 1, lwd = 2.5), silent = TRUE)
  
  # derivative shading for threshold interval
  points(ci.final$response[ci[, "lower_ci2"] > 0] ~ ci.final$driver[ci[, "lower_ci2" ] > 0],
         col = gray(0.4), pch = 16, cex = 0.8)
  points(ci.final$response[ci[, "upper_ci2"] < 0] ~ ci.final$driver[ci[, "upper_ci2"] < 0],
         col = gray(0.4), pch = 16, cex = 0.8)
 
  # add box
  box()
  dev.off()
}
