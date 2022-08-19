# analyse-ks-results.R
# Copyright (C) 2016, 2017 flossCoder
#
# This file is part of largeDevCorrelation.
#
# largeDevCorrelation is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# largeDevCorrelation is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

require("gplots")
require("minpack.lm")

#' Define a basic gaussian distribution with one peak.
#' 
#' @param x The abscissa value.
#' @param a The amplitude of the peak.
#' @param m The mean value of the gaussian distribution.
#' @param s The standard deviation of the gaussian distribution.
#' 
#' @return a*exp(-1/2*((x-m)/s)^2)
gaussian1 <- function(x, a, m, s) {
  return(a*exp(-1/2*((x-m)/s)^2))
}

#' Define a simple parabola f(x) = a * (x + b)^2 + c.
#' 
#' @param x The abscissa value.
#' @param a The amplitude and direction.
#' @param b The shift on abscissa direction.
#' @param c The shift on ordinate direction.
#' 
#' @return a * (x + b)^2 + c
par <- function(x, a, b, c) {
  return(a * (x + b)^2 + c)
}

#' Calculate the chi^2 / df for a given fit result.
#' 
#' @param fitResult The output of the fit.
#' 
#' @return The chi^2 / df statistics.
chiSqrdf <- function(fitResult) {
  chi2 <- sum(fitResult$m$resid()^2)
  df <- summary(fitResult)$df[2]
  return(chi2 / df)
}

#' Fit an peak gaussian to the given data set.
#' 
#' @param data The given data set (this one has to be analysed).
#' @param directory A valid directory, where the plot of the data set plus the two fits should be
#'        saved.
#' 
#' @return A list containing the statistical properties of the fit and the minimal position.
doFit <- function(data, directory) {
  # delete all non finite results
  redX <- data[((data[, 1] > 1.5) & (data[, 1] < 4) & is.finite(data[, 2])), 1]
  redY <- data[((data[, 1] > 1.5) & (data[, 1] < 4) & is.finite(data[, 2])), 2]
  
  # prepare starting conditions
  as <- min(redY)
  ms <- min(redX[(redY == as)])
  ss <- 0.5
  
  # do the one-peak gaussian fit
  gaussianFit <- nlsLM(redY ~ gaussian1(redX, a, m, s),
                       start = list(a = as, m = ms, s = ss),
                       control = nls.lm.control(maxiter = 1024))
  
  # do the parabola fit
  parFit <- nlsLM(redY ~ par(redX, a, b, c),
                  start = list(a = 1, b = as, c = ms),
                  control = nls.lm.control(maxiter = 1024))
  
  # plot the result and save it
  pdf(paste(directory, "ks-fit.pdf", sep = ""))
  plotCI(redX, redY, gap = 0) # plot histogram
  x <- c(1:10000) / 10000 * max(redX)
  auxG <- summary(gaussianFit)
  yG <- gaussian1(x, auxG$coefficients[1], auxG$coefficients[2], auxG$coefficients[3])
  lines(x, yG, col = "blue")
  points(x[yG == min(yG[is.finite(yG)])], min(yG[is.finite(yG)]), col = "blue")
  auxPar <- summary(parFit)
  yPar <- par(x, auxPar$coefficients[1], auxPar$coefficients[2], auxPar$coefficients[3])
  lines(x, yPar, col = "red")
  points(x[yPar == min(yPar[is.finite(yPar)])], min(yPar[is.finite(yPar)]), col = "red")
  dev.off()
  
  chi2dfgaussianFit <- chiSqrdf(gaussianFit)
  chi2dfparFit <- chiSqrdf(parFit)
  
  return(list(gaussianFit, chi2dfgaussianFit, x[yG == min(yG)], parFit, chi2dfparFit, x[yPar == min(yPar)]))
}

# load the data we want to analyse
data <- as.matrix(read.table(paste(directory, "/", param, "/", "ks-statistics.dat", sep = ""), sep = " ", header = FALSE))

# do the fitting business
result <- doFit(cbind(data[, 2], data[, 5]), paste(directory, "/", param, "/", sep = ""))

# save the R image
save.image(paste(directory, "/", param, "/", "ks-statistics.RData", sep = ""))

# add the result of the current configuration to the list of results
input <- strsplit(param, "_")[[1]]

numberOfVertices <- as.numeric(gsub("-", ".", input[length(input) - 2]))
cinter <- as.numeric(gsub("-", ".", input[length(input)]))

cat(sprintf("%s %s %s %s %s %s\n",
            toString(cinter),
            toString(format(result[[3]], digits=22)), # minimum cintra of log(p)-ks for gaussian fit
            toString(format(0.3 * summary(result[[1]])$coefficients[3], digits=22)), # error of fit
            toString(format(min(data[(data[, 5] != 0), 2]), digits=22)), # the minimum cintra, where log(p) != 0
            toString(format(data[data[, 2] == 1.5, 5], digits=22)), # the log(p)-value, where cintra = 1.5
            toString(format(min(data[data[, 4] < 0.05, 2]), digits=22))), # the cintra-value, where the p-value first drops below 0.05
    file = paste(directory, "/ks-result_", numberOfVertices, ".dat", sep = ""),
    append = TRUE)
