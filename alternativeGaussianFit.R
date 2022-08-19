# alternativeGaussianFit.R
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

#' Define a gaussian distribution with two peaks.
#' 
#' @param x The abscissa value.
#' @param a1 The amplitude of the first peak.
#' @param m1 The mean value of the first gaussian distribution.
#' @param s1 The standard deviation of first the gaussian distribution.
#' @param a2 The amplitude of the second peak.
#' @param m2 The mean value of the secong gaussian distribution.
#' @param s2 The standard deviation of second the gaussian distribution.
#' 
#' @return a1*exp(-1/2*((x1-m1)/s1)^2) + a2*exp(-1/2*((x2-m2)/s2)^2)
gaussian2 <- function(x, a1, m1, s1, a2, m2, s2) {
  return(gaussian1(x, a1, m1, s1) + gaussian1(x, a2, m2, s2))
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

#' Fit an one peak gaussian to the given histogram using the errorbars of the histogram as weights.
#' 
#' @param histogram The given histogram (this one has to be analysed).
#' 
#' @return A list containing the statistical properties of the fit.
doGaussianFit <- function(histogram) {
  # prepare starting conditions
  as1 <- max(histogram[, 2])
  ms1 <- histogram[(histogram[, 2] == as1), 1]
  ss1 <- 0.5
  
  # do the one-peak gaussian fit
  gaussianFit1 <- nlsLM(histogram[, 2] ~ gaussian1(histogram[, 1], a, m, s),
                        start = list(a = as1, m = ms1, s = ss1),
                        weights = histogram[, 3],
                        control = nls.lm.control(maxiter = 1024))
  
  chi2dfgaussianFit1 <- chiSqrdf(gaussianFit1)
  
  return(list(gaussianFit1, chi2dfgaussianFit1))
}

#' Fit an one peak gaussian and a two peak gaussian to the given histogram using the errorbars of the
#' histogram as weights.
#' 
#' @param histogram The given histogram (this one has to be analysed).
#' @param directory A valid directory, where the plot of the histogram data plus the two fits should be
#'        saved.
#' 
#' @return A list containing the statistical properties of the two fits and the difference matrix with the first fit.
doFitsOfDifferents <- function(histogram, directory) {
  firstFit <- doGaussianFit(histogram)
  # calculate the difference of the fit and the data
  firstFitPrediction <- predict(firstFit[[1]])
  difference <- histogram
  difference[, 2] <- histogram[, 2] - firstFitPrediction
  secondFit <- doGaussianFit(difference)
  secondFitPrediction <- predict(secondFit[[1]])
  
  pdf(paste(directory, "test-gaussian-fit.pdf", sep = ""))
  ymin = min(histogram[, 2], firstFitPrediction, difference[, 2], secondFitPrediction)
  ymax = max(histogram[, 2], firstFitPrediction, difference[, 2], secondFitPrediction)
  plotCI(histogram[, 1], histogram[, 2], uiw = histogram[, 3], ylim = c(ymin, ymax), gap = 0) # plot histogram
  lines(histogram[, 1], firstFitPrediction, col = "blue")
  lines(histogram[, 1], difference[, 2], col = "red")
  lines(histogram[, 1], secondFitPrediction, col = "green")
  legend("topright",
         c("data", "gaussian fit", "difference", "difference fit"),
         col = c("black", "blue", "red", "green"),
         lwd = c(1, 1, 1, 1),
         lty = c(NA, 1, 1, 1),
         pch = c(1, NA, NA, NA))
  dev.off()
  
  return(list(firstFit, difference, secondFit))
}

# for none scientific file processing => needed for file names
options(scipen = 999)

dataMetropolis <- NA
dataWangLandau <- NA

source(arg)

# open the resulting histogram
histogram <- as.matrix(read.table(paste(directory, "hist_", numberOfVertices, "_result.dat", sep = ""),
                                  sep = " ", header = FALSE))

# make sure, the histogram is normalized
histogram[, 2] <- histogram[, 2] / sum(histogram[, 2])

# do the fitting business
fitResult <- doFitsOfDifferents(histogram, directory)

# Save all results as R workspace:
save.image(paste(directory, "ana-gaussianTest.RData", sep = ""))
