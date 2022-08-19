# gaussianFits.R
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

#' Fit an one peak gaussian and a two peak gaussian to the given histogram using the errorbars of the
#' histogram as weights.
#' 
#' @param histogram The given histogram (this one has to be analysed).
#' @param directory A valid directory, where the plot of the histogram data plus the two fits should be
#'        saved.
#' 
#' @return A list containing the statistical properties of the two fits.
doFits <- function(histogram, directory) {
  # prepare starting conditions
  as1 <- max(histogram[, 2])
  ms1 <- histogram[(histogram[, 2] == as1), 1]
  ss1 <- 0.5
  ms2 <- NA
  as2 <- NA
  ss2 <- 0.8
  
  # do the one-peak gaussian fit
  gaussianFit1 <- nlsLM(histogram[, 2] ~ gaussian1(histogram[, 1], a, m, s),
                        start = list(a = as1, m = ms1, s = ss1),
                        weights = histogram[, 3],
                        control = nls.lm.control(maxiter = 1024))
  
  gaussianFit2 <- NA
  chi2dfgaussianFit2 <- Inf
  # do the two-peak gaussian fit
  for (ms2 in c(1:dim(histogram)[1])) {
    # try different ms2
    as2 <- histogram[ms2, 2]
    
    tryCatch(
      {
        gaussianFit2C <- nlsLM(histogram[, 2] ~ gaussian2(histogram[, 1], a1, m1, s1, a2, m2, s2),
                               start = list(a1 = as1, m1 = ms1, s1 = ss1, a2 = as2, m2 = ms2, s2 = ss2),
                               weights = histogram[, 3],
                               control = nls.lm.control(maxiter = 1024))
        chi2dfgaussianFit2C <- chiSqrdf(gaussianFit2C)
        if (chi2dfgaussianFit2C < chi2dfgaussianFit2) {
          # the current fit is the best => save it
          gaussianFit2 <- gaussianFit2C
          chi2dfgaussianFit2 <- chi2dfgaussianFit2C
        }
      }, error = function(e) print(paste(ms2, "failed"))
    )
  }
  
  # plot the result and save it
  pdf(paste(directory, "gaussian-fits.pdf", sep = ""))
  plotCI(histogram[, 1], histogram[, 2], uiw = histogram[, 3], gap = 0, legend.text = "data") # plot histogram
  lines(histogram[, 1], predict(gaussianFit1), col = "blue", legend.text = "gaussian 1-peak")
  lines(histogram[, 1], predict(gaussianFit2), col = "red", legend.text = "gaussian 2-peak")
  #legend("best")
  dev.off()
  
  chi2dfgaussianFit1 <- chiSqrdf(gaussianFit1)
  
  return(list(gaussianFit1, gaussianFit2, chi2dfgaussianFit1, chi2dfgaussianFit2))
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

input <- strsplit(param, "_")[[1]]

cinter <- as.numeric(gsub("-", ".", input[length(input) - 1]))
cintra <- as.numeric(gsub("-", ".", input[length(input)]))

# do the fitting business
fitResult <- doFits(histogram, directory)

# Save all results as R workspace:
save.image(paste(directory, "ana-gaussianFits.RData", sep = ""))
