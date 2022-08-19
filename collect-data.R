# collect-data.R
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

require("minpack.lm")

#' Open a certain file.
#' 
#' @param directory Where to find the data?
#' @param filename Which file should be opened?
#' 
#' @return matrix containing the data.
openData <- function(directory, filename) {
  return(as.matrix(read.table(paste(directory, "/", filename, sep=""), sep = " ", header = FALSE)))
}

#' Save the given data.
#' 
#' @param directory Where the data have to be saved.
#' @param filename Which should the file have? Be aware of overwriting existing files.
#' @param data Matrix with stuff to save.
saveData <- function(directory, filename, data) {
  write.table(data, file = paste(directory, "/", filename, sep=""), sep = " ", row.names = FALSE, col.names = FALSE)
}

#' Calculate exponential growth.
#' 
#' @param x value
#' @param a parameter
#' @param b parameter
#' @param c parameter
#' 
#' @return y value (a * x^(-b) * c)
expGrowth <- function(x, a, b, c) {
  return(a * x^(-b) + c)
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

# Read the data for all system sizes and save the minimal log(p) plus a fit.
sizes <- c(18, 25, 35, 40, 50, 60, 70, 75, 100)
directory <- "/home/stephan/Schreibtisch/ks-results"
configs <- c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)

data <- list()
fitResults <- list()
fitChiSqrdf <- list()

# Open all required data.
for (size in c(1:length(sizes))) {
  data[[size]] <- openData(directory, paste("ks-result_", sizes[size], ".dat", sep = ""))
}

# Do the main work for all configurations and system sizes.
for (config in c(1:length(configs))) {
  minLogP <- matrix(0, length(sizes), 3)
  for (size in c(1:length(sizes))) {
    minLogP[size, 1] <- sizes[size]
    currentData <- data[[size]]
    minLogP[size, 2] <- currentData[currentData[, 1] == configs[config], 2]
    minLogP[size, 3] <- currentData[currentData[, 1] == configs[config], 3]
  }
  
  # save the data file
  saveData(directory, paste("minLogP_", configs[config], ".dat", sep = ""), minLogP)
  
  # do the fitting business
  aI <- 0
  bI <- 0
  cI <- 0
  fitResults[[size]] <- nlsLM(minLogP[, 2] ~ expGrowth(minLogP[, 1], a, b, c),
                              start = list(a = aI, b = bI, c = cI),
                              control = nls.lm.control(maxiter = 1024))
  
  # save the fit results as text files
  x <- c(0:1000)/1000 * (max(minLogP[,1]) - min(minLogP[, 1])) + min(minLogP[, 1])
  coef <- summary(fitResults[[size]])$coefficients
  y <- expGrowth(x, coef[1], coef[2], coef[3])
  
  saveData(directory, paste("minLogP_fit_", configs[config], ".dat", sep = ""), cbind(x, y))
  
  # calculate chiSqrdf
  fitChiSqrdf[[size]] <- chiSqrdf(fitResults[[size]])
}

save.image(paste(directory, "/minLogP-fits.RData", sep = ""))
