# kstest.R
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

require("caTools")

# for none scientific file processing => needed for file names
options(scipen = 999)

dataMetropolis <- NA
dataWangLandau <- NA

source(paste(dirname(sys.frame(1)$ofile), "/rate.R", sep = ""))

source(arg)

#' This function fills up the given histogram in case the extreme large deviation tails have not been sampled.
#' 
#' @param histogram The given histogram.
#' @param numberOfVertices The number of vertices of the corresponding graph.
#' 
#' @return The filled up histogram.
fillUpZeroValues <- function(histogram, numberOfVertices) {
  if (dim(histogram)[1] == numberOfVertices) {
    return(histogram) # everything is fine, nothing to do
  }
  
  # Repair tail around one, if required.
  if (histogram[1, 1] > 1) {
    aux <- matrix(0, nrow = histogram[1, 1] - 1, ncol = 3)
    aux[, 1] <- c(1:(histogram[1, 1] - 1))
    histogram <- rbind(aux, histogram)
  }
  
  # Repair tail around the max size, if required.
  if (histogram[dim(histogram)[1], 1] < numberOfVertices) {
    aux <- matrix(0, nrow = numberOfVertices - histogram[dim(histogram)[1], 1], ncol = 3)
    aux[, 1] <- c((histogram[dim(histogram)[1], 1] + 1):numberOfVertices)
    histogram <- rbind(histogram, aux)
  }
  
  return(histogram)
}

#' This function performs a Kolmogorov-Smirnov test and saves the result.
#' 
#' @param experiment The experimental data of the SBM.
#' @param er The result for the corresponding ER graphs.
#' @param name A sensefull name for saving the results.
#' 
#' @return The result of the KS-test
calulateKSTest <- function(experiment, er, name) {
  res <- ks.test(experiment, er, alternative = "two.sided")
  writeLines(capture.output(print(res)), con = paste(directory, "ana-ks_", name, ".ks", sep = ""))
  
  return(res)
}

#' Calculate the distance between the two given data sets. This function even works in cases the x values
#' are not exactly identical.
#' 
#' @param er The result for the corresponding ER graphs.
#' @param experiment The experimental data of the SBM.
#' 
#' @return The distance.
distance <- function(er, experiment) {
  result <- matrix(0, nrow = dim(experiment)[1], ncol = 2)
  result[, 1] <- experiment[, 1]
  
  for (i in 1:dim(experiment)[1]) {
    # obtain the best fitting index of the x-coordinate 
    difference <- abs(er[, 1] - experiment[i, 1])
    minimum <- min(difference)
    index <- which(difference == minimum, arr.ind = TRUE)
    
    # calculate the difference of the y-coordinate
    result[i, 2] <- (er[index, 2] - experiment[i, 2])
  }
  
  return(result)
}

#' Calculate the cumulative distribution function (CDF) = integral of -inf to x over histogram(x)
#' of the given PDF.
#' 
#' @param histogram Given PDF.
#' 
#' @return The CDF of the given PDF.
calculateCDF <- function(histogram) {
  cdf <- matrix(0, nrow = dim(histogram)[1], ncol = 3)
  cdf[, 1] <- histogram[, 1]
  cdf[1, 2] <- histogram[1, 2]
  cdf[1, 3] <- histogram[1, 2]
  for (i in c(2:dim(histogram)[1])) {
    cdf[i, 2] <- cdf[(i - 1), 2] + histogram[i, 2]
    cdf[i, 3] <- cdf[(i - 1), 3] + histogram[i, 3]
  }
  return(cdf)
}

input <- strsplit(param, "_")[[1]] # dirER contains simulation results for the corresponding ER graphs

cinter <- as.numeric(gsub("-", ".", input[length(input) - 1]))
cintra <- as.numeric(gsub("-", ".", input[length(input)]))

c <- (cinter + cintra) / 2.0 # calculate the corresponding er connectivity

erAna <- doCalculation(numberOfVertices, c, 20)[[2]] # calculate the analyticaly rate function

# load all required data

# load the experimental rate function of the stochastic blockmodel
data <- as.matrix(read.table(paste(directory, "hist_", numberOfVertices, "_rate-function.dat", sep = ""),
                             sep = " ", header = FALSE))

# load the experimental rate function of the corresponding ER graph
er <- as.matrix(read.table(paste(dirER, "/er_non_digraph_", numberOfVertices, "_", gsub("\\.", "-", c),
                                   "/", "hist_", numberOfVertices, "_rate-function.dat", sep = ""),
                             sep = " ", header = FALSE))

# load the experimental result of the stochastic blockmodel
dataResult <- as.matrix(read.table(paste(directory, "hist_", numberOfVertices, "_result.dat", sep = ""),
                                   sep = " ", header = FALSE))

# fill up not determined large deviation tail with zero (if required)
dataResult <- fillUpZeroValues(dataResult, numberOfVertices)

# load the experimental result of the corresponding ER graph
erResult <- as.matrix(read.table(paste(dirER, "/er_non_digraph_", numberOfVertices, "_", gsub("\\.", "-", c),
                                       "/", "hist_", numberOfVertices, "_result.dat", sep = ""),
                                 sep = " ", header = FALSE))

# fill up not determined large deviation tail with zero (if required)
erResult <- fillUpZeroValues(erResult, numberOfVertices)

# normalize for security reasons (to prevent numerically errors)
dataResult[, 2] <- dataResult[, 2] / sum(dataResult[, 2])
erResult[, 2] <- erResult[, 2] / sum(erResult[, 2])

# load the analytical results
#phi <- as.matrix(read.table(ana, sep = " ", header = FALSE))
# the analytical results are NOT normalized => normalize them
#phi[, 2] <- phi[, 2] / sum(phi[, 2])
# save the normalized results
#write.table(phi, file = paste(directory, "phi.dat", sep = ""), row.names = FALSE, col.names = FALSE)

# KS-Test

# calculate the CDF for the SBM and the ER distripution

dataCDF <- calculateCDF(dataResult)
erCDF <- calculateCDF(erResult)

# make the KS test

#ksResult <- calulateKSTest(dataCDF, erCDF, param)
ksResult <- calulateKSTest(dataResult, erResult, param)

# compare experimental rate functions

# calculate the difference between the two rate functions
difference <- er
difference[, 2] <- er[, 2] - data[, 2]

# compare experimental probability functions

# calculate the difference between the two rate functions
differenceResult <- erResult
differenceResult[, 2] <- erResult[, 2] - dataResult[, 2]

# plot the two rate functions
pdf(paste(directory, "rate-functions.pdf", sep = ""))
plot(er, xlab = "s", ylab = "er(s)",
     xlim = c(min(er[is.finite(er[, 1]), 1], data[is.finite(data[, 1]), 1]), max(er[is.finite(er[, 1]), 1], data[is.finite(data[, 1]), 1])),
     ylim = c(min(er[is.finite(er[, 2]), 2], data[is.finite(data[, 2]), 2]), max(er[is.finite(er[, 2]), 2], data[is.finite(data[, 2]), 2])),
     col = "red")
points(data, col = "blue")
dev.off()

# plot the two rate functions plut the difference
pdf(paste(directory, "rate-difference.pdf", sep = ""))
plot(er, xlab = "s", ylab = "er(s)",
     xlim = c(min(er[is.finite(er[, 1]), 1], data[is.finite(data[, 1]), 1], difference[is.finite(difference[, 1]), 1]), max(er[is.finite(er[, 1]), 1], data[is.finite(data[, 1]), 1], difference[is.finite(difference[, 1]), 1])),
     ylim = c(min(er[is.finite(er[, 2]), 2], data[is.finite(data[, 2]), 2], difference[is.finite(difference[, 2]), 2]), max(er[is.finite(er[, 2]), 2], data[is.finite(data[, 2]), 2], difference[is.finite(difference[, 2]), 2])),
     col = "red")
points(data, col = "blue")
points(difference, col = "green")
dev.off()

pdf(paste(directory, "rate-only-difference.pdf", sep = ""))
plot(difference, xlab = "s", ylab = "er(s)")
dev.off()

# plot the two result functions
pdf(paste(directory, "result-functions.pdf", sep = ""))
plot(erResult, xlab = "s", ylab = "P(s)",
     xlim = c(min(erResult[is.finite(erResult[, 1]), 1], dataResult[is.finite(dataResult[, 1]), 1]), max(erResult[is.finite(erResult[, 1]), 1], dataResult[is.finite(dataResult[, 1]), 1])),
     ylim = c(min(erResult[is.finite(erResult[, 2]), 2], dataResult[is.finite(dataResult[, 2]), 2]), max(erResult[is.finite(erResult[, 2]), 2], dataResult[is.finite(dataResult[, 2]), 2])),
     col = "red")
points(dataResult, col = "blue")
dev.off()

# plot the two result functions plut the difference
pdf(paste(directory, "result-difference.pdf", sep = ""))
plot(erResult, xlab = "s", ylab = "P(s)",
     xlim = c(min(erResult[is.finite(erResult[, 1]), 1], dataResult[is.finite(dataResult[, 1]), 1], differenceResult[is.finite(differenceResult[, 1]), 1]), max(erResult[is.finite(erResult[, 1]), 1], dataResult[is.finite(dataResult[, 1]), 1], differenceResult[is.finite(differenceResult[, 1]), 1])),
     ylim = c(min(erResult[is.finite(erResult[, 2]), 2], dataResult[is.finite(dataResult[, 2]), 2], differenceResult[is.finite(differenceResult[, 2]), 2]), max(erResult[is.finite(erResult[, 2]), 2], dataResult[is.finite(dataResult[, 2]), 2], differenceResult[is.finite(differenceResult[, 2]), 2])),
     col = "red")
points(dataResult, col = "blue")
points(differenceResult, col = "green")
dev.off()

# Calculate the sum of the differences, the sum of the absolut differences, the mean of the differences
# and the mean of the absolut differences.
sumDiff <- sum(difference[, 2] * is.finite(difference[, 2]), na.rm = TRUE)
sumAbsDiff <- sum(abs(difference[, 2]) * is.finite(difference[, 2]), na.rm = TRUE)
meanDiff <- mean(difference[, 2] * is.finite(difference[, 2]), na.rm = TRUE)
meanAbsDiff <- sum(abs(difference[, 2]) * is.finite(difference[, 2]), na.rm = TRUE)

# Calculate the sum of the differences, the sum of the absolut differences, the mean of the differences
# and the mean of the absolut differences.
sumDiffResult <- sum(differenceResult[, 2] * is.finite(differenceResult[, 2]), na.rm = TRUE)
sumAbsDiffResult <- sum(abs(differenceResult[, 2]) * is.finite(differenceResult[, 2]), na.rm = TRUE)
meanDiffResult <- mean(differenceResult[, 2] * is.finite(differenceResult[, 2]), na.rm = TRUE)
meanAbsDiffResult <- sum(abs(differenceResult[, 2]) * is.finite(differenceResult[, 2]), na.rm = TRUE)

# Calculate the integrals of the two data sets.
integralData <- sum(data[, 2] * is.finite(data[, 2]), na.rm = TRUE)
integralER <- sum(er[, 2] * is.finite(er[, 2]), na.rm = TRUE)

# Calculate the integrals of the two data sets.
integralDataResult <- sum(dataResult[, 2] * is.finite(dataResult[, 2]), na.rm = TRUE)
integralERResult <- sum(erResult[, 2] * is.finite(erResult[, 2]), na.rm = TRUE)

# Calculate the difference of the two integrals
integralDifference <- (integralData - integralER)
absIntegralDifference <- abs(integralDifference)

# Calculate the difference of the two integrals
integralDifferenceResult <- (integralDataResult - integralERResult)
absIntegralDifferenceResult <- abs(integralDifferenceResult)

# compare experimental sbm with analytical er

# calculate the difference between the two rate functions
differenceAna <- distance(erAna, data)

# plot the two rate functions
pdf(paste(directory, "rate-functions_ana.pdf", sep = ""))
plot(er, xlab = "s", ylab = "er(s)",
     xlim = c(min(er[is.finite(er[, 1]), 1], data[is.finite(data[, 1]), 1]), max(er[is.finite(er[, 1]), 1], data[is.finite(data[, 1]), 1])),
     ylim = c(min(er[is.finite(er[, 2]), 2], data[is.finite(data[, 2]), 2]), max(er[is.finite(er[, 2]), 2], data[is.finite(data[, 2]), 2])),
     col = "red")
points(data, col = "blue")
dev.off()

# plot the two rate functions plut the difference
pdf(paste(directory, "rate-difference_ana.pdf", sep = ""))
plot(er, xlab = "s", ylab = "er(s)",
     xlim = c(min(erAna[is.finite(erAna[, 1]), 1], data[is.finite(data[, 1]), 1]), max(erAna[is.finite(erAna[, 1]), 1], data[is.finite(data[, 1]), 1])),
     ylim = c(min(erAna[is.finite(erAna[, 2]), 2], data[is.finite(data[, 2]), 2]), max(erAna[is.finite(erAna[, 2]), 2], data[is.finite(data[, 2]), 2])),
     col = "red")
points(data, col = "blue")
points(difference, col = "green")
dev.off()

# Calculate the sum of the differences, the sum of the absolut differences, the mean of the differences
# and the mean of the absolut differences.
sumDiffAna <- sum(differenceAna[, 2] * is.finite(differenceAna[, 2]), na.rm = TRUE)
sumAbsDiffAna <- sum(abs(differenceAna[, 2]) * is.finite(differenceAna[, 2]), na.rm = TRUE)
meanDiffAna <- mean(differenceAna[, 2] * is.finite(differenceAna[, 2]), na.rm = TRUE)
meanAbsDiffAna <- sum(abs(differenceAna[, 2]) * is.finite(differenceAna[, 2]), na.rm = TRUE)

# Calculate the integrals of the two data sets.
integralERAna <- sum(erAna[, 2] * is.finite(erAna[, 2]), na.rm = TRUE)

# Calculate the difference of the two integrals
integralDifferenceAna <- (integralData - integralERAna)
absIntegralDifferenceAna <- abs(integralDifferenceAna)

# Save all results as R workspace:
save.image(paste(directory, "ana-ks.RData", sep = ""))
