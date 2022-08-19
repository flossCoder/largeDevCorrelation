# aggregateKS.R
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

# Load the list of experiments
listOfExperiments <- as.matrix(read.table(paste(directoryResults, "/listOfExperiments", sep = ""), sep = " ", header = FALSE))
# Load the parameter list
paramList <- as.matrix(read.table(paste(directoryResults, "/paramList", sep = ""), sep = " ", header = FALSE))
# define some colors
color <- c("black", "blue", "green", "red", "darkblue", "darkgoldenrod4", "darkgrey", "cadetblue1", "darkorchid")

# for all parameter pairs (cinter and cintra) aggregate the ks data
for (cinterIndex in c(1:dim(paramList)[1])) {
  resultDir <- paste(directoryResults, "/", gsub("\\.", "-", sprintf("%.2f", paramList[cinterIndex, 1])), sep = "")
  for (cintraIndex in c(1:dim(paramList)[1])) {
    resultSubDir <- paste(resultDir, "/", gsub("\\.", "-", sprintf("%.2f", paramList[cinterIndex, 1])), "_", gsub("\\.", "-", sprintf("%.2f", paramList[cintraIndex, 2])), sep = "")
    ksData <- matrix(NA, length(listOfExperiments), 3)
    
    # for all experiments
    resultPDF <- list()
    minresultPDF <- -Inf
    maxresultPDF <- Inf
    resultLDFKT <- list()
    minresultLDFKT <- -Inf
    maxresultLDFKT <- Inf
    for (expIndex in c(1:length(listOfExperiments))) {
      ksData[expIndex, 1] <- listOfExperiments[expIndex]
      load(paste(rootDirectory, "/sbm_", listOfExperiments[expIndex], "/sbm_non_digraph_", listOfExperiments[expIndex], "_2_", gsub("\\.", "-", sprintf("%.2f", paramList[cinterIndex, 1])), "/sbm_non_digraph_", listOfExperiments[expIndex], "_2_", gsub("\\.", "-", sprintf("%.2f", paramList[cinterIndex, 1])), "_", gsub("\\.", "-", sprintf("%.2f", paramList[cintraIndex, 2])), "/ana-ks.RData", sep = ""))
      ksData[expIndex, 2] <- ksResult$p.value
      ksData[expIndex, 3] <- log(ksResult$p.value)
      # save the PDF of the current configuration in the resulting directory
      resultPDF[[length(resultPDF) + 1]] <- read.table(paste(rootDirectory, "/sbm_", listOfExperiments[expIndex], "/sbm_non_digraph_", listOfExperiments[expIndex], "_2_", gsub("\\.", "-", sprintf("%.2f", paramList[cinterIndex, 1])), "/sbm_non_digraph_", listOfExperiments[expIndex], "_2_", gsub("\\.", "-", sprintf("%.2f", paramList[cinterIndex, 1])), "_", gsub("\\.", "-", sprintf("%.2f", paramList[cintraIndex, 2])), "/hist_", listOfExperiments[expIndex], "_result.dat", sep = ""), header = FALSE)
      write.table(resultPDF[[length(resultPDF)]], paste(resultSubDir, "/hist_", listOfExperiments[expIndex], "_result.dat", sep = ""), col.names = FALSE, row.names = FALSE)
      if ((min(resultPDF[[length(resultPDF)]][is.finite(resultPDF[[length(resultPDF)]][, 2]), 2]) >= minresultPDF)) {
        minresultPDF <- min(resultPDF[[length(resultPDF)]][is.finite(resultPDF[[length(resultPDF)]][, 2]), 2])
      }
      if ((max(resultPDF[[length(resultPDF)]][is.finite(resultPDF[[length(resultPDF)]][, 2]), 2]) <= maxresultPDF)) {
        maxresultPDF <- max(resultPDF[[length(resultPDF)]][is.finite(resultPDF[[length(resultPDF)]][, 2]), 2])
      }
      # save the large-deviation rate function of the current configuration in the resulting directory
      resultLDFKT[[length(resultLDFKT) + 1]] <- read.table(paste(rootDirectory, "/sbm_", listOfExperiments[expIndex], "/sbm_non_digraph_", listOfExperiments[expIndex], "_2_", gsub("\\.", "-", sprintf("%.2f", paramList[cinterIndex, 1])), "/sbm_non_digraph_", listOfExperiments[expIndex], "_2_", gsub("\\.", "-", sprintf("%.2f", paramList[cinterIndex, 1])), "_", gsub("\\.", "-", sprintf("%.2f", paramList[cintraIndex, 2])), "/hist_", listOfExperiments[expIndex], "_rate-function.dat", sep = ""), header = FALSE)
      write.table(resultLDFKT[[length(resultLDFKT)]], paste(resultSubDir, "/hist_", listOfExperiments[expIndex], "_rate-function.dat", sep = ""), col.names = FALSE, row.names = FALSE)
      if ((min(resultLDFKT[[length(resultLDFKT)]][is.finite(resultLDFKT[[length(resultLDFKT)]][, 2]), 2]) >= minresultLDFKT)) {
        minresultLDFKT <- min(resultLDFKT[[length(resultLDFKT)]][is.finite(resultLDFKT[[length(resultLDFKT)]][, 2]), 2])
      }
      if ((max(resultLDFKT[[length(resultLDFKT)]][is.finite(resultLDFKT[[length(resultLDFKT)]][, 2]), 2]) <= maxresultLDFKT)) {
        maxresultLDFKT <- max(resultLDFKT[[length(resultLDFKT)]][is.finite(resultLDFKT[[length(resultLDFKT)]][, 2]), 2])
      }
    }
    # plot the rate function
    pdf(paste(resultSubDir, "/rate-function.pdf", sep = ""))
    par(oma=c(0, 0, 0, 3))
    plot(resultLDFKT[[1]][, 1], resultLDFKT[[1]][, 2], ylim = c(minresultLDFKT, maxresultLDFKT), xlab = "s", ylab = "PHI(s)")
    for (i in c(2:length(resultLDFKT))) {
      points(resultLDFKT[[i]][, 1], resultLDFKT[[i]][, 2], col = color[i], pch = i)
    }
    legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,
           c(paste(listOfExperiments)),
           col = color[c(1:length(resultLDFKT))],
           pch = c(1:length(resultLDFKT)))
    dev.off()
    
    write.table(ksData, paste(resultSubDir, "/ks-p-values.dat", sep = ""), col.names = FALSE, row.names = FALSE)
    
    # do the fitting business
    aI <- 0
    bI <- -1
    cI <- 0
    
    fitPVal <- NA
    fitPValChiSqrdf <- NA
    tryCatch(
      {
        fitPVal <- nlsLM(ksData[, 2] ~ expGrowth(ksData[, 1], a, b, c),
                         start = list(a = aI, b = bI, c = cI),
                         control = nls.lm.control(maxiter = 1024))
        fitPValChiSqrdf <- chiSqrdf(fitPVal)
      }, error = function(e) print("fit failed")
    )
    fitlogPVal <- NA
    fitlogPValChiSqrdf <- NA
    tryCatch(
      {
        fitlogPVal <- nlsLM(ksData[, 3] ~ expGrowth(ksData[, 1], a, b, c),
                            start = list(a = aI, b = bI, c = cI),
                            control = nls.lm.control(maxiter = 1024))
        fitlogPValChiSqrdf <- chiSqrdf(fitlogPVal)
      }, error = function(e) print("fit failed")
    )
    
    # save the fit results as text file and as image
    x <- c(0:1000) / 1000 * (max(ksData[, 1]) - min(ksData[, 1])) + min(ksData[, 1])
    coefPVal <- NA
    yPVal <- NA
    if (!is.null(dim(fitPVal))) {
      coefPVal <- summary(fitPVal)$coefficients
      yPVal <- expGrowth(x, coefPVal[1], coefPVal[2], coefPVal[3])
      write.table(cbind(x, yPVal), paste(resultSubDir, "/ks-p-fit.dat", sep = ""), col.names = FALSE, row.names = FALSE)
    }
    coeflogPVal <- NA
    ylogPVal <- NA
    if (!is.null(dim(fitlogPVal))) {
      coeflogPVal <- summary(fitlogPVal)$coefficients
      ylogPVal <- expGrowth(x, coeflogPVal[1], coeflogPVal[2], coeflogPVal[3])
      write.table(cbind(x, ylogPVal), paste(resultSubDir, "/ks-log-p-fit.dat", sep = ""), col.names = FALSE, row.names = FALSE)
    }
    
    if (is.null(dim(fitPVal)) && is.null(dim(fitlogPVal))) {
      save(file = paste(resultSubDir, "/ks-p-fit.RData", sep = ""), ksData)
    } else if (!is.null(dim(fitPVal)) && is.null(dim(fitlogPVal))) {
      save(file = paste(resultSubDir, "/ks-p-fit.RData", sep = ""), ksData, fitPVal, fitPValChiSqrdf, x, yPVal)
    } else if (is.null(dim(fitPVal)) && !is.null(dim(fitlogPVal))) {
      save(file = paste(resultSubDir, "/ks-p-fit.RData", sep = ""), ksData, fitlogPVal, fitlogPValChiSqrdf, x, ylogPVal)
    } else {
      save(file = paste(resultSubDir, "/ks-p-fit.RData", sep = ""), ksData, fitPVal, fitPValChiSqrdf, fitlogPVal, fitlogPValChiSqrdf, x, yPVal, ylogPVal)
    }
    
    # plot the p-value results
    pdf(paste(resultSubDir, "/ks-p-fit.pdf", sep = ""))
    if (!is.na(yPVal)) {
      plot(ksData[, 1], ksData[, 2], xlab = "N", ylab = "p-value", ylim = c(min(ksData[, 2], yPVal), max(ksData[, 2], yPVal)))
      lines(x, yPVal, col = "red")
    } else {
      plot(ksData[, 1], ksData[, 2], xlab = "N", ylab = "p-value")
    }
    dev.off()
    # plot the log p-value results
    pdf(paste(resultSubDir, "/ks-log-p-fit.pdf", sep = ""))
    if (!is.na(ylogPVal)) {
      plot(ksData[, 1], ksData[, 3], xlab = "N", ylab = "log(p)-value", ylim = c(min(ksData[, 3], ylogPVal), max(ksData[, 3], ylogPVal)))
      lines(x, ylogPVal, col = "red")
    } else {
      plot(ksData[, 1], ksData[, 3], xlab = "N", ylab = "log(p)-value")
    }
    dev.off()
  }
  
  # plot the log(p) values of the ks-test
  resultslogKS <- list()
  minVal <- -Inf
  for (expIndex in c(1:length(listOfExperiments))) {
    currentKS <- as.matrix(read.table(paste(rootDirectory, "/sbm_", listOfExperiments[expIndex], "/sbm_non_digraph_", listOfExperiments[expIndex], "_2_", gsub("\\.", "-", sprintf("%.2f", paramList[cinterIndex, 1])), "/ks-statistics.dat", sep = ""), sep = " ", header = FALSE))
    write.table(currentKS, paste(directoryResults, "/", gsub("\\.", "-", sprintf("%.2f", paramList[cinterIndex, 1])), "/ks-statistics_", listOfExperiments[expIndex], ".dat", sep = ""), col.names = FALSE, row.names = FALSE)
    resultslogKS[[length(resultslogKS) + 1]] <- currentKS
    if ((min(currentKS[is.finite(currentKS[, 5]), 5]) >= minVal)) {
      minVal <- min(currentKS[is.finite(currentKS[, 5]), 5])
    }
  }
  
  # plot the data
  pdf(paste(directoryResults, "/", gsub("\\.", "-", sprintf("%.2f", paramList[cinterIndex, 1])), "/ks-statistics_", gsub("\\.", "-", sprintf("%.2f", paramList[cinterIndex, 1])), ".pdf", sep = ""))
  plot(resultslogKS[[1]][, 2], resultslogKS[[1]][, 5], xlab = "cintra", ylab = "log(p)-value", ylim = c(minVal, 0), col = color[1], pch = 1)
  for (i in c(2:length(resultslogKS))) {
    points(resultslogKS[[i]][, 2], resultslogKS[[i]][, 5], col = color[i], pch = i)
  }
  legend("bottomleft",
         c(paste(listOfExperiments)),
         col = color[c(1:length(resultslogKS))],
         pch = c(1:length(resultslogKS)))
  dev.off()
}

# collect potentially indicators
ksResults <- list()
minKSVal <- -Inf
maxKSVal <- Inf
for (expIndex in c(1:length(listOfExperiments))) {
  ksResult <- as.matrix(read.table(paste(rootDirectory, "/sbm_", listOfExperiments[expIndex], "/ks-result_", listOfExperiments[expIndex], ".dat", sep = ""), sep = " ", header = FALSE))
  write.table(ksResult, paste(directoryResults, "/ks-result_", listOfExperiments[expIndex], ".dat", sep = ""), col.names = FALSE, row.names = FALSE)
  ksResults[[length(ksResults) + 1]] <- ksResult
  if (min(ksResult[is.finite(ksResult[, 2]), 2]) >= minKSVal) {
    minKSVal <- min(ksResult[is.finite(ksResult[, 2]), 2])
  }
  if (max(ksResult[is.finite(ksResult[, 2]), 2]) <= maxKSVal) {
    maxKSVal <- max(ksResult[is.finite(ksResult[, 2]), 2])
  }
  if(all(!is.finite(ksResult[, 2]))) {
    print("all ks non finite")
  }
  
  kMeansResults <- as.matrix(read.table(paste(rootDirectory, "/sbm_", listOfExperiments[expIndex], "/", "kMeans.dat", sep = ""), sep = " ", header = FALSE))
  write.table(kMeansResults, paste(directoryResults, "/kMeans_", listOfExperiments[expIndex], ".dat", sep = ""), col.names = FALSE, row.names = FALSE)
  
  pdf(paste(directoryResults, "/kMeans_", listOfExperiments[expIndex], ".pdf", sep = ""))
  matplot(kMeansResults[, 1], kMeansResults[, 2:8], xlab = "cinter", ylab = "cintra", col = c("black", "blue", "green", "red", "darkblue", "darkgoldenrod4", "darkgrey"), pch = c(0:6), lwd = c(1, 1, 1, 1, 1, 1, 1), lty = c(1:7))
  legend("topright",
         c("center", "sum", "sum < 0.99", "simple 1", "kmeans 1", "simple 2", "kmeans 2"),
         col = c("black", "blue", "green", "red", "darkblue", "darkgoldenrod4", "darkgrey"),
         pch = c(0:6),
         lwd = c(1, 1, 1, 1, 1, 1, 1),
         lty = c(1:7))
  dev.off()
}

pdf(paste(directoryResults, "/ks-result.pdf", sep = ""))
plot(ksResults[[1]][, 1], ksResults[[1]][, 2], ylim = c(minKSVal, maxKSVal), xlab = "cinter", ylab = "cintra")
for (i in c(2:length(ksResults))) {
  points(ksResults[[i]][, 1], ksResults[[i]][, 2], col = color[i], pch = i)
}
legend("bottomleft",
       c(paste(listOfExperiments)),
       col = color[c(1:length(ksResults))],
       pch = c(1:length(ksResults)))
dev.off()

# process the minimal KS results to fit an exponential through the data
for (cinterIndex in c(1:dim(paramList)[1])) {
  cinterminLogPKSResult <- matrix(0, length(listOfExperiments), 3)
  for (expIndex in c(1:length(listOfExperiments))) {
    cinterminLogPKSResult[expIndex, 1] <- listOfExperiments[expIndex]
    current <- ksResults[[expIndex]]
    cinterminLogPKSResult[expIndex, 2] <- current[current[, 1] == paramList[cinterIndex], 2]
    cinterminLogPKSResult[expIndex, 3] <- current[current[, 1] == paramList[cinterIndex], 3]
  }
  # save this view
  write.table(cinterminLogPKSResult, paste(directoryResults, "/", gsub("\\.", "-", sprintf("%.2f", paramList[cinterIndex, 1])), "/minLogP_", gsub("\\.", "-", sprintf("%.2f", paramList[cinterIndex, 1])), ".dat", sep = ""), col.names = FALSE, row.names = FALSE)
  
  # do the fitting business
  aI <- 0
  bI <- 0
  cI <- 0
  fitResults <- nlsLM(cinterminLogPKSResult[, 2] ~ expGrowth(cinterminLogPKSResult[, 1], a, b, c),
                      start = list(a = aI, b = bI, c = cI),
                      control = nls.lm.control(maxiter = 1024))
  
  # save the fit results as text files
  x <- c(0:1000)/1000 * (max(cinterminLogPKSResult[, 1]) - min(cinterminLogPKSResult[, 1])) + min(cinterminLogPKSResult[, 1])
  coef <- summary(fitResults)$coefficients
  y <- expGrowth(x, coef[1], coef[2], coef[3])
  
  write.table(cbind(x, y), paste(directoryResults, "/", gsub("\\.", "-", sprintf("%.2f", paramList[cinterIndex, 1])), "/minLogP-fit_", gsub("\\.", "-", sprintf("%.2f", paramList[cinterIndex, 1])), ".dat", sep = ""), col.names = FALSE, row.names = FALSE)
  
  # calculate chiSqrdf
  fitChiSqrdf <- chiSqrdf(fitResults)
  
  save(file = paste(directoryResults, "/", gsub("\\.", "-", sprintf("%.2f", paramList[cinterIndex, 1])), "/minLogP-fit_", gsub("\\.", "-", sprintf("%.2f", paramList[cinterIndex, 1])), ".RData", sep = ""), fitResults, fitChiSqrdf, x, coef, y)
  
  # plot the results
  pdf(paste(directoryResults, "/", gsub("\\.", "-", sprintf("%.2f", paramList[cinterIndex, 1])), "/minLogP-fit_", gsub("\\.", "-", sprintf("%.2f", paramList[cinterIndex, 1])), ".pdf", sep = ""))
  plot(cinterminLogPKSResult[, c(1, 2)], xlab = "N", ylab = "cintra", ylim = c(min(cinterminLogPKSResult[, 2], y), max(cinterminLogPKSResult[, 2], y)), col = "black", pch = 1, lwd = 1, lty = 1)
  lines(x, y, col = "blue", pch = NA, lwd = 1, lty = 1)
  legend("topright",
         c(paste("cinter = ", paramList[cinterIndex, 1], sep = ""), "a * N^(-b) + c fit"),
         col = c("black", "blue"),
         pch = c(1, NA),
         lwd = c(1, 1),
         lty = c(NA, 1))
  dev.off()
}
