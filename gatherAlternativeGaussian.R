# gatherAlternativeGaussian.R
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

load(paste(directory2, "/", param, "/ana-gaussianTest.RData", sep = ""))

input <- strsplit(param, "_")[[1]]

cinter <- as.numeric(gsub("-", ".", input[length(input) - 1]))
cintra <- as.numeric(gsub("-", ".", input[length(input)]))

firstFitPrediction <- predict(fitResult[[1]][[1]])
secondFitPrediction <- predict(fitResult[[3]][[1]])

maxFirst <- histogram[firstFitPrediction == max(firstFitPrediction), 1]
maxSecond <- histogram[secondFitPrediction == max(secondFitPrediction), 1]

maxfFP <- max(firstFitPrediction)
maxsFP <- max(secondFitPrediction)

firstCoefs <- summary(fitResult[[1]][[1]])$coefficients
secondCoefs <- summary(fitResult[[3]][[1]])$coefficients

maxHy <- max(histogram[, 2])
maxHx <- histogram[histogram[, 2] == maxHy, 1]

cat(sprintf("%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n",
            toString(cinter),
            toString(cintra),
            toString((cintra + cinter) / 2.0), # average connectivity c
            toString(maxFirst), # x of max first prediction
            toString(max(firstFitPrediction)), # max first prediction
            toString(maxSecond), # x of max second prediction
            toString(max(secondFitPrediction)), # max second prediction
            toString(maxFirst - maxSecond), # difference of the two maximum x-positions
            toString(maxFirst / maxSecond), # ratio of the two maximum x-positions
            toString(maxfFP - maxsFP), # difference of the two maximum y-positions
            toString(maxfFP / maxsFP), # ratio of the two maximum y-positions
            toString((maxfFP - maxsFP) / (maxfFP + maxsFP)),
            toString(firstCoefs[1]), # a-Fitparameter of the first fit
            toString(firstCoefs[2]), # m-Fitparameter of the first fit
            toString(firstCoefs[3]), # s-Fitparameter of the first fit
            toString(secondCoefs[1]), # a-Fitparameter of the second fit
            toString(secondCoefs[2]), # m-Fitparameter of the second fit
            toString(secondCoefs[3]), # s-Fitparameter of the second fit
            toString(maxHx), # x-value of the maximum of the histogram data
            toString(maxHy)), # y-value of the maximum of the histogram data
    file = paste(directory2, "/gaussianTest.dat", sep = ""),
    append = TRUE)
