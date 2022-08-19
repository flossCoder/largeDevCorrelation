# intersectionPointLD.R
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

#' Calculate the intersection points of the er and sbm large-deviation rate-function.
#' 
#' @param er The large-deviation rate-function of the ER graph.
#' @param sbm The large-deviation rate-function of the SBM graph.
#' 
#' @return A list containing all intersection points of the two large-deviation rate-functions.
calculcateIntersectionPoints <- function(er, sbm) {
  # calculate the difference of the the two large deviation rate functions
  diff <- er[, 2] - sbm[, 2]
  diff2 <- diff
  
  diff <- diff[-1] # remove the first line from diff
  diff2 <- diff2[-length(diff2)] # remove the last line form diff
  
  # determine the sign of the n'th * (n - 1)'th element
  sign <- diff * diff2
  
  # determine a boolean vector, with TRUE: intersection happened and FALSE: no intersection
  bool <- !(is.na(sign < 0) | !(sign < 0))
  
  # return the x-values of the intersection points
  return(er[c(FALSE, bool), 1])
}

# prepare some variables, that will be needed in the following
input <- strsplit(param, "_")[[1]] # dirER contains simulation results for the corresponding ER graphs

cinter <- as.numeric(gsub("-", ".", input[length(input) - 1]))
cintra <- as.numeric(gsub("-", ".", input[length(input)]))

c <- (cinter + cintra) / 2.0 # calculate the corresponding er connectivity

numberOfVertices <- as.numeric(input[length(input) - 3])

# open the large-deviation rate-functions for the analysis
# load the experimental rate function of the stochastic blockmodel
sbm <- as.matrix(read.table(paste(directory, "/", param, "/hist_", numberOfVertices, "_rate-function.dat", sep = ""),
                             sep = " ", header = FALSE))

# load the experimental rate function of the corresponding ER graph
er <- as.matrix(read.table(paste(dirER, "/er_non_digraph_", numberOfVertices, "_", gsub("\\.", "-", c),
                                 "/", "hist_", numberOfVertices, "_rate-function.dat", sep = ""),
                           sep = " ", header = FALSE))

# calculate the intersection
intersection <- calculcateIntersectionPoints(er, sbm)

# calculate the closest intersection point to the minimum of the sbm function
minSBM <- sbm[sbm[, 2] == min(sbm[is.finite(sbm[, 2]), 2]), 1]
minDifferenceSBM <- min(abs(intersection - minSBM))
closestIPSBM <- intersection[abs(intersection - minSBM) == minDifferenceSBM]

# calculate the closest intersection point to the minimum of the er function
minER <- er[er[, 2] == min(er[is.finite(er[, 2]), 2]), 1]
minDifferenceER <- min(abs(intersection - minER))
closestIPER <- intersection[abs(intersection - minER) == minDifferenceER]

cat(sprintf("%s %s %s %s %s\n",
            toString(cinter),
            toString(cintra),
            toString(c), # average connectivity c
            toString(minSBM),#closestIPSBM
            toString(minER)),#closestIPER
    file = paste(directory, "/ld-rf-intersectionPoints.dat", sep = ""),
    append = TRUE)
