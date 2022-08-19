# correlation.R
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


# for none scientific file processing => needed for file names
options(scipen = 999)

dataMetropolis <- NA
dataWangLandau <- NA

source(arg)

# open simple sampling
data <- as.matrix(read.table(paste(directory, "ss_", numberOfVertices, "_", numberOfGraphs, ".dat", sep=""),
                             sep = " ", header = FALSE))

# append all results obtained by the Metropolis-Algorithm
if (!is.null(dim(dataMetropolis))) {
  for (i in 1:dim(dataMetropolis)[1]) {
    data <- rbind(data, as.matrix(read.table(paste(directory, "is_", numberOfVertices, "_",
                                                   dataMetropolis[i, 2], "_",
                                                   dataMetropolis[i, 1], ".dat", sep=""),
                                             sep = " ", header = FALSE)))
  }
}

# append all results obtained by the Wang-Landau-Algorithm
if (!is.null(dim(dataWangLandau))) {
  for (i in 1:dim(dataWangLandau)[1]) {
    data <- rbind(data, as.matrix(read.table(paste(directory, "wl_", numberOfVertices, "_",
                                                   dataWangLandau[i, 1], "_",
                                                   dataWangLandau[i, 2], ".dat", sep=""),
                                             sep = " ", header = FALSE)))
  }
}

#' This function performs a chi-square test, saves the result plus a contingency table and returns both.
#' 
#' @param x First variable.
#' @param y Second variable.
#' @param name A sensefull name for saving the results.
#' 
#' @return [[contingency-table, chi-square-result, cor(x, y)]]
calulateChisqTest <- function(x, y, name) {
  # contingency table
  tab <- table(x, y)
  write.table(tab, file = paste(directory, "ana-cor_", name, ".out", sep = ""), col.names = FALSE, row.names = FALSE)

  # calculate chi-square test
  corRes <- chisq.test(tab)
  writeLines(capture.output(print(corRes)), con = paste(directory, "ana-cor_", name, ".cor", sep = ""))
  
  # calculate ch-square test with simulation of the p-values
  corResMC <- chisq.test(tab, simulate.p.value = TRUE, B = 1000)
  writeLines(capture.output(print(corResMC)), con = paste(directory, "ana-cor-mc_", name, ".cor", sep = ""))
  
  return(list(tab, corRes, corResMC, cor(x, y), qchisq(0.95, corRes$parameter)))
}

# Compare largest components vs. number of components:
lc_vs_nc <- calulateChisqTest(data[, 2], data[, 3], "lc-vs-nc")

# Compare largest components vs. number of edges:
lc_vs_ne <- calulateChisqTest(data[, 2], data[, 4], "lc-vs-ne")

# Compare number of components vs number of edges:
nc_vs_ne <- calulateChisqTest(data[, 3], data[, 4], "nc-vs-ne")

# Save all results as R workspace:
save.image(paste(directory, "ana-cor_chi-sqr.RData", sep = ""))
