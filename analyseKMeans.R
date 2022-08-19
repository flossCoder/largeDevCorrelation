# analysekMeans.R
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

cinter <- c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)

result <- matrix(NA, length(cinter), 8)
# for all cinter connectivities
for (cinterIndex in c(1:length(cinter))) {
  # load the resulting data
  data <- as.matrix(read.table(paste(directory, "sbm_non_digraph_", numberOfVertices, "_2_", gsub("\\.", "-", sprintf("%.2f", cinter[cinterIndex])), "/kMeans-results.dat", sep = ""),
                               sep = " ", header = FALSE))
  result[cinterIndex, 1] <- cinter[cinterIndex]
  result[cinterIndex, 2] <- min(data[data[, 6] != 1, 2]) # threshold according to the cluster center position fit
  result[cinterIndex, 3] <- min(data[data[, 7] != 1, 2]) # absolute sum of cluster differences != 1
  result[cinterIndex, 4] <- min(data[data[, 7] < 0.99, 2]) # absolute sum of cluster differences < 0.99
  result[cinterIndex, 5] <- min(data[data[, 9] == 1, 2]) # threshold according to the simple heuristics
  result[cinterIndex, 6] <- min(data[data[, 10] == 1, 2]) # threshold according to the k-means heuristics
  result[cinterIndex, 7] <- data[(max(which(data[, 9] == 0, arr.ind = TRUE)) + 1), 2] # threshold according to the simple heuristics
  result[cinterIndex, 8] <- data[(max(which(data[, 10] == 0, arr.ind = TRUE)) + 1), 2] # threshold according to the k-means heuristics
}

write.table(result, file = paste(directory, "/kMeans.dat", sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)
