# gatherkMeans.R
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

load(paste(directory2, "/", param, "/kMeans.RData", sep = ""))

input <- strsplit(param, "_")[[1]]

cinter <- as.numeric(gsub("-", ".", input[length(input) - 1]))
cintra <- as.numeric(gsub("-", ".", input[length(input)]))

maxHy <- max(histogram[, 2])
maxHx <- histogram[histogram[, 2] == maxHy, 1]

cat(sprintf("%s %s %s %s %s %s %s %s %s %s\n",
            toString(cinter),
            toString(cintra),
            toString((cintra + cinter) / 2.0), # average connectivity c
            toString(maxHx), # x-value of the maximum of the histogram data
            toString(maxHy), # y-value of the maximum of the histogram data
            toString(cor(predict(fitRes), abs(clusters[, 4] - clusters[, 6]))), # correlation between linear fit of abs x-cluster center differences and the original one
            toString(absSumClusterDiff), # absolute sum of cluster differences
            toString(absMeanClusterDiff), # absolute mean of cluster differences
            toString(twoSimple * 1.0),
            toString(twoKMeans * 1.0)),
    file = paste(directory2, "/kMeans-results.dat", sep = ""),
    append = TRUE)
