# analyse-slope.R
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

#' Calculate the derivative of the given data.
#' 
#' @param x The abscissa of the data.
#' @param y The ordinate of the data.
#' 
#' @return A matrix containing the derivative of the given function.
derivative <- function(x, y) {
  return(diff(y) / diff(x))
  
}

# Load the data, that should be analysed.
data <- as.matrix(read.table(paste(directory, "/sum-statistics.dat", sep = ""), sep = " ", header = FALSE))

# Just to make sure everything is fine, order the data according to the intrablockconnectivity.
data[order(data[, 2]),]

# Calculate the result.
derivData <- cbind(data[-1, 1:3], derivative(data[, 2], data[, 5]), derivative(data[, 2], data[, 9]))

# Save the result in text format.
write.table(derivData, file=
              paste(directory, "/derivative-sum-statistics.dat", sep = ""),
            row.names=FALSE, col.names=FALSE)
