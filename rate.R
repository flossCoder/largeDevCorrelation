# rate.R
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

#' Calculate the m for the rate function.
#' 
#' @param y The given value to calculate.
#' 
#' @return log(1.0 - exp(-y))
m <- function(y) {
  return(log(1.0 - exp(-y)))
}

#' Calculate the sup of s such that s / (1 - k * s) - (1 - exp(-c * s)) == 0
#' 
#' @param k The parameter.
#' @param c The connectivity of the graph.
#' 
#' @return sk
sk <- function(k, c) {
  if (k == 0) {
    return(1)
  }
  
  ds <- 0.00001
  lower <- 0.0
  upper <- 1.0 / k - ds
  
   while ((upper - lower) > ds) { # just to prevent infinity loops
     mid <- (upper + lower) / 2.0
     val <- (function(k, c, s) {
       return(s / (1 - k * s) - (1 - exp(-c * s)))
     })(k, c, mid)
     if (val > 0) {
       upper <- mid
     } else {
       lower <- mid
     }
   }
  return(mid)
}

#' Determine the value of k, such that s_k < s < s_k-1.
#' 
#' @param c The connectivity of the graph.
#' @param s The size of the largest component (s in [0, 1]) => divided by number of vertices
#' 
#' @return The value of k.
calcK <- function(c, s) {
  # s0 is per definition 1
  sKp <- 1
  k <- 0
  while (k < 1000000) { # just to prevent infinity loops
    k <- k + 1
    # calculate sk
    sK <- sk(k, c)
    if ((sK < s) && (s <= sKp)) {
      # case the condition holds => return k
      return(k)
    }
    # sk is sk-1 in the next step
    sKp <- sK
  }
  return(k)
}

#' Calculate the analytical rate function of er graphs.
#' 
#' @param c The connectivity of the graph.
#' @param s The size of the largest component (s in [0, 1]) => divided by number of vertices
#' @param k The k parameter.
#' 
#' @return The analytically determined value of the rate function.
phi <- function(c, s, k) {
  a1 <- -k * s * m(c * s)
  a2 <- k * s * log(s)
  a3 <- (1 - k * s) * log(1 - k * s)
  a4 <- k * c * s
  a5 <- -k * (k + 1) * c * s * s / 2.0
  
  return(a1 + a2 + a3 + a4 + a5)
}

#' Calculate the analytically rate function of ER-graphs.
#' 
#' @param numberOfVertices The number of vertices of the graph.
#' @param c The connectivity of the graph.
#' @param kMax The maximal k parameter.
#' 
#' @return A list containing the k-values and the resulting rate function.
doCalculation <- function(numberOfVertices, c, kMax) {
  x <- matrix(0, ncol = 2, nrow = (kMax + 1))
  res <- matrix(0, ncol = 2, nrow = numberOfVertices - 2)
  
  for (i in 1:(kMax + 1)) {
    x[i, 1] <- (i - 1)
    x[i, 2] <- sk((i - 1), c)
  }
  
  step <- 1
  s <- x[(kMax + 1), 2]
  k <- kMax
  while (s < 1) {
    while (s > x[k, 2]) {
      k <- k - 1
    }
    if (step > dim(res)[1]) {
      res <- rbind(res, c(s, phi(c, s, k)))
    } else {
      res[step, 1] <- s
      res[step, 2] <- phi(c, s, k)
    }
    step <- step + 1
    s <- s + 0.01
  }
  
  return(list(x, res))
}
