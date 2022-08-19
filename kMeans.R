# kMeans.R
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
require("cluster")

#' Calculate the cluster according to the simple clustering rule:
#' Assume a new cluster, when the current cluster is seperated in x-direction from the next point.
#' 
#' @param histogram The given histogram.
#' 
#' @return The indexes of the clusters.
calculateClusters <- function(histogram) {
  result <- list()
  cur <- 2
  currentCluster <- matrix(histogram[1, 1])
  while (cur < (dim(histogram)[1] - 1)) {
    if (histogram[cur, 1] == (histogram[cur - 1, 1] + 1)) {
      currentCluster[length(currentCluster) + 1] <- histogram[cur, 1] # add current node to cluster
    } else {
      result[[length(result) + 1]] <- currentCluster
      currentCluster <- matrix(histogram[cur, 1])
    }
    
    cur <- cur + 1
  }
  result[[length(result) + 1]] <- currentCluster
  
  return(result)
}

#' Generate the tree of the given histogram, where on each level the node with the smallest
#' probability is removed. Save in each step the indexes of the different clusters.
#' 
#' @param histogram The given histogram.
#' 
#' @return The tree containing for each level the indexes of the different clusters according to
#'         the simple clustering method.
generateTree <- function(histogram) {
  tree <- list()
  while (dim(histogram)[1] >= 3) {
    # delete the smallest probability from the histogram
    minP <- min(histogram[, 2])
    index <- min(which(histogram[, 2] == minP, arr.ind = TRUE))
    histogram <- histogram[-index,]
    tree[[length(tree) + 1]] <- calculateClusters(histogram)
  }
  
  return(tree)
}

#' Determine the minimum cluster size for the given set of clusters.
#' 
#' @param clusters The given set of clusters.
#' 
#' @return The size of the smallest cluster.
smallestCluster <- function(clusters) {
  minClusterSize <- Inf
  for (cluster in clusters) {
    if (length(cluster) <= minClusterSize) {
      minClusterSize <- length(cluster)
    }
  }
  return(minClusterSize)
}

#' Auxiliary function to aggregate the indexes of a cluster set for further processing.
#' 
#' @param clusters The given set of clusters.
#' 
#' @return A list containing all cluster indexes.
getIndexesOfAllClusters <- function(clusters) {
  index <- c()
  for (cluster in clusters) {
    index <- c(index, cluster)
  }
  return(index)
}

#' Run the k-means clustering algorithm on the given histogram with the given indexes (tree-level).
#' 
#' @param clusters The given set of clusters.
#' @param histogram The given histogram.
#' 
#' @return The k-Means result object.
runKMeans <- function(clusters, histogram) {
  index <- getIndexesOfAllClusters(clusters)
  return(kmeans(histogram[index, 1:2], 2, nstart = 1000))
}

#' Run the k-medoids clustering algorithm on the given histogram with the given indexes (tree-level).
#' 
#' @param clusters The given set of clusters.
#' @param histogram The given histogram.
#' 
#' @return The pam result object.
runKMedoids <- function(clusters, histogram) {
  index <- getIndexesOfAllClusters(clusters)
  return(pam(histogram[index, 1:2], 2))
}

#' Analyse the given cluster-tree.
#' 
#' @param The tree containing for each level the indexes of the different clusters according to
#'        the simple clustering method.
#' @param histogram The given histogram.
#' 
#' @return A matrix containing several information about the current clusters obtained by k-Means /
#'         k-Medoids.
analyseTree <- function(tree, histogram) {
  clusters <- matrix(NA, length(tree) - 3, 13) # first col: # clusters, second col: smallest clustersize
  for (t in c(1:(length(tree) - 3))) {
    clusters[t, 1] <- t
    clusters[t, 2] <- length(tree[[t]]) # number of clusters
    clusters[t, 3] <- smallestCluster(tree[[t]])
    currentKMeans <- runKMeans(tree[[t]], histogram)
    clusters[t, 4] <- currentKMeans$centers[1]
    clusters[t, 5] <- currentKMeans$centers[3]
    clusters[t, 6] <- currentKMeans$centers[2]
    clusters[t, 7] <- currentKMeans$centers[4]
    clusters[t, 8] <- currentKMeans$betweenss
    clusters[t, 9] <- currentKMeans$tot.withinss
    currentKMedoids <- runKMedoids(tree[[t]], histogram)
    clusters[t, 10] <- currentKMedoids$medoids[1]
    clusters[t, 11] <- currentKMedoids$medoids[3]
    clusters[t, 12] <- currentKMedoids$medoids[2]
    clusters[t, 13] <- currentKMedoids$medoids[4]
  }
  
  return(clusters)
}

# for none scientific file processing => needed for file names
options(scipen = 999)

dataMetropolis <- NA
dataWangLandau <- NA

source(arg)

# open the resulting histogram
histogram <- as.matrix(read.table(paste(directory, "hist_", numberOfVertices, "_result.dat", sep = ""),
                                  sep = " ", header = FALSE))

# make sure, the histogram is normalized
histogram[, 2] <- histogram[, 2] / sum(histogram[, 2])

tree <- generateTree(histogram)
clusters <- analyseTree(tree, histogram)
twoClusters <- clusters[, 2] == 2

#fitRes <- lm(clusters[, 1]~abs(clusters[, 4]-clusters[, 6]))

fitRes <- nlsLM(abs(clusters[, 4] - clusters[, 6]) ~ m * clusters[, 1] + b,
                start = list(b = 1, m = -1),
                control = nls.lm.control(maxiter = 1024))

slope <- (abs(clusters[1, 4] - clusters[1, 6]) - abs(clusters[dim(clusters)[1], 4] - clusters[dim(clusters)[1], 6])) / (clusters[1, 1] - clusters[dim(clusters)[1], 1])
intercept <- abs(clusters[1, 4] - clusters[1, 6]) - slope#clusters[1, 1] * (1 + slope)

difference <- abs(slope * clusters[, 1] + intercept - abs(clusters[, 4] - clusters[, 6]))
ind <- which(difference != min(difference), arr.ind = TRUE)
decision <- sum(ind) / length(ind)

pdf(paste(directory, "kMeans.pdf", sep = ""))
plot(histogram[, 1], histogram[, 2])
twoKMeans <- FALSE
if (is.na(decision)) {
  # one peak phase
  y <- max(histogram[, 2])
  x <- histogram[histogram[, 2] == y, 1]
  points(x, y, col = "blue")
} else {
  # two peak phase
  twoKMeans <- TRUE
  points(clusters[ind, 4], clusters[ind, 5], col = "blue")
  points(clusters[ind, 6], clusters[ind, 7], col = "blue")
}
twoSimple <- FALSE
if (length(clusters[twoClusters, 3]) > 0 && max(clusters[twoClusters, 3]) > 1) {
  index <- clusters[max(which(max(clusters[twoClusters, 3]) == clusters[, 3] & (twoClusters), arr.ind = TRUE)), 1]
  c1 <- histogram[tree[[index]][[1]],]
  y1 <- max(c1[, 2])
  x1 <- c1[(c1[, 2] == y1), 1]
  c2 <- histogram[tree[[index]][[2]],]
  y2 <- max(c2[, 2])
  x2 <- c2[(c2[, 2] == y2), 1]
  points(x1, y1, col="red")
  points(x2, y2, col="red")
  twoSimple <- TRUE
}
dev.off()

# determine clusters for histogram and calculate theire mean difference
theClusters <- kmeans(histogram[, 1:2], 2, nstart = 1000)$cluster
absSumClusterDiff <- abs(sum(histogram[theClusters == 1, 2]) - sum(histogram[theClusters == 2, 2]))
absMeanClusterDiff <- abs(mean(histogram[theClusters == 1, 2]) - mean(histogram[theClusters == 2, 2]))

write.table(clusters, file = paste(directory, "kMeans-analysis.dat", sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)

# Save all results as R workspace:
save.image(paste(directory, "kMeans.RData", sep = ""))
