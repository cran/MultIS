#' Calculate the bw index
#'
#' @param distance Distance or Dis-Similarity Matrix
#' @param clusters The clustering to evaluate.
#' @param bwBalance The balance [0, 1] between inner cluster similarity (Compactness) and the similarity between clusters (Separation).
#'                A balance value < 1 increases the importance of Compactness, whereas a value > 1 increases the importance of Separation.
#' @param indCluster If true, the bw value for all individual clusters is returned
#' @return A score that describes how well the clustering fits the data.
bw <- function(distance, clusters, bwBalance = 1.0, indCluster = FALSE) {

  clusterIDs <- unique(clusters)

  s <- sapply(1:nrow(distance), function(i) {

    iC <- clusters == clusters[i]
    iC[i] <- FALSE
    iC <- which(iC) # innerCluster Elements

    oC <- 1:nrow(distance)
    oC <- oC[!(oC %in% c(i, iC))]

    if(length(iC) == 0) return(0)

    a <- mean(distance[i, iC])

    b <- min(sapply(clusterIDs[!(clusterIDs %in% clusters[i])], function(e) {
      mean(distance[i, clusters == e])
    }))

    return((b - bwBalance * a) / max(bwBalance * a, b))
  })

  ifelse(indCluster,
         return(sapply(clusterIDs, function(e) {mean(s[clusters == e])})),
         return(mean(s)))

}
