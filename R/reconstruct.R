#' Apply a clustering algorithm recursively to a given time course.
#'
#' @export
#' @param readouts The time course for which to find clusters.
#' @param method Either "kmedoids", "kmeans" or any string permitted as a method for stats::hclust.
#' @param sim A similarity matrix used with all methods except "kmeans".
#' @param splitSimilarity Similarity Threshold. If any two elements within a cluster are below this threhsold, another split is initiated.
#' @param combineSimilarity After Splitting, a combination phase is activated. If any two elements between two clusters have a similarity higher than this threshold, the cluster are combined.
#' @param useSilhouette If TRUE, silhouette is used to define number of cluster during splitting, otherwise cluster are always splitted into two new clusters.
#' @param clusterObj If TRUE, a clusterObject with the readouts, similarity and clustering is returned.
#' @return A matrix with two columns: "Clone" and "Barcode" or if clusterObj = TRUE a cluster object, which can be used to plot the clustering.
reconstructRecursive <- function(readouts,
                                 method = "kmedoids",
                                 sim = MultIS::getSimilarityMatrix(readouts = readouts,
                                                                   upper = TRUE),
                                 splitSimilarity = .7,
                                 combineSimilarity = .9,
                                 useSilhouette = TRUE,
                                 clusterObj = FALSE) {


  # delete NA values to have a proper clustering and additionally delete out the 0 only rows in readouts
  sim <- sim[,!apply(is.na(sim), MARGIN = 2, all)]
  sim <- sim[!apply(is.na(sim), MARGIN = 1, all),]
  readouts <- readouts[!apply(readouts == 0, MARGIN = 1, all),]

  # start with one cluster
  clusters <- list(colnames(sim))
  minSimilarities <- min(sim) # minimal similarity in cluster

  # recursively split until all clusters consists of a minimal inner cluster similiarity > splitSimilarity
  while (any(minSimilarities < splitSimilarity)) {

    i <- which(minSimilarities < splitSimilarity)[1]
    cluster <- clusters[[i]]
    clusters[[i]] <- NULL

    smCluster <- sim[cluster, cluster]
    dataCluster <- readouts[cluster,]

    if (length(cluster) == 2) { # direct split
      mapping <- matrix(c('1', '2', cluster), nrow = 2, dimnames = list(NULL, c('Clone', 'Barcode')))
    } else {

      nrCluster <- 2
      if (useSilhouette) {
        bestNrCluster <- MultIS::findBestNrCluster(data = dataCluster,
                                                   sim = smCluster,
                                                   method.reconstruction = method,
                                                   returnAll = TRUE)
        nrCluster <- as.integer(names(which.max(bestNrCluster)))
      }

      mapping <- MultIS::reconstruct(readouts = dataCluster,
                                     sim = smCluster,
                                     method = method,
                                     targetCommunities = nrCluster,
                                     clusterObj = TRUE)$mapping
    }

    clusters <- c(clusters, lapply(unique(mapping[, 'Clone']), function(clone) {
      mapping[mapping[, 'Clone'] == clone, 'Barcode']
    }))

    # read out the min similarities
    minSimilarities <- unlist(lapply(clusters, function(cluster) {
      min(sim[cluster, cluster])
    }))
  }

  mapping <- as.matrix(do.call('rbind', lapply(1:length(clusters), function(i) {
    data.frame(as.character(i), clusters[[i]])
  })))

  colnames(mapping) <- c('Clone', 'Barcode')
  rownames(mapping) <- mapping[, 'Barcode']

  mapping <- mapping[rownames(readouts),]
  rownames(mapping) <- NULL

  # calculate the means for all cluster and put together the ones with a similarity > combineSimilarity
  clusters <- do.call('rbind', lapply(as.character(sort(as.numeric(unique(mapping[, 'Clone'])))), function(cluster) {
    colMeans(readouts[mapping[mapping[, 'Clone'] == cluster, 'Barcode'], , drop = FALSE])
  }))
  rownames(clusters) <- 1:nrow(clusters)
  sm <- getSimilarityMatrix(readouts = clusters) > combineSimilarity

  finalClusters <- list()
  for (i in 1:ncol(sm)) {
    # check whether i or any of elements is already in a cluster
    cI <- unlist(lapply(finalClusters, function(e) {
      any(which(sm[,i]) %in% e)
    }))
    if (is.null(cI) | all(cI == FALSE)) {
      cI <- length(finalClusters) + 1
      finalClusters[[cI]] <- as.integer(which(sm[,i]))
      next
    }
    cI <- which(cI)[1]

    finalClusters[[cI]] <- sort(unique(c(finalClusters[[cI]], as.integer(which(sm[,i])))))
  }
  length_fC <- Inf
  while (length(finalClusters) != length_fC) {
    length_fC <- length(finalClusters)

    for (i in 1:length(finalClusters)) {
      combines <- unlist(lapply(finalClusters, function(e) {
        any(finalClusters[[i]] %in% e)
      }))
      if (sum(combines) < 2) next
      finalClusters[[length(finalClusters) + 1]] <- sort(unique(unlist(finalClusters[which(combines)])))
      finalClusters <- finalClusters[-which(combines)]
      break
    }
  }

  mapping_tmp <- mapping
  for (i in 1:length(finalClusters)) {
    mapping[mapping_tmp[, 'Clone'] %in% as.character(finalClusters[[i]]), 'Clone'] <- as.character(i)
  }

  if (clusterObj) {
    clusterObj <- list(readouts = readouts,
                       sim = sim,
                       mapping = mapping)
    class(clusterObj) <- c(class(clusterObj), "clusterObj")

    return(clusterObj)
  }

  return(mapping)
}


#' Apply a clustering algorithm to a given time course.
#'
#' @export
#' @param readouts The time course for which to find clusters.
#' @param targetCommunities The number of clusters to cluster for.
#' @param method Either "kmedoids", "kmeans" or any string permitted as a method for stats::hclust.
#' @param sim A similarity matrix used with all methods except "kmeans".
#' @param clusterObj If TRUE, a clusterObject with the readouts, similarity and clustering is returned.
#' @return A matrix with two columns: "Clone" and "Barcode" or if clusterObj = TRUE a cluster object, which can be used to plot the clustering.
reconstruct <- function(readouts,
                        targetCommunities,
                        method = "kmedoids",
                        sim = MultIS::getSimilarityMatrix(readouts = readouts,
                                                          upper = TRUE),
                        clusterObj = FALSE) {


  # delete NA values to have a proper clustering and additionally delete out the 0 only rows in readouts
  sim <- sim[,!apply(is.na(sim), MARGIN = 2, all)]
  sim <- sim[!apply(is.na(sim), MARGIN = 1, all),]
  readouts <- readouts[!apply(readouts == 0, MARGIN = 1, all),]

  diag(sim) <- 0

  if (method == "kmedoids") {
    # TODO: inline reconstructKmedoid here
    labM <- reconstructKMedoid(readouts = readouts,
                               targetCommunities = targetCommunities,
                               sim = sim)
  } else if (method == "kmeans") {
    clus <- stats::kmeans(x = readouts, centers = targetCommunities)$cluster
  } else {
    d <- stats::as.dist(max(sim) - sim, diag = TRUE, upper = TRUE)
    hc <- stats::hclust(d = d, method = method)
    clus <- stats::cutree(hc, k = targetCommunities)
  }

  if (method != "kmedoids")
    labM <- matrix(
      c(clus, names(clus)),
      nrow = length(clus),
      ncol = 2,
      dimnames = list(c(), c("Clone", "Barcode")))

  if (clusterObj) {
    clusterObj <- list(readouts = readouts,
                       sim = sim,
                       mapping = labM)
    class(clusterObj) <- c(class(clusterObj), "clusterObj")

    return(clusterObj)
  }

  return(labM)
}

#' Calculate the k-medoids clustering for a given time course.
#'
#' @param readouts The time course for which to find clusters.
#' @param targetCommunities The number of clusters to cluster for.
#' @param sim A similarity matrix for the time course.
#' @return A matrix with two columns: "Clone" and "Barcode".
reconstructKMedoid <- function(readouts,
                               targetCommunities,
                               sim = MultIS::getSimilarityMatrix(readouts = readouts,
                                                                 self = 0,
                                                                 upper = TRUE)) {
  m <- sim
  m <- m / max(m, na.rm = TRUE)
  if (sum(m, na.rm = TRUE) <= 0)
    return(NULL)
  diag(m) <- 1

  m <- stats::as.dist(1 - m)

  clus <- cluster::pam(m, targetCommunities)

  labM <- matrix(
    c(clus$clustering, rownames(readouts)),
    nrow = length(clus$clustering),
    ncol = 2,
    dimnames = list(c(), c("Clone", "Barcode"))
  )

  return(labM)
}
