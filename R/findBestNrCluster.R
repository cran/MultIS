#' Finds the best number of clusters according to silhouette
#'
#' @export
#' @param data The barcode data in a matrix.
#' @param sim A similarity matrix.
#' @param method.reconstruction The clustering method to use.
#' @param method.evaluation The evaluation method to use.
#' @param report Whether the current progress should be reported. Note that this will not work if parallel is set to TRUE.
#' @param parallel Whether the clustering should be performed in parallel.
#' @param best The method to use to determine the best clustering.
#' @param returnAll Wheter to return the silhouette score for all clusterings.
#' @param ... passed params to evaluating clustering
#' @return The R^2 value for rows is1 and is2 in matrix dat
findBestNrCluster <- function(data,
                              sim,
                              method.reconstruction = "kmedoids",
                              method.evaluation = "silhouette",
                              report = FALSE,
                              parallel = FALSE,
                              best = max,
                              returnAll = FALSE,
                              ...) {
  if (parallel)
    dox <- foreach::`%dopar%`
  else
    dox <- foreach::`%do%`

  # delete NA values to have a proper clustering and additionally delete out the 0 only rows in readouts
  sim <- sim[,!apply(is.na(sim), MARGIN = 2, all)]
  sim <- sim[!apply(is.na(sim), MARGIN = 1, all),]
  data <- data[!apply(data == 0, MARGIN = 1, all),]

  #diag(sim) <- 0

  k <- 0
  ks <- dox(foreach::foreach(k = 2:(nrow(sim) - 1)), {
    if (report)
      print(sprintf("Trying %d of max %d clusters", k, nrow(sim) - 1))
    currVal <- tryCatch({
      currRec <- MultIS::reconstruct(readouts = data, targetCommunities = k, sim = sim, method = method.reconstruction)
      currVal <- MultIS::evaluateClustering(readouts = data, clustering = currRec, sim = sim, method = method.evaluation, ...)
    }, error = function(x) { NA })

    return(currVal)
  })
  names(ks) <- 2:(nrow(sim) - 1)

  ks <- unlist(ks)
  m <- as.integer(names(ks)[ks == best(ks)])

  if (returnAll)
    return(ks)
  else
    return(m[1])
}
