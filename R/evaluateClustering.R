#' Evaluate a clustering using the given method
#'
#' @export
#' @param readouts The readouts the clustering and similarity matrix are based on.
#' @param clustering The clustering to evaluate.
#' @param sim The similarity matrix, this clustering is based on.
#' @param method The method to evaluate the given clustering. This might be one of "silhouette", "sdindex", "ptbiserial", "dunn", "bw", or "custom'.
#' @param customEval A custom function to be run for evaluating a clustering. Only used with method "custom".
#' @param ... Further arguments that are passed to a custom function.
#' @return A score that describes how well the clustering fits the data.
evaluateClustering <- function(readouts, clustering, sim, method, customEval = NULL, ...) {
  methods <- c("silhouette", "sdindex", "ptbiserial", "dunn", "bw", "custom")

  if (!is.character(method))
    stop("The method must be given as a string.")
  if (!method %in% methods)
    stop(paste0("The method must be one of: ", paste(methods, collapse = ", ")))

  ret <- switch(method,
               "silhouette" = evaluateClusteringSilhouette(readouts, clustering, sim),
               "sdindex"    = evaluateClusteringSDindex(readouts, clustering, sim),
               "ptbiserial" = evaluateClusteringPtBiserial(readouts, clustering, sim),
               "dunn"       = evaluateClusteringDunn(readouts, clustering, sim),
               "bw"         = evaluateClusteringBw(readouts, clustering, sim, ...),
               "custom"     = evaluateClusteringCustom(readouts, clustering, sim, customEval, ...)

  )

  return(ret)
}

#' Evaluate a clustering using the silhouette index
#'
#' @param readouts The readouts the clustering and similarity matrix are based on.
#' @param clustering The clustering to evaluate.
#' @param sim The similarity matrix, this clustering is based on.
#' @return A score that describes how well the clustering fits the data.
evaluateClusteringSilhouette <- function(readouts, clustering, sim) {
  dsim <- max(sim) - sim
  sil <- cluster::silhouette(as.integer(clustering[,"Clone"]), dsim)
  ssil <- summary(sil)
  return(ssil$avg.width)
}

#' Evaluate a clustering using the SD-index
#'
#' @param readouts The readouts the clustering and similarity matrix are based on.
#' @param clustering The clustering to evaluate.
#' @param sim The similarity matrix, this clustering is based on.
#' @return A score that describes how well the clustering fits the data.
evaluateClusteringSDindex <- function(readouts, clustering, sim) {
  # Code from BCA package
  SD.clv <- function(x, clus, alpha) {
    if(!is.data.frame(x)) x <- as.data.frame(x)
    scatt <- clv::clv.Scatt(x, clus)
    dis <- clv::clv.Dis(scatt$cluster.center)
    SD <- clv::clv.SD(scatt$Scatt, dis, alfa=alpha)
    return(SD)
  }
  cl <- as.integer(clustering[,"Clone"])
  score <- SD.clv(x = readouts, clus = cl, alpha = clv::clv.Dis(readouts))

  return(score)
}

#' Evaluate a clustering using the point-biserial index
#'
#' @param readouts The readouts the clustering and similarity matrix are based on.
#' @param clustering The clustering to evaluate.
#' @param sim The similarity matrix, this clustering is based on.
#' @return A score that describes how well the clustering fits the data.
evaluateClusteringPtBiserial <- function(readouts, clustering, sim) {
  clones <- as.integer(clustering[,"Clone"])
  dsim <- max(sim, na.rm = TRUE) - sim

  combs <- t(utils::combn(x = 1:nrow(readouts), m = 2))
  x <- apply(combs, MARGIN = 1, function(x) { dsim[x[1], x[2]] })
  #y <- apply(combs, MARGIN = 1, function(x, y) { as.integer(clones[x] != clones[y]) })
  y <- apply(combs, MARGIN = 1, function(x) { as.integer(clones[x[1]] != clones[x[2]]) })

  score <- ltm::biserial.cor(x = x, y = y, level = 2)

  return(score)
}

#' Evaluate a clustering using the dunn index
#'
#' @param readouts The readouts the clustering and similarity matrix are based on.
#' @param clustering The clustering to evaluate.
#' @param sim The similarity matrix, this clustering is based on.
#' @return A score that describes how well the clustering fits the data.
evaluateClusteringDunn <- function(readouts, clustering, sim) {
  clones <- as.integer(clustering[,"Clone"])
  dsim <- max(sim, na.rm = TRUE) - sim

  score <- clValid::dunn(distance = dsim, clusters = clones)

  return(score)
}

#' Evaluate a clustering using the bw index
#'
#' @param readouts The readouts the clustering and similarity matrix are based on.
#' @param clustering The clustering to evaluate.
#' @param sim The similarity matrix, this clustering is based on.
#' @param ... Further arguments that are passed to the bw function.
#' @return A score that describes how well the clustering fits the data.
evaluateClusteringBw <- function(readouts, clustering, sim, ...) {
  clones <- as.integer(clustering[,"Clone"])
  dsim <- max(sim, na.rm = TRUE) - sim

  score <- bw(distance = dsim, clusters = clones, ...)

  return(score)
}

#' Evaluate a clustering using a custom evaluation function
#'
#' @param readouts The readouts the clustering and similarity matrix are based on.
#' @param clustering The clustering to evaluate.
#' @param sim The similarity matrix, this clustering is based on.
#' @param customEval The custom function to be run for evaluating a clustering.
#' @param ... Further arguments that are passed to the custom function.
#' @return A score that describes how well the clustering fits the data.
evaluateClusteringCustom <- function(readouts, clustering, sim, customEval, ...) {
  return(customEval(readouts, clustering, sim, ...))
}
