## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
  # cache = T,
  # autodep = T
)

set.seed(42)

require(MultIS)

## ----eval=FALSE, include=FALSE------------------------------------------------
#  # Example was generated using the following code
#  simData <- simulate(roCompartments = 1,
#                      tps = seq(1, 2*365, 60),
#                      nrClones = 7,
#                      meanNrBarcodes = 6,
#                      cloneAmpSigma = 0.4,
#                      readOutNoise = 0.2,
#                      simulate.clones.params = list(
#                        nrClones = 7, tps = seq(0, 2 * 365, 60),
#                        prolif = list(type = "nDistributedFixed", nDistributed.mean = 0.3,
#                                      nDistributed.sigma = 0, withLF = TRUE, CC = 10000),
#                        diff = list(type = "nDistributedFixed", nDistributed.mean = 0.2,
#                                    nDistributed.sigma = 0, withLF = FALSE),
#                        initdist = list(type = "equal", equal.equc = (0.3/0.2) * 100))
#  )
#  save(simData, file = 'example.RData', version = 2)

## -----------------------------------------------------------------------------
dat <- read.table(file = "example_readouts.csv",
                  sep = ",", header = TRUE, row.names = 1, check.names = FALSE)
dat <- as.matrix(dat)
class(dat) <- c(class(dat), "timeseries")

## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(dat[1:10,], row.names = TRUE, digits = 2)

## ---- fig.width=6, fig.height=4, fig.align="center"---------------------------
plot(dat)

## ----QS-Filtering, fig.width=6, fig.height=4, fig.align="center"--------------
filteredDat <- MultIS::filter_atTP_biggestN(dat, at = "720", n = 10L)
plot(filteredDat)

## ----QS-rSquareSim, warning=F, fig.width=12, fig.height=4, fig.align="center"----
p1 <- MultIS::plotRsquare(dat, "19", "20") +
  ggplot2::annotate("text", 120, 550,
                    label = paste("R^2 ==", round(
                      summary(stats::lm(y ~ 0 + x, data = data.frame(
                        x = dat["19", ],
                        y = dat["20", ])))$r.squared, 3)), parse = T)
p2 <- MultIS::plotRsquare(dat, "20", "21") +
  ggplot2::annotate("text", 500, 250,
                    label = paste("R^2 ==", round(
                      summary(stats::lm(y ~ 0 + x, data = data.frame(
                        x = dat["20", ],
                        y = dat["21", ])))$r.squared, 3)), parse = T)
gridExtra::grid.arrange(p1, p2, ncol = 2)

## ----QS-similarityMatrix, warning=F-------------------------------------------
similarityMatrix <- MultIS::getSimilarityMatrix(dat, parallel = FALSE)

## ----QS-similarityMatrixHeatmap, warning=F, fig.width=7.2, fig.height=6, fig.align="center"----
plot(similarityMatrix)

## ----QS-clusteringC3, warning=F, fig.width=12, fig.height=6, fig.align="center"----
clusterObjC2 <- MultIS::reconstruct(readouts = dat,
                                    targetCommunities = 2,
                                    method = "kmedoids",
                                    clusterObj = T,
                                    sim = similarityMatrix)
clusterObjC4 <- MultIS::reconstruct(readouts = dat,
                                    targetCommunities = 4,
                                    method = "kmedoids",
                                    clusterObj = T,
                                    sim = similarityMatrix)
p1 <- plot(clusterObjC2)
p2 <- plot(clusterObjC4)

gridExtra::grid.arrange(p1, p2, ncol = 2)

## -----------------------------------------------------------------------------
bestNrCluster <- MultIS::findBestNrCluster(
  data = dat,
  sim = similarityMatrix,
  method.reconstruction = "kmedoids",
  method.evaluation = "silhouette",
  returnAll = TRUE)
plotDf <- data.frame(
  k = as.numeric(names(bestNrCluster)),
  score = as.numeric(bestNrCluster)
)
ggplot2::ggplot(plotDf, ggplot2::aes(x = k, y = score, group = 1)) +
  ggplot2::geom_line() +
  ggplot2::geom_point(ggplot2::aes(col = (score == max(score)))) +
  ggplot2::scale_color_manual(values = c("TRUE" = "#FF0000", "FALSE" = "#000000")) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "none")

## ----QS-Silhouette, warning=F, fig.width=7.2, fig.height=6, fig.align="center"----
bestNrCluster <- MultIS::findBestNrCluster(
  data = dat,
  sim = similarityMatrix,
  method.reconstruction = "kmedoids",
  method.evaluation = "silhouette",
  returnAll = FALSE)

## -----------------------------------------------------------------------------
clusterObjBest <- MultIS::reconstruct(
  readouts = dat,
  targetCommunities = bestNrCluster,
  method = "kmedoids",
  clusterObj = TRUE,
  sim = similarityMatrix)
plot(clusterObjBest)

## -----------------------------------------------------------------------------
load("example.RData")

## -----------------------------------------------------------------------------
str(simData, max.level = 1)

## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(simData$barcodeReadouts[1:10,], digits = 2, row.names = TRUE)

## ----QS-Bushman-Clone-Readouts, fig.width=12, fig.height=8, fig.align="center"----
p1 <- plot(simData$cloneCounts) + ggplot2::ggtitle("Basic clonal simulation")
p2 <- plot(simData$cloneReadouts) + ggplot2:: ggtitle("Added clonal differences")
p3 <- plot(simData$barcodeCounts) + ggplot2::ggtitle("Superimposition of integration sites")
p4 <- plot(simData$barcodeReadouts) + ggplot2::ggtitle("Added measurement noise")
gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

## ----QS-Mappings, results="asis", echo=F--------------------------------------
mapping <- data.frame(Clone = unique(simData$mapping[,"Clone"]))
mapping$Barcodes <- sapply(mapping$Clone, function(e) {
                      paste(summary(simData$mapping[simData$mapping[, "Clone"] == e, "Barcode"])[c("Min.", "Max.")], collapse = " - ")
                    })
knitr::kable(mapping)

## ----QS-ARI, warning=F, fig.width=6, fig.height=4, fig.align="center"---------
aris <- sapply(3:12, function(k) {
  clusterObj <- MultIS::reconstruct(simData$barcodeReadouts,
                                    targetCommunities = k,
                                    clusterObj = T,
                                    sim = similarityMatrix)
  mclust::adjustedRandIndex(clusterObj$mapping[,"Clone"], 
                            simData$mapping[,"Clone"])  
})
arisDF <- data.frame(
  k = 3:12,
  ARI = aris,
  stringsAsFactors = F
)
ggplot2::ggplot(arisDF, ggplot2::aes(x = k, y = ARI, colour = col)) +
  ggplot2::geom_line(colour = "black") +
  ggplot2::geom_point(size = 4, ggplot2::aes(color = (ARI == max(ARI)))) +
  ggplot2::scale_color_manual(values = c("TRUE" = "#FF0000", "FALSE" = "#000000")) +
  ggplot2::scale_x_continuous(breaks = 3:12) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "none",
                 text = ggplot2::element_text(size = 16))

