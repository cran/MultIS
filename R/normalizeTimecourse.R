normalizeTimecourse <- function(barcodeReadouts, rec, recFirst = FALSE, reduceClones = TRUE) {
  norm.bcr <- barcodeReadouts
  clones <- unique(rec[,"Clone"])
  barcodes <- rec[,"Barcode"]

  if (reduceClones) {
    norm.clones <- matrix(NA, nrow = length(clones), ncol = ncol(norm.bcr), dimnames = list(clones, colnames(norm.bcr)))

    for (c in clones) {
      rel.is <- rec[rec[,"Clone"] == c, "Barcode"]
      rel.clone <- norm.bcr[rel.is,,drop = FALSE]
      sum.clone <- colSums(rel.clone, na.rm = TRUE)
      mean.clone <- sum.clone/nrow(rel.clone)
      norm.clones[c,] <- mean.clone
    }
  } else {
    norm.clones <- matrix(NA, nrow = length(barcodes), ncol = ncol(norm.bcr), dimnames = list(barcodes, colnames(norm.bcr)))

    for (c in clones) {
      rel.is <- rec[rec[,"Clone"] == c, "Barcode"]
      rel.clone <- norm.bcr[rel.is,,drop = FALSE]
      # total clone should have the mean size
      mean.clone <- colMeans(rel.clone, na.rm = TRUE)
      norm.clones[rel.is,] <- t((t(rel.clone)/colSums(rel.clone, na.rm = TRUE)) * mean.clone)
    }
  }

  mean.is.per.clone <- nrow(rec)/length(clones)

  rest <- rownames(barcodeReadouts)[!rownames(barcodeReadouts) %in% rec[,"Barcode"]]
  norm.bcr <- norm.bcr[rest,,drop = FALSE]
  norm.bcr <- norm.bcr/mean.is.per.clone
  norm.bcr <- rbind(norm.clones, norm.bcr)

  if (!recFirst && !reduceClones)
    norm.bcr <- norm.bcr[rownames(barcodeReadouts),]

  return(norm.bcr)
}

