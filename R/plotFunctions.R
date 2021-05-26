# ggplotColours <- function(n = 6, h = c(0, 360) + 15){
#   if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
#   grDevices::hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
# }

ggplotColors <-
ggplotColours <- function(n = 6, h = c(0, 360) + 15, l = c(65, 65)) {
  if (length(n) == 1) {
    if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
    grDevices::hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
  } else {
    if ((diff(h) %% 360) < 1)
      h[2] <- h[2] - 360/length(n)
    base_h <- seq(h[1], h[2], length = length(n))
    unlist(lapply(1:length(n), function(i) {
      h <- base_h[i]
      if (n[i] == 1) {
        ls <- (l[1] + l[2]) / 2
      } else {
        ls <- seq(l[1], l[2], length.out = n[i])
      }

      grDevices::hcl(h = h, c = 100, l = ls)
    }))
  }
}

weightedSpringModel <- function(readouts, mapping, gt,
                                sim = getSimilarityMatrix(readouts, self = NA, upper = FALSE, parallel = FALSE),
                                recPal = NULL, clonePal = NULL, lineColor = "#009900FF", seed = 4711L) {
  net <- igraph::graph_from_adjacency_matrix(adjmatrix = sim,
                                             mode = "undirected",
                                             diag = FALSE,
                                             weighted = TRUE)
  set.seed(seed)
  layout <- igraph::layout_with_fr(net, grid = "nogrid")

  v <- data.frame(
    x = layout[,1],
    y = layout[,2],
    IS = colnames(sim)
  )

  if (!missing(mapping) && !is.null(mapping)) {
    v$Rec <- plyr::mapvalues(v$IS,
                             mapping[,"Barcode"],
                             mapping[,"Clone"],
                             warn_missing = FALSE)
    v$Rec <- as.factor(v$Rec)
    if (is.null(recPal)) {
      recPal <- ggplotColours(n = length(levels(v$Rec)))
      names(recPal) <- levels(v$Rec)
    }
  } else {
    v$Rec <- as.factor(0)
    recPal <- c("0" = "grey50")
  }

  if (!missing(gt) && !is.null(gt)) {
    v$Clone <- plyr::mapvalues(v$IS,
                               gt[,"Barcode"],
                               gt[,"Clone"])
    v$Clone <- as.factor(v$Clone)
    if (is.null(clonePal)) {
      clonePal <- ggplotColours(n = length(levels(v$Clone)))
      names(clonePal) <- levels(v$Clone)
    }
  } else {
    v$Clone <- as.factor(0)
    clonePal <- c("0" = "grey25")
  }

  # TODO: This apparently gives the same order... find a way to make this more robust
  el <- igraph::get.edgelist(graph = net, names = FALSE)
  ew <- igraph::get.edge.attribute(net, "weight")

  e <- data.frame(
    from = el[,1],
    to = el[,2],
    Weight = ew**4,
    Width = 1
  )
  e$Weight <- (e$Weight - min(e$Weight))/(max(e$Weight) - min(e$Weight))

  e$from.x <- v$x[e$from]
  e$from.y <- v$y[e$from]
  e$to.x <- v$x[e$to]
  e$to.y <- v$y[e$to]

  x_min <- min(c(e$from.x, e$to.x))
  x_max <- max(c(e$from.x, e$to.x))
  y_min <- min(c(e$from.y, e$to.y))
  y_max <- max(c(e$from.y, e$to.y))

  # cut <- max(sort(abs(e$Weight), decreasing = TRUE)[2 * nrow(v)], quantile(abs(e$Weight), 0.75))
  # avgW = (abs(e$Weight) - cut)/(max(e$Weight) - cut)
  # avgW[avgW < 0] = 0
  #
  # esize <- 15 * exp(-nrow(v)/90) + 1
  # edge.width <- avgW * (esize - 1) + 1
  # edge.width[edge.width < 1] = 1

  p <- ggplot2::ggplot(v, ggplot2::aes_string(x = "x", y = "y", fill = "Rec", col = "Clone")) +
    ggplot2::geom_segment(data = e,
                          ggplot2::aes_string(x = "from.x", y = "from.y", xend = "to.x",
                                              yend = "to.y", size = "3", alpha = "Weight"),
                          color = lineColor, inherit.aes = FALSE) +
    ggplot2::geom_point(size = 12, shape = 21, stroke = 3) +
    ggplot2::geom_text(ggplot2::aes_string(label = "IS"), color = "#000000FF") +
    ggplot2::scale_fill_manual(values = recPal) +
    ggplot2::scale_color_manual(values = clonePal) +
    ggplot2::scale_size_identity() +
    ggplot2::scale_alpha_identity() +
    ggplot2::theme_void() +
    ggplot2::expand_limits(x = c(x_min - 0.12, x_max + 0.12), y = c(y_min - 0.12, y_max + 0.12))

  return(p)
}

#' @importFrom dplyr %>%
#' @importFrom rlang .data
bushmanPlot <- function(bc_dat, aes = NULL, col = NULL, only = NULL, rec = NULL, time = NULL, facet = NULL) {
  bc_dat[is.na(bc_dat)] <- 0
  # bc_dat <- bc_dat[rowSums(bc_dat) != 0,]
  bc_dat <- t(t(bc_dat)/colSums(bc_dat))

  # colnames(bc_dat) <- seq_len(ncol(bc_dat))

  if (!is.null(rec)) {
    rn <- rownames(bc_dat)
    others <- rn[!(rn %in% unique(rec[,"Barcode"]))]
    others <- colSums(bc_dat[others,])
    bc_dat <- rbind(bc_dat[unique(rec[,"Barcode"]),], others)
  }

  m <- reshape2::melt(bc_dat[rowSums(bc_dat) != 0,])
  colnames(m) <- c("IS", "Measurement", "Contribution")
  m$IS <- factor(m$IS, ordered = FALSE, levels = unique(m$IS))

  if (is.null(time)) {
    m$Time <- factor(m$Measurement, ordered = TRUE)
  } else if (is.function(time)) {
    m$Time <- time(m$Measurement)
    x_breaks <- unique(m$Time)
    x_labels <- as.character(x_breaks)
  } else {
    stop("time needs to be either a function or NULL")
  }

  if (!is.null(facet)) {
    if (is.function(facet)) {
      m$Facet <- facet(m$Measurement)
      # Optimization: in each facet, remove ISs that have 0 reads and thus would not be drawn
      m <- as.data.frame(m %>%
                           dplyr::group_by(.data$IS, .data$Facet) %>%
                           dplyr::filter(sum(.data$Contribution) > 0) %>%
                           dplyr::ungroup())
      #m <- m[m$Contribution != 0,]

      for (f in unique(m$Facet))
        if (length(unique(subset(m, .data$Facet == f)$Time)) < 2)
          m <- m[m$Facet != f,]

      m$Facet <- factor(m$Facet, ordered = FALSE)
    } else {
      stop("facet needs to be either a function or NULL")
    }
  }

  if (is.null(aes))
    p <- ggplot2::ggplot(m)
  else
    p <- ggplot2::ggplot(m, aes)

  p <- p + ggplot2::geom_area(ggplot2::aes_string(group = "IS", fill = "IS", x = "Time", y = "Contribution"))

  if (!is.null(col)) {
    if (length(col) == nrow(bc_dat)) {
      cols <- col
    } else {
      cols <- rep("gray50", nrow(bc_dat))
      names(cols) <- rownames(bc_dat)
      cols[names(col)] <- col
    }

    p <- p + ggplot2::scale_fill_manual(values = cols)
  } else if (!is.null(only)) {
    cols <- rep("gray50", nrow(bc_dat))
    cols_fc <- ggplotColours(n = nrow(bc_dat))
    repl <- rownames(bc_dat) %in% only
    cols[repl] <- cols_fc[repl]
    names(cols) <- rownames(bc_dat)

    p <- p + ggplot2::scale_fill_manual(values = cols)
  } else if (!is.null(rec)) {
    cols <- rep("gray50", nrow(bc_dat))
    names(cols) <- rownames(bc_dat)

    clones <- rec[,"Clone"]
    cols[rec[,"Barcode"]] <- clones

    cols <- plyr::mapvalues(cols,
                            unique(clones),
                            ggplotColours(n = length(unique(clones))),
                            warn_missing = FALSE)

    p <- p + ggplot2::scale_fill_manual(values = cols)
  }

  if (is.integer(m$Time))
    p <- p + ggplot2::scale_x_continuous(expand = c(0, 0), breaks = x_breaks, labels = x_labels) +
      ggplot2::geom_vline(ggplot2::aes_string(xintercept = "Time"))
  else if (is.factor(m$Time))
    p <- p + ggplot2::scale_x_discrete(expand = c(0, 0))
  p <- p + ggplot2::scale_y_continuous(expand = c(0, 0))

  if (!is.null(m$Facet))
    p <- p + ggplot2::facet_wrap(~Facet, ncol = 1, strip.position = "right")

  return(p)
}

.time_fn <- function(x) { as.integer(unlist(lapply(strsplit(as.character(x), "_"), utils::head, 1)))}
.comp_fn <- function(x) {unlist(lapply(strsplit( as.character(x), "_"), function(y) { paste(utils::tail(y, -1), collapse = "_")} ))}

bushmanPlotSplitCompartment <- function(timeFn = .time_fn, compFn = .comp_fn, ...) {
  return(bushmanPlot(time = timeFn, facet = compFn, ...))
}

#' @importFrom rlang .data
linePlotSplitClone <- function(bd, rec, order = NULL, mapping = NULL, sim = NULL, silhouetteValues = !is.null(sim), singletons = TRUE, zero.values = TRUE) {
  if (!is.matrix(bd))
    bd <- as.matrix(bd)

  m <- reshape2::melt(bd)
  colnames(m) <- c("IS", "Time", "Count")
  m$Clone <- plyr::mapvalues(m$IS,
                             rec[,"Barcode"],
                             as.integer(rec[,"Clone"]))

  if (!is.null(order))
    m$IS <- factor(x = m$IS, levels = order)
  else
    m$IS <- factor(x = m$IS)

  m$Time <- as.factor(m$Time)
  m$Clone <- as.factor(m$Clone)

  if (silhouetteValues) {
    sil <- cluster::silhouette(as.numeric(rec[,"Clone"]), max(sim) - sim)
    ssil <- summary(sil)
    from <- names(ssil$clus.avg.widths)
    cloneSize <- table(rec[,"Clone"])
    cloneSize <- cloneSize[as.character(sort(as.integer(names(cloneSize))))]
    cloneSize <- as.vector(cloneSize)
    to <- sprintf("%s: %.3f [%d IS]", names(ssil$clus.avg.widths), ssil$clus.avg.widths, cloneSize)
    m$Clone <- plyr::mapvalues(x = m$Clone,
                               from = from,
                               to = to)
  }

  if (!singletons)
    m <- m[!endsWith(as.character(m$Clone), "[1 IS]"),, drop = FALSE]

  if (!zero.values)
    m[m$Count == 0, "Count"] <- NA

  am <- c(mapping, ggplot2::aes_string(x = "Time", y = "Count", color = "Clone", group = "IS"))
  class(am) <- "uneval"
  p <- ggplot2::ggplot(m, mapping = am) +
    ggplot2::geom_path() +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~.data$Clone)

  return(p)
}

#' Plots time series data, which consists of multiple measurements over time / place of different (cols)
#' clones / integration sites (rows).
#'
#' @export
#' @param x The data to plot.
#' @param ... Further arguments are ignored.
#' @return A ggplot object, which can be used to further indiviualize or to plot directly.
plot.timeseries <- function(x, ...) {
  return(bushmanPlot(x) + ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none",
                   axis.title.x = ggplot2::element_blank(),
                   text = ggplot2::element_text(size = 16),
                   plot.margin = ggplot2::margin(.5, .5, .5, .5, "cm")))
}

#' Plots R^2 of two integration sites
#'
#' @export
#' @param dat The matrix that holds the values
#' @param is1 The name of the first row
#' @param is2 The name of the second row
#' @return A ggplot object, which can be used to further indiviualize or to plot directly.
plotRsquare <- function(dat, is1, is2) {
  return(ggplot2::ggplot(as.data.frame(t(as.data.frame(dat[c(is1, is2),], row.names = c("Count1", "Count2")))),
         ggplot2::aes_string(x = "Count1", y = "Count2")) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ 0 + x, se = FALSE) +
    ggplot2::geom_point() +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 75, hjust = 1),
                   legend.position = "none",
                   axis.title.x = ggplot2::element_blank(),
                   text = ggplot2::element_text(size = 16),
                   plot.margin = ggplot2::margin(.5, .5, .5, .5, "cm")) +
    ggplot2::xlab(bquote('Abundance ('~1^st~' IS)')) +
    ggplot2::ylab(bquote('Abundance ('~2^nd~' IS)')))
}

#' Plots the similarity of integration sites
#'
#' @export
#' @param x The matrix that holds the similarity values
#' @param na.rm whether NA values should be deleted beforehand
#' @param ... Further arguments are ignored.
#' @return A ggplot object, which can be used to further indiviualize or to plot directly.
plot.ISSimilarity <- function(x, na.rm = TRUE, ...) {
  mDat <- reshape2::melt(x)
  colnames(mDat) <- c("IS1", "IS2", "Similarity")

  if (na.rm) mDat <- mDat[!is.na(mDat$Similarity),]

  mDat$IS1 <- as.factor(as.numeric(mDat$IS1))
  mDat$IS2 <- as.factor(as.numeric(mDat$IS2))

  return(ggplot2::ggplot(mDat,
                  ggplot2::aes_string(x = "IS1", y = "IS2", fill = "Similarity")) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colours = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")), space = "Lab")(100)) +
#    ggplot2::scale_x_discrete(expand = c(0, 0)) +
#    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 75, hjust = 1),
                   text = ggplot2::element_text(size = 16),
                   axis.text = ggplot2::element_text(size = 7)))
}

#' Plots the clustering based on a clustering object
#'
#' @export
#' @param x The clustering object.
#' @param ... Further arguments are ignored.
#' @return A ggplot object, which can be used to further indiviualize or to plot directly.
plot.clusterObj <- function(x, ...) {
  return(weightedSpringModel(readouts = x$readouts,
                      mapping = x$mapping,
                      gt = x$mapping,
                      sim = x$sim) +
    ggplot2::guides(col = FALSE) +
    ggplot2::theme(legend.position = "none"))
}

stackedBarPlot <- function(bcDat, comp_fn = NULL, rec = NULL) {
	bcDat[is.na(bcDat)] <- 0

	if (!is.null(rec)) {
		bcDat_cl <- do.call(rbind,
				 lapply(unique(rec[,"Clone"]), function(x) {
						colSums(bcDat[rec[rec[, "Clone"] == x, "Barcode"],,drop=FALSE])
				 }))
		bcDat_ot <- colSums(bcDat[!(rownames(bcDat) %in% rec[,"Barcode"]),])
		bcDat <- rbind(bcDat_cl, bcDat_ot)
		rownames(bcDat) <- c(unique(rec[,"Clone"]), "others")
	}

	bcDat <- t(t(bcDat)/colSums(bcDat))

	mdat <- reshape2::melt(as.matrix(bcDat))
	colnames(mdat) <- c("IS", "Measurement", "Value")

	if (!is.null(comp_fn)) {
		mdat$Facet <- comp_fn(mdat$Measurement)
		ggplot2::ggplot(mdat, ggplot2::aes_string(x = "Measurement", y = "Value", fill = "IS")) +
			ggplot2::geom_bar(position = "stack", stat = "identity") +
			ggplot2::facet_grid(~Facet, scales = "free_x") +
			#ggplot2::scale_fill_manual(values = clone_cols) +
			ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 75, hjust = 1)) +
			ggplot2::ylab("Relative abundance")
	} else {
		ggplot2::ggplot(mdat, ggplot2::aes_string(x = "Measurement", y = "Value", fill = "IS")) +
			ggplot2::geom_bar(position = "stack", stat = "identity") +
			ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 75, hjust = 1)) +
			ggplot2::ylab("Relative abundance")
	}
}

