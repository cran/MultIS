#' Converts a matrix to relative abundances
#'
#' @export
#' @param dat The matrix to convert
#' @return The matrix with all columns in percent
convert_columnwise_percent <- function(dat) {
  clss <- class(dat)
  cols <- colSums(dat) != 0
  if (length(cols) == 0)
    return(dat)
  dat[, cols] <- t(t(dat[, cols, drop = FALSE])/colSums(dat[, cols, drop = FALSE], na.rm = TRUE))
  class(dat) <- clss
  return(dat)
}

#' Filters for ISs that have a minimum occurence
#'
#' @export
#' @param dat The readout matrix to filter.
#' @param at A filter for the columns. Only (partially) matching columns are considered, though all will be returned.
#' @param min The minumum with which an IS has to occurr. This could be either absolute or relative reads. If `at` matches multiple columns, the rowSum will be used.
#' @return A matrix with only the ISs that occurr with a minimum at the selected measurements.
filter_atTP_min <- function(dat, at = "168", min = 0.02) {
  stopifnot(length(min) == 1)
  if (!is.numeric(min)) {
    warning("min is not numeric, converting now")
    min <- as.numeric(min)
  }
  clss <- class(dat)
  tp <- dat[,grep(paste0("^.*", at, ".*$"), colnames(dat)), drop = FALSE] # given timepoint
  tpGTmin <- rowSums(tp >= min, na.rm = TRUE) > 0              # at least min at the last time point
  dat <- dat[tpGTmin,,drop = FALSE]
  class(dat) <- clss
  return(dat)                                     # the data for the previous filters
}

#' Filters for the n biggest ISs
#'
#' @export
#' @param dat The readout matrix to filter.
#' @param at A filter for the columns. Only (partially) matching columns are considered, though all will be returned.
#' @param n The number of biggest ISs to return. If `at` matches multiple columns, the rowSum will be used to rank the ISs. For ties, more than `n` ISs can be returned.
#' @return A matrix with only the n biggest ISs at the selected measurements.
filter_atTP_biggestN <- function(dat, at = "168", n = 50) {
  stopifnot(length(n) == 1)
  if (!is.integer(n)) {
    warning("n is not an integer, converting now")
    n <- as.integer(n)
  }
  clss <- class(dat)
  tp <- dat[,grep(paste0("^.*", at, ".*$"), colnames(dat)), drop = FALSE] # given timepoint
  if (ncol(tp) > 1) {
    tp <- rowSums(tp)
    dat <- dat[rank(tp) > length(tp) - n,]
  } else
    dat <- dat[rank(tp) > nrow(tp) - n,]

  class(dat) <- clss
  return(dat)
}

#' Filters for a minimum number of time points/measurements
#'
#' @export
#' @param dat The readout matrix to filter.
#' @param min The minimum number of measurements where an IS needs to have a value that is not 0 or NA.
#' @return A matrix with only ISs that have more than `min` columns that are not 0 or NA.
filter_nrTP_min <- function(dat, min = 6) {
  stopifnot(length(min) == 1)
  if (!is.integer(min)) {
    warning("min is not an integer, converting now")
    min <- as.integer(min)
  }
  clss <- class(dat)
  dat <- dat[rowSums(!is.na(dat), na.rm = TRUE) >= min,, drop = FALSE]
  dat <- dat[rowSums(dat != 0, na.rm = TRUE) >= min,, drop = FALSE]
  class(dat) <- clss
  return(dat)
}

#' Filters a vector of names and returns the shortest common prefix.
#'
#' @export
#' @param names The vector of names to filter.
#' @param by A regexp that splits the string. The default filters by special characters. A split by character can be achieved by using "." as the regexp.
#' @return The names shortened to the shortest prefix (in chunks defined by the regexp) where all names are unique.
filter_names <- function(names, by = "[_():]|[^_():]*") {
  names <- regmatches(names, gregexpr(by, names))

  maxLen <- max(unlist(lapply(names, length)))
  for (shrtst in 1:maxLen) {
    # might pad names that are too short with NA
    tmp <- lapply(names, `[`, 1:shrtst)
    if (!anyDuplicated(tmp))
      break
  }

  return(lapply(tmp, function(x) { paste(stats::na.omit(x), collapse = "") }))
}

#' Shortens the rownames of a readout matrix to the shortest distinct prefix
#'
#' @export
#' @param dat The readout matrix for which the names should be filtered.
#' @param by The regexp used to split the names.
#' @return A matrix with the names filtered to the shortest unique prefix.
#' @seealso filter.names
filter_ISnames <- function(dat, by = "[_():]|[^_():]*") {
  rownames(dat) <- filter_names(rownames(dat), by = by)
  return(dat)
}

#' Filters for columns containing a certain substring.
#'
#' @export
#' @param dat The readout matrix to filter.
#' @param match The substring that columns must match.
#' @return A readout matrix that only contains the columns whose names contain the substring.
filter_match <- function(dat, match = "E2P11") {
  clss <- class(dat)
  dat <- dat[,grep(paste0(".*", match, ".*"), colnames(dat))]
  class(dat) <- clss
  return(dat)
}

#' Splits a vector of strings by a given regexp, selects and rearranges the parts and joins them again
#'
#' @export
#' @param dat The readout matrix to filter.
#' @param elems The elements to select. They are rearrange in the order that is given via this argument.
#' @param by The string used for splitting the names of the columns.
#' @return A matrix where the names of the columns are split by the given string, rearranged and again joined by the string.
filter_measurement_names <- function(dat, elems = c(1, 3), by = "_") {
  clss <- class(dat)
  oldCN <- colnames(dat)
  newCN <- unlist(lapply(oldCN, function(x) { paste(strsplit(x, split = by)[[1]][elems], collapse = by) }))
  colnames(dat) <- newCN
  class(dat) <- clss
  return(dat)
}

#' Combines columns that have the same name. The columns are joined additively.
#'
#' @export
#' @param dat The readout matrix to filter.
#' @param pre.norm Whether to normalize columns before joining them.
#' @param post.norm Whether to normalize comumns after they are joined.
#' @return A matrix in which columns that had the same name are added and (possibly) normalized.
filter_combine_measurements <- function(dat, pre.norm = TRUE, post.norm = TRUE) {
  clss <- class(dat)

  nm <- colnames(dat)
  ret <- matrix(data = NA, nrow = nrow(dat), ncol = length(unique(nm)), dimnames = list(rownames(dat), unique(nm)))

  for (ucn in unique(nm)) {
    cols <- nm == ucn
    ss <- dat[,cols, drop = FALSE]
    if (pre.norm)
      ss <- t(t(ss)/colSums(ss, na.rm = TRUE))
    ret[,ucn] <- rowSums(ss, na.rm = TRUE)
  }

  if (post.norm)
    ret <- t(t(ret)/colSums(ret, na.rm = TRUE))

  class(ret) <- clss
  return(ret)
}

#' Removes rows that only contain 0 or NA.
#'
#' @export
#' @param dat The readout matrix to filter.
#' @return A matrix where rows that where only 0 or NA are filtered out.
filter_zero_rows <- function(dat) {
  clss <- class(dat)
  dat <- dat[rowSums(dat, na.rm = TRUE) != 0,,drop = FALSE]
  class(dat) <- clss
  return(dat)
}

#' Removes columns that only contain 0 or NA.
#'
#' @export
#' @param dat The readout matrix to filter.
#' @return A matrix where columns that where only 0 or NA are filtered out.
filter_zero_columns <- function(dat) {
  clss <- class(dat)
  dat <- dat[,colSums(dat, na.rm = TRUE) != 0,drop = FALSE]
  class(dat) <- clss
  return(dat)
}
