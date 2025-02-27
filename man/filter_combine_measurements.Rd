% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_functions.R
\name{filter_combine_measurements}
\alias{filter_combine_measurements}
\title{Combines columns that have the same name. The columns are joined additively.}
\usage{
filter_combine_measurements(dat, pre_norm = TRUE, post_norm = TRUE)
}
\arguments{
\item{dat}{The readout matrix to filter.}

\item{pre_norm}{Whether to normalize columns before joining them.}

\item{post_norm}{Whether to normalize columns after they are joined.}
}
\value{
A matrix in which columns that had the same name are added and
        (possibly) normalized.
}
\description{
Combines columns that have the same name. The columns are joined additively.
}
