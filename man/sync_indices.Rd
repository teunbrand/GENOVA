% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sync_indices.R
\name{sync_indices}
\alias{sync_indices}
\title{Synchronising indices across experiments.}
\usage{
sync_indices(explist)
}
\arguments{
\item{explist}{A \code{list} of GENOVA \code{contacts} objects of the same
resolution.}
}
\value{
A \code{list} of GENOVA \code{contacts} objects
}
\description{
Particularly when Hi-C data is loaded from different sources, it occurs that
the assigned bin indices don't match between datasets. This function
re-indexes a list of experiments to a shared set of bin indices.
}
\details{
Chromosome naming conventions are expected to be equal among input,
  e.g. all \code{contacts} objects use \code{"chr1"}, \code{"chr2"} etc. or
  all \code{contacts} objects use \code{"1"}, \code{"2"} etc.

  The first \code{contacts} object is used as a template. This should only
  matter when the different experiments have slightly different end-positions
  at e.g. chromosome ends, in which case the ends of the first
  \code{contacts} object is used..
}
\examples{
\dontrun{
# Data loaded from Hi-C Pro
exp1 <- load_contacts(signal_path = "exp1_10kb_iced.matrix",
                      indices_path = "exp1_10kb_abs.bed",
                      sample_name = "exp1")
# Data loaded from Juicer
exp2 <- load_contacts("exp2_10kb.cooler",
                      balancing = TRUE,
                      sample_name = "exp2")
synched <- sync_indices(list(exp1, exp2))
}
}
