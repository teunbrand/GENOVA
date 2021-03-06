% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_anchors.R
\name{anchors_shift}
\alias{anchors_shift}
\title{Shift anchors}
\usage{
anchors_shift(IDX, anchors, rel_pos, shift = 1)
}
\arguments{
\item{IDX}{The indices slot of a GENOVA \code{contacts} object.}

\item{anchors}{A \code{matrix} with two columns containing pre-computed
anchor indices. See the \code{\link[GENOVA]{anchors}} documentation.}

\item{rel_pos}{An \code{integer} vector indicating relative positions in
bins, indicating ranges around anchors to lookup.}

\item{shift}{An \code{integer} of length 1 indicating how many bins the
anchors should be shifted.}
}
\value{
A \code{anchors} object with two columns in \code{matrix} format.
}
\description{
Shifts anchors upstream by a specified amount, unless they become out of
bounds. In that case, shift the anchors downstream.
}
\details{
The resulting matrix contains indices to the Hi-C matrix slot in the
  GENOVA experiment.

  An index is considered out of bounds when that index plus the shift size
  and maximum relative position would belong to a different chromosome.
}
\seealso{
\code{\link[GENOVA]{PESCAn}} for context and
  \code{\link[GENOVA]{anchors}}. \code{\link[GENOVA]{anchors_filter_oob}} for
  general out of bounds filtering of anchors.
}
