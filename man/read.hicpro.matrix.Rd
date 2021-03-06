% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_contacts.R
\name{read.hicpro.matrix}
\alias{read.hicpro.matrix}
\title{Read hicpro three column matrix format.}
\usage{
read.hicpro.matrix(file, scale_bp = 1e+09)
}
\arguments{
\item{file}{Full path to file.}

\item{scale_bp}{Normalising factor. Set to NULL to skip scale_bp.}
}
\value{
A data.table with normalised counts.
}
\description{
This function loads a HiC-pro file as a matrix. It assumes a three-column
layout: bin1, bin2 and score. All scors are normalised to contacts per *scale_bp*
total contacts.
}
