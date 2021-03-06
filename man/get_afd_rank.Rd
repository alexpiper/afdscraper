% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download.R
\name{get_afd_rank}
\alias{get_afd_rank}
\title{Get the taxonomic rank of a query on AFD}
\usage{
get_afd_rank(taxa, quiet = FALSE)
}
\arguments{
\item{taxa}{(Character) A query taxonomic name, or vector of taxonomic names}

\item{quiet}{(Logical) Option to print warnings to console}
}
\value{
A character vector of the same length as \code{taxa} containing the taxonomic rank/s of the query taxa.
}
\description{
Get the taxonomic rank of a query on AFD
}
