% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download.R
\name{check_afd_name}
\alias{check_afd_name}
\title{Check if name is valid on the Australian Faunal Directory}
\usage{
check_afd_name(taxa, quiet = FALSE)
}
\arguments{
\item{taxa}{(Character) A query taxonomic name, or vector of taxonomic names}

\item{quiet}{(Logical) Option to print warnings to console}
}
\value{
(character) A character vector of the same length as \code{taxa} containing the valid taxonomic name/s on AFD, or NULL if the name was not found.
}
\description{
Check if name is valid on the Australian Faunal Directory
}
