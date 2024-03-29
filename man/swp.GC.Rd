% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DNAanalysis_additional-fxns.R
\name{swp.GC}
\alias{swp.GC}
\title{Sliding Widow Plot- GC content}
\usage{
swp.GC(windowsize, inputseq)
}
\arguments{
\item{windowsize}{number of nucleotide sequences in a single step}

\item{inputseq}{a single vector of individual character}
}
\description{
This function allows you create a sliding window plot of GC content within a specified sequence; originally created by Avril Coghlan, from Little Book of R
}
\examples{
swp.GC(100, FISHsequence)
GC
}
\keyword{swp}
