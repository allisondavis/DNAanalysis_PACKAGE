% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DNAanalysis_additional-fxns.R
\name{swp.CG}
\alias{swp.CG}
\title{Sliding Widow Plot-CG content}
\usage{
swp.CG(windowsize, inputseq)
}
\arguments{
\item{windowsize}{number of nucleotide sequences in a single step}

\item{inputseq}{a single vector of individual character}
}
\description{
This function allows you create a sliding window plot of CG content within a specified sequence; originally created by Avril Coghlan, from Little Book of R and modified with our CG function
}
\examples{
swp.CG(100, FISHsequence)
CG
}
\keyword{swp}
