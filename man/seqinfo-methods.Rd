% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ChromExperiment.R
\name{seqinfo-methods}
\alias{seqinfo-methods}
\alias{seqinfo,ChromExperiment-method}
\alias{seqinfo<-,ChromExperiment-method}
\alias{seqlevels,ChromExperiment-method}
\alias{seqlevels<-,ChromExperiment-method}
\alias{seqnames,ChromExperiment-method}
\alias{seqnames<-,ChromExperiment-method}
\alias{seqlengths,ChromExperiment-method}
\alias{seqlengths<-,ChromExperiment-method}
\alias{genome,ChromExperiment-method}
\alias{genome<-,ChromExperiment-method}
\alias{isCircular,ChromExperiment-method}
\alias{isCircular<-,ChromExperiment-method}
\title{Access and modify sequence information for ChromExperiment}
\usage{
\S4method{seqinfo}{ChromExperiment}(x)

\S4method{seqinfo}{ChromExperiment}(x) <- value

\S4method{seqlevels}{ChromExperiment}(x)

\S4method{seqlevels}{ChromExperiment}(x) <- value

\S4method{seqnames}{ChromExperiment}(x)

\S4method{seqnames}{ChromExperiment}(x) <- value

\S4method{seqlengths}{ChromExperiment}(x)

\S4method{seqlengths}{ChromExperiment}(x) <- value

\S4method{genome}{ChromExperiment}(x)

\S4method{genome}{ChromExperiment}(x) <- value

\S4method{isCircular}{ChromExperiment}(x)

\S4method{isCircular}{ChromExperiment}(x) <- value
}
\arguments{
\item{x}{A \code{\link{ChromExperiment}} object.}

\item{value}{A \code{\link[GenomeInfoDb:Seqinfo-class]{Seqinfo}} object or name of a UCSC genome to
store in the \code{\link{ChromExperiment}}}
}
\description{
Methods for accessing and modifying the \code{\link[GenomeInfoDb:Seqinfo-class]{Seqinfo}} object
stored in a \code{\link{ChromExperiment}} object.
}
\details{
These methods are intended to be consistent with methods for
\code{\link[Signac:ChromatinAssay-class]{ChromatinAssay}} in the \pkg{Signac} package.
}
\seealso{
\itemize{
\item \code{\link[GenomeInfoDb]{seqinfo}} in the \pkg{GenomeInfoDb} package.
\item \code{\link[Signac]{seqinfo-methods}}
}
}
