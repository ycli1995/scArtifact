% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/SingleCellMultiExperiment.R
\name{defaultExp}
\alias{defaultExp}
\alias{defaultExp<-}
\alias{defaultExp,SingleCellMultiExperiment-method}
\alias{defaultExp<-,SingleCellMultiExperiment,character-method}
\title{Default experiment in SingleCellMultiExperiment}
\usage{
defaultExp(x, ...)

defaultExp(x, ...) <- value

\S4method{defaultExp}{SingleCellMultiExperiment}(x, return.exp = FALSE, ...)

\S4method{defaultExp}{SingleCellMultiExperiment,character}(x, ...) <- value
}
\arguments{
\item{x}{A \code{\link{SingleCellMultiExperiment}}}

\item{...}{Arguments passed to other metheds.}

\item{value}{Update the default experiment strings. Use \code{NA} to manipulate
the \code{\link[SingleCellExperiment]{reducedDims}} and \code{\link[SingleCellExperiment]{colPairs}} stored in
\code{\link{int_SCE}}.}

\item{return.exp}{If \code{TRUE}, will return the default experiment object, else
will only return the strings specifying the default experiment. If
\code{defaultExp} is \code{NA}, this argument will be ignored.}
}
\description{
Get or set the default experiment. This is used to control and switch the
experiment to be manipulated.
}
