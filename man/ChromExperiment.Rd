% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R, R/ChromExperiment.R
\docType{class}
\name{ChromExperiment-class}
\alias{ChromExperiment-class}
\alias{ChromExperiment}
\alias{show,ChromExperiment-method}
\title{The ChromExperiment class}
\usage{
ChromExperiment(
  ...,
  sep = c("-", "-"),
  ranges = NULL,
  genome = NULL,
  annotations = NULL,
  fragments = NULL
)

\S4method{show}{ChromExperiment}(object)
}
\arguments{
\item{...}{Arguments passed to the constructor function of
\code{\link[SingleCellExperiment]{SingleCellExperiment}} to fill the slots of the base class.}

\item{sep}{Separators to use for strings encoding genomic coordinates. The
first element is used to separate the chromosome from the coordinates, the
second element is used to separate the start from end coordinate. Only used
if \code{ranges} is \code{NULL}.}

\item{ranges}{A set of \code{\link[GenomicRanges:GRanges-class]{GRanges}} corresponding to the rows
(peaks) of the input matrix.}

\item{genome}{A \code{\link[GenomeInfoDb:Seqinfo-class]{Seqinfo}} object containing basic information
about the genome used. Alternatively, the name of a UCSC genome
can be provided and the sequence information will be downloaded from UCSC.}

\item{annotations}{A \code{\link[GenomicRanges:GRanges-class]{GRanges}} containing annotations for
the genome used}

\item{fragments}{Fragments data for the input matrix. Can be one of the
following:
\itemize{
\item Tabix-indexed fragment files (*.tsv.gz from 10x Genomics).
\item Directories created by \code{\link[BPCells]{write_fragments_dir}}.
\item HDF5 files created by \code{\link[BPCells]{write_fragments_hdf5}}.
\item One or a list of \code{\link[BPCells:IterableFragments-methods]{IterableFragments}} objects. Note that
if the \code{\link[BPCells:IterableFragments-methods]{IterableFragments}} is created by opening the 10x
\verb{*.tsv.gz} using \code{\link[BPCells]{open_fragments_10x}} directly, the cell
names cannot be accessed or modified by \code{\link[BPCells]{cellNames}}.
}}

\item{object}{A \code{ChromExperiment} object}
}
\value{
A \code{ChromExperiment} object
}
\description{
The \code{ChromExperiment} class is designed to represent single-cell
chromatin accessibility (such as ATAC-seq) data using \pkg{BPCells}. It
inherits from the \code{\link[SingleCellExperiment]{SingleCellExperiment}} class and is used in
the same manner. In addition, it supports storage of barcoded and aligned
fragments information via \code{\link{fragments}}, in which the
\code{\link[BPCells:IterableFragments-methods]{IterableFragments}} objects are used as the underlying data.
}
\details{
In this class, rows should represent genomic coordinates (e.g., peaks) while
columns represent samples generated from single cells.

The extra arguments in the constructor (e.g., \code{sep}, \code{ranges}, \code{genome},
\code{annotations} and \code{\link{fragments}}) represent the main extensions
implemented in the \code{ChromExperiment} class. Readers are referred to
the specific documentation pages for more details.
}
