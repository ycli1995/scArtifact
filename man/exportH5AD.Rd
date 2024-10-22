% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/exportH5AD.R
\name{exportH5AD}
\alias{exportH5AD}
\alias{prepH5AD}
\alias{exportH5MU}
\alias{prepH5MU}
\alias{prepH5AD.SingleCellExperiment}
\alias{prepH5AD.Seurat}
\alias{exportH5AD.SingleCellExperiment}
\alias{exportH5AD.Seurat}
\title{Write an object to H5AD file}
\usage{
exportH5AD(object, ...)

exportH5AD(object, ...)

prepH5AD(object, ...)

exportH5MU(object, ...)

prepH5MU(object, ...)

\method{prepH5AD}{SingleCellExperiment}(
  object,
  X = NULL,
  layers = NULL,
  obsm = NULL,
  obsp = NULL,
  varp = NULL,
  uns = NULL,
  raw = NULL,
  ...
)

\method{prepH5AD}{Seurat}(
  object,
  assay = NULL,
  X = NULL,
  layers = NULL,
  obsm = NULL,
  obsp = NULL,
  uns = NULL,
  raw = NULL
)

\method{exportH5AD}{SingleCellExperiment}(
  object,
  file,
  name = "/",
  overwrite = FALSE,
  gzip_level = 0L,
  verbose = TRUE,
  ...
)

\method{exportH5AD}{Seurat}(
  object,
  file,
  name = "/",
  overwrite = FALSE,
  gzip_level = 0L,
  verbose = TRUE,
  assay = NULL,
  X = NULL,
  layers = NULL,
  obsm = NULL,
  obsp = NULL,
  uns = NULL,
  ...
)
}
\arguments{
\item{object}{An object}

\item{...}{Arguments passed to other metheds. Mainly to \code{prepH5AD} and \code{prepH5MU}}

\item{X}{Which \code{\link[SummarizedExperiment]{assay}} to be used as
\code{adata.X}. Default is the first one.}

\item{layers}{Which \code{\link[SummarizedExperiment]{assays}} to be used in
\code{adata.layers}. Default will be all assays except the \code{X} one.}

\item{obsm}{Which \code{\link[SingleCellExperiment]{reducedDim}} to be used
in \code{adata.obsm}. Default will be all reduction data.}

\item{obsp}{Which \code{\link[SingleCellExperiment]{colPair}} to be used in
\code{adata.obsp}.}

\item{varp}{Which \code{\link[SingleCellExperiment]{rowPair}} to be used in
\code{adata.varp}.}

\item{uns}{Which \code{\link[S4Vectors]{metadata}} to be used in \code{adata.uns}.}

\item{raw}{Which \code{\link[SingleCellExperiment]{altExp}} to be used in
\code{adata.raw}. The \code{nrows()} of \code{raw} must be larger than the \code{mainExp}.}

\item{file}{Path to the H5AD file.}

\item{name}{Name of the HDF5 group where to save \code{object}.}

\item{overwrite}{Whether or not to overwrite the existing HDF5 link.}

\item{gzip_level}{Enable zipping at the level given here.}

\item{verbose}{Print progress.}

\item{assayt}{Which Seurat \code{assay} to use. If \code{NULL}, use
\code{DefaultAssay(object)}.}
}
\description{
Write an object to H5AD file

Write an object to H5AD (H5MU) file
}
