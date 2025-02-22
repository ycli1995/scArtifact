% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/merge.R
\name{SCE-merge}
\alias{SCE-merge}
\alias{merge}
\alias{merge,SingleCellExperiment,SingleCellExperiment-method}
\alias{merge,SingleCellExperiment,list_OR_List-method}
\alias{merge,ChromExperiment,ChromExperiment-method}
\alias{merge,ChromExperiment,list_OR_List-method}
\alias{merge,SingleCellMultiExperiment,SingleCellMultiExperiment-method}
\alias{merge,SingleCellMultiExperiment,list_OR_List-method}
\title{Merge single-cell data objects}
\usage{
\S4method{merge}{SingleCellExperiment,SingleCellExperiment}(
  x,
  y,
  assays = NULL,
  reducedDims = NULL,
  altExps = NULL,
  altExps.assays = NULL,
  label = NULL,
  add.idx = NULL,
  idx.collapse = "_",
  verbose = TRUE,
  ...
)

\S4method{merge}{SingleCellExperiment,list_OR_List}(
  x,
  y,
  assays = NULL,
  reducedDims = NULL,
  altExps = NULL,
  altExps.assays = NULL,
  label = NULL,
  add.idx = NULL,
  idx.collapse = "_",
  verbose = TRUE,
  ...
)

\S4method{merge}{ChromExperiment,ChromExperiment}(
  x,
  y,
  assays = NULL,
  reducedDims = NULL,
  altExps = NULL,
  altExps.assays = NULL,
  label = NULL,
  add.idx = NULL,
  idx.collapse = "_",
  verbose = TRUE,
  ...
)

\S4method{merge}{ChromExperiment,list_OR_List}(
  x,
  y,
  assays = NULL,
  reducedDims = NULL,
  altExps = NULL,
  altExps.assays = NULL,
  label = NULL,
  add.idx = NULL,
  idx.collapse = "_",
  verbose = TRUE,
  ...
)

\S4method{merge}{SingleCellMultiExperiment,SingleCellMultiExperiment}(
  x,
  y,
  experiments = NULL,
  global.reducedDims = NULL,
  assays = NULL,
  reducedDims = NULL,
  altExps = NULL,
  altExps.assays = NULL,
  label = NULL,
  add.idx = NULL,
  idx.collapse = "_",
  verbose = TRUE,
  ...
)

\S4method{merge}{SingleCellMultiExperiment,list_OR_List}(
  x,
  y,
  experiments = NULL,
  global.reducedDims = NULL,
  assays = NULL,
  reducedDims = NULL,
  altExps = NULL,
  altExps.assays = NULL,
  label = NULL,
  add.idx = NULL,
  idx.collapse = "_",
  verbose = TRUE,
  ...
)
}
\arguments{
\item{x}{A \code{SingleCellExperiment} or \code{SingleCellMultiExperiment}
object.}

\item{y}{A single object or a list of objects. The class of each object
should be the same as \code{x}.}

\item{assays}{Names of \code{\link[SummarizedExperiment]{assays}} to be
merged. If set to \code{NULL}, common assays will be automatically searched
from input objects. If set to \code{FALSE}, no assay will be merged.}

\item{reducedDims}{Names of \code{\link[SingleCellExperiment]{reducedDims}}
to be merged. If set to \code{NULL}, common reductions will be automatically
searched from input objects. If set to \code{FALSE}, no reduction will be
merged.}

\item{altExps}{Names of \code{\link[SingleCellExperiment]{altExps}} to be
merged. If set to \code{NULL}, common alternative experiments will be
automatically searched from input objects. If set to \code{FALSE}, no
\code{altExps} will be merged.}

\item{altExps.assays}{Names of \code{assays} to be merged in \code{altExps}.}

\item{label}{Add a new column to \code{\link[SummarizedExperiment]{colData}}
of each input to hold the batch information.}

\item{add.idx}{A character vector of \code{length(x = c(x, y))}; appends the
corresponding values to the start of each objects' cell names.}

\item{idx.collapse}{String to collapse \code{add.idx} and cell names of each
input object.}

\item{verbose}{Display progress.}

\item{...}{Arguments passed to other methods.}

\item{experiments}{For \code{SingleCellMultiExperiment}, names of
\code{\link[MultiAssayExperiment]{experiments}} to be merged. If set to
\code{NULL}, all experiments will be used.}

\item{global.reducedDims}{For \code{SingleCellMultiExperiment}, names of the
global \code{reducedDims} to be merged.}
}
\value{
A merged object
}
\description{
Methods to combine multiple single-cell data objects. The behaviors of
merging are the same as \code{\link[SeuratObject]{merge}} function from
\pkg{Seurat}.
}
\examples{
# BPCExperiment
bpce1 <- load_example_sce()
bpce2 <- load_example_sce(dataset = "unsorted")
bpce <- merge(bpce1, bpce2, add.idx = c("sorted", "unsorted"))
bpce

# ChromExperiment
cbpce1 <- load_example_csce()
cbpce2 <- load_example_csce(dataset = "unsorted")
cbpce <- merge(cbpce1, cbpce2, add.idx = c("sorted", "unsorted"))
cbpce

# SingleCellMultiExperiment
scme1 <- load_example_scme()
scme2 <- load_example_scme(dataset = "unsorted")
scme <- merge(scme1, scme2)
scme

}
