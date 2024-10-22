% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/ObjFile-io.R, R/dumpObj.R
\name{ObjFile-IO}
\alias{ObjFile-IO}
\alias{prepInfo}
\alias{writeObjFile}
\alias{readObjFile}
\alias{prepInfo.default}
\alias{prepInfo.list}
\alias{prepInfo.SimpleList}
\alias{prepInfo.data.frame}
\alias{prepInfo.HDF5Array}
\alias{prepInfo.IterableMatrix}
\alias{prepInfo.IterableFragments}
\alias{prepInfo.GRanges}
\alias{prepInfo.SingleCellMultiExperiment}
\alias{prepInfo.StdAssay}
\alias{prepInfo.Graph}
\alias{prepInfo.DimReduc}
\alias{prepInfo.Seurat}
\title{Read and write the OBJECT.json file}
\usage{
prepInfo(object, ...)

writeObjFile(path, class, extra = list())

readObjFile(path)

\method{prepInfo}{default}(object, ...)

\method{prepInfo}{list}(object, ...)

\method{prepInfo}{SimpleList}(object, ...)

\method{prepInfo}{data.frame}(object, file.name = "DataFrame.csv", idx.name = "_index", ...)

\method{prepInfo}{HDF5Array}(
  object,
  file.name = "matrix.h5",
  name = "matrix",
  with.dimnames = TRUE,
  ...
)

\method{prepInfo}{IterableMatrix}(object, dir.name = "matrix", ...)

\method{prepInfo}{IterableFragments}(object, dir.name = "fragments", ...)

\method{prepInfo}{GRanges}(object, file.name = "GRanges.gtf", ...)

\method{prepInfo}{SingleCellMultiExperiment}(object, ...)

\method{prepInfo}{StdAssay}(object, ...)

\method{prepInfo}{Graph}(object, ...)

\method{prepInfo}{DimReduc}(object, ...)

\method{prepInfo}{Seurat}(object, ...)
}
\arguments{
\item{object}{An object}

\item{...}{Arguments passed to other metheds.}

\item{path}{Path to the directory representing an object.}

\item{class}{String specifying the data class of the object.}

\item{extra}{Named list containing extra metadata to be written.}
}
\description{
The \code{OBJECT.json} file provides the metadata of the object that is
represented by the corresponding directory.
}
