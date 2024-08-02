#' @importFrom methods as callGeneric callNextMethod getMethod is selectMethod
#' setGeneric setMethod slot slotNames
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S3 Generics ##################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Backend path ################################################################

#' Fetch the file/directory path information for an on-disk object
#'
#' Get the path and the link (for HDF5-backed) name for an on-disk object.
#'
#' @param object A disk-backed object.
#' @param ... `r .dot_param`
#'
#' @return A \code{\link{data.frame}} with three columns:
#' \itemize{
#' \item `path`: Absolute path of the backend file/directory.
#' \item `type`: Type of the backend file/directory (dir, h5, ...).
#' \item `group` : Name of HDF5 group in the backend file.
#' }
#'
#' @details
#' This function can be useful to check the file path of an object to avoid IO
#' conflicts.
#'
#' @rdname getPath
#' @export getPath
getPath <- function(object, ...) {
  UseMethod(generic = "getPath", object = object)
}

## Merge matrices ##############################################################

#' Merge Single-cell data matrices
#'
#' Single-cell data matrices are mostly in cell-by-feature formation. In this
#' function, matrices are merged through \code{\link{cbind}}. The union of rows
#' for input matrices are used.
#'
#' @param x A matrix object
#' @param y One or more matrices of the same class or coercible to the same
#' class as `x`
#' @param ... `r .dot_param`
#'
#' @note
#' The column names along all the input matrices must be unique, or an error
#' will be raised.
#'
#' @seealso The manner of `mergeMatrices()` is actually the same as
#' \pkg{Seurat}'s \code{\link[SeuratObject]{merge}}
#'
#' @return A single matrix of type `class(x)`
#'
#' @rdname mergeMatrices
#' @export mergeMatrices
mergeMatrices <- function(x, y, ...) {
  UseMethod(generic = "mergeMatrices", object = x)
}

## BPCells #####################################################################

#' Open fragments file
#'
#' @param x One of the followings:
#' \itemize{
#' \item A `.tsv.gz` fragments file (must has the corresponding `.tbi` index).
#' \item An HDF5 file created by `r .doc_links("write_fragments_hdf5")`.
#' \item A directory created by `r .doc_links("write_fragments_dir")`.
#' \item A `r .doc_links("Fragments")` object from \pkg{Signac}.
#' }
#' @param group The link name of HDF5 group, where the fragments data are
#' stored. Default is 'fragments'.
#' @param ... `r .dot_param`
#'
#' @return An \code{IterableFragments} object
#'
#' @seealso
#' \itemize{
#' \item `r .doc_links("open_fragments_hdf5")`
#' \item `r .doc_links("open_fragments_dir")`
#' \item `r .doc_links("open_fragments_10x")`
#' }
#'
#' @rdname openFragments
#' @export openFragments
openFragments <- function(x, group = NULL, ...) {
  UseMethod(generic = "openFragments", object = x)
}

## Convert to h5ad #############################################################

#' Write an object to H5AD file
#'
#' @param object An object
#' @param file Path of the H5AD file
#' @param ... `r .dot_param`
#'
#' @rdname exportH5AD
#' @export exportH5AD
exportH5AD <- function(object, file, ...) {
  UseMethod(generic = "exportH5AD", object = object)
}


## Format object ###############################################################

#' Format an object
#'
#' @param object An object
#' @param ... `r .dot_param`
#'
#' @rdname formatObj
#' @export formatObj
formatObj <- function(object, ...) {
  UseMethod(generic = "formatObj", object = object)
}

## Diet object #################################################################

#' Slim down a SingleCellExperiment-like object
#'
#' Keep only certain aspects of the `r .doc_links("SingleCellExperiment")`
#' object.
#'
#' @param object A `r .doc_links("SingleCellExperiment")`-like object
#' @param ... `r .dot_param`
#'
#' @rdname dietObj
#' @export dietObj
dietObj <- function(object, ...) {
  UseMethod(generic = "dietObj", object = object)
}

## IO for each part of an object ###############################################

#' IO for different data structures
#'
#' Read or write assorted R objects into appropriate on-disk representations.
#'
#' @param object An SCE-like object.
#' @param ... `r .dot_param`
#'
#' @name Obj-IO
NULL

#' @rdname Obj-IO
#' @export writeObj
writeObj <- function(object, ...) {
  UseMethod(generic = "writeObj", object = object)
}

#' @rdname Obj-IO
#' @export reloadObj
reloadObj <- function(object, ...) {
  UseMethod(generic = "reloadObj", object = object)
}

#' @rdname Obj-IO
#' @export prepObj
prepObj <- function(object, ...) {
  UseMethod(generic = "prepObj", object = object)
}

#' @rdname Obj-IO
#' @export dumpObj
dumpObj <- function(object, ...) {
  UseMethod(generic = "dumpObj", object = object)
}

#' Read and write the OBJECT.json file
#'
#' The `OBJECT.json` file provides the metadata of the object that is
#' represented by the corresponding directory.
#'
#' @param object An object
#' @param ... `r .dot_param`
#'
#' @name ObjFile-IO
#' @export prepInfo
prepInfo <- function(object, ...) {
  UseMethod(generic = "prepInfo", object = object)
}

## coerce ######################################################################

#' Read and write the OBJECT.json file
#'
#' The `OBJECT.json` file provides the metadata of the object that is
#' represented by the corresponding directory.
#'
#' @param object An object
#' @param ... `r .dot_param`
#'
#' @name coerce-objects
NULL

#' @rdname coerce-objects
#' @export toSCE
toSCE <- function(object, ...) {
  UseMethod(generic = "toSCE", object = object)
}

#' @rdname coerce-objects
#' @export toSCME
toSCME <- function(object, ...) {
  UseMethod(generic = "toSCME", object = object)
}

#' @rdname coerce-objects
#' @export toSeurat
toSeurat <- function(object, ...) {
  UseMethod(generic = "toSeurat", object = object)
}

## Export to H5AD and H5MU #####################################################

#' Write an object to H5AD (H5MU) file
#'
#' @param object An object
#' @param file Path of the H5AD (H5MU) file
#' @param ... `r .dot_param` Mainly to `prepH5AD` and `prepH5MU`
#'
#' @name exportH5AD
NULL

#' @rdname exportH5AD
#' @export exportH5AD
exportH5AD <- function(object, ...) {
  UseMethod(generic = "exportH5AD", object = object)
}

#' @rdname exportH5AD
#' @export exportH5AD
prepH5AD <- function(object, ...) {
  UseMethod(generic = "prepH5AD", object = object)
}

#' @rdname exportH5AD
#' @export exportH5MU
exportH5MU <- function(object, ...) {
  UseMethod(generic = "exportH5MU", object = object)
}

#' @rdname exportH5AD
#' @export exportH5MU
prepH5MU <- function(object, ...) {
  UseMethod(generic = "prepH5MU", object = object)
}

## ChromatinAssay5 #############################################################

as.ChromatinAssay5 <- function(x, ...) {
  UseMethod(generic = "as.ChromatinAssay5", object = x)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 Generics ##################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Rename functions ############################################################

#' @export
#' @rdname set-dimnames
setGeneric(
  "setRownames",
  function(x, new.names, ...) standardGeneric("setRownames")
)

#' @export
#' @rdname set-dimnames
setGeneric(
  "setColnames",
  function(x, new.names, ...) standardGeneric("setColnames")
)

## SingleCellExperiment ########################################################

#' Fetch cellular data
#'
#' Fetch data for a set of observations (columns) in an object
#'
#' @param object A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' @param ... Arguments to be passed to other methods
#'
#' @return
#' A \code{\link[base]{data.frame}} with cells as rows and \code{vars} data as
#' columns.
#'
#' @export
#' @rdname fetchData
setGeneric(
  "fetchData",
  function(object, ...) standardGeneric("fetchData")
)

#' @export
#' @rdname fetchData
setGeneric(
  "variableFeatures",
  function(object, ...) standardGeneric("variableFeatures")
)

#' @export
#' @rdname fetchData
setGeneric(
  "variableFeatures<-",
  function(object, ..., value) standardGeneric("variableFeatures<-")
)

## ChromatinExperiment #########################################################

#' @export
#' @rdname fragments
setGeneric(
  "fragments",
  function(x, ...) standardGeneric("fragments")
)

#' @export
#' @rdname fragments
setGeneric(
  "fragments<-",
  function(x, ..., value) standardGeneric("fragments<-")
)

#' @export
#' @rdname annotations
setGeneric(
  "annotations",
  function(x, ...) standardGeneric("annotations")
)

#' @export
#' @rdname annotations
setGeneric(
  "annotations<-",
  function(x, ..., value) standardGeneric("annotations<-")
)

## MultiAssayExperiment ########################################################

#' @export
#' @rdname experiment
setGeneric(
  "experiment",
  function(x, e, ...) standardGeneric("experiment")
)

## SingleCellMultiExperiment ###################################################

#' @export
#' @rdname SingleCellMultiExperiment-internal
setGeneric(
  "int_SCE",
  function(x, ...) standardGeneric("int_SCE")
)

#' @export
#' @rdname SingleCellMultiExperiment-internal
setGeneric(
  "int_SCE<-",
  function(x, ..., value) standardGeneric("int_SCE<-")
)

#' @export
#' @rdname defaultExp
setGeneric(
  "defaultExp",
  function(x, ...) standardGeneric("defaultExp")
)

#' @export
#' @rdname defaultExp
setGeneric(
  "defaultExp<-",
  function(x, ..., value) standardGeneric("defaultExp<-")
)

