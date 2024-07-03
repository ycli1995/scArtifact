
#' @rdname Obj-IO
#' @export
#' @method reloadObj default
reloadObj.default <- function(object, path, verbose = TRUE, ...) {
  return(object)
}

#' @rdname Obj-IO
#' @export
#' @method reloadObj list
reloadObj.list <- function(object, path, verbose = TRUE, ...) {
  nms <- names(object)
  for (i in seq_along(object)) {
    nm <- gsub("\\s|\\/", "_", nms[i])
    filepath <- file.path(path, nm)
    verboseMsg("Reload '", nms[i], "' from ", filepath)
    object[[i]] <- reloadObj(object[[i]], path = filepath, verbose = verbose)
  }
  return(object)
}

#' @importClassesFrom S4Vectors SimpleList
#'
#' @rdname Obj-IO
#' @export
#' @method reloadObj SimpleList
reloadObj.SimpleList <- function(object, path, verbose = TRUE, ...) {
  old_func <- getS3method(f = "reloadObj", class = "list")
  old_func(object, path = path, verbose = verbose, ...)
}

#' @rdname Obj-IO
#' @export
#' @method reloadObj data.frame
reloadObj.data.frame <- function(object, path, verbose = TRUE, ...) {
  return(object)
}

#' @importClassesFrom S4Vectors DataFrame
#'
#' @rdname Obj-IO
#' @export
#' @method reloadObj DataFrame
reloadObj.DataFrame <- function(object, path, verbose = TRUE, ...) {
  return(object)
}

#' @rdname Obj-IO
#' @export
#' @method reloadObj matrix
reloadObj.matrix <- function(
    object,
    path,
    verbose = TRUE,
    to.mem = FALSE,
    ...
) {
  if (to.mem) {
    return(object)
  }
  old_func <- getS3method(f = "reloadObj", class = "DelayedArray")
  old_func(object = object, path = path, verbose = verbose, ...)
}

#' @importClassesFrom DelayedArray DelayedArray
#'
#' @rdname Obj-IO
#' @export
#' @method reloadObj matrix
reloadObj.DelayedArray <- function(object, path, verbose = TRUE, ...) {
  verboseMsg("Reload HDF5Array from ", path)
  extra <- readObjFile(path)
  HDF5Array(
    filepath = file.path(path, extra$file_name),
    name = extra$name,
    ...
  )
}

#' @rdname Obj-IO
#' @export
#' @method reloadObj IterableMatrix
reloadObj.IterableMatrix <- function(object, path, verbose = TRUE, ...) {
  verboseMsg("Reload MatrixDir from ", path)
  extra <- readObjFile(path)
  open_matrix_dir(dir = file.path(path, extra$dir_name), ...)
}

#' @rdname Obj-IO
#' @export
#' @method reloadObj IterableFragments
reloadObj.IterableFragments <- function(object, path, verbose = TRUE, ...) {
  verboseMsg("Reload FragmentsDir from ", path)
  extra <- readObjFile(path)
  open_fragments_dir(dir = file.path(path, extra$dir_name), ...)
}

#' @rdname Obj-IO
#' @export
#' @method reloadObj dgCMatrix
reloadObj.dgCMatrix <- function(
    object,
    path,
    verbose = TRUE,
    to.mem = FALSE,
    ...
) {
  if (to.mem) {
    return(object)
  }
  old_func <- getS3method(f = "reloadObj", class = "IterableMatrix")
  old_func(object = object, path = path, verbose = verbose, ...)
}

#' @importFrom SummarizedExperiment assays assays<-
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#'
#' @rdname Obj-IO
#' @export
#' @method reloadObj RangedSummarizedExperiment
reloadObj.RangedSummarizedExperiment <- function(
    object,
    path,
    verbose = TRUE,
    ...
) {
  verboseMsg("Reload ", class(object)[1], " from ", path)
  assays(object) <- reloadObj(
    assays(object),
    path = file.path(path, "assays"),
    verbose = verbose
  )
  return(object)
}

#' @importFrom SingleCellExperiment altExps altExps<-
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#'
#' @rdname Obj-IO
#' @export
#' @method reloadObj SingleCellExperiment
reloadObj.SingleCellExperiment <- function(object, path, verbose = TRUE, ...) {
  old_func <- getS3method(f = "reloadObj", class = "RangedSummarizedExperiment")
  object <- old_func(object, path = path, verbose = verbose, ...)
  altExps(object) <- reloadObj(
    altExps(object),
    path = file.path(path, "altExps"),
    verbose = verbose
  )
  int_metadata(object)[[.path_key]] <- file_path_as_absolute(path)
  return(object)
}

#' @rdname Obj-IO
#' @export
#' @method reloadObj ChromExperiment
reloadObj.ChromExperiment <- function(object, path, verbose = TRUE, ...) {
  old_func <- getS3method(f = "reloadObj", class = "SingleCellExperiment")
  object <- old_func(object, path = path, verbose = verbose, ...)
  fragments(object) <- reloadObj(
    fragments(object),
    path = file.path(path, "fragments"),
    verbose = verbose
  )
  return(object)
}

#' @importFrom MultiAssayExperiment experiments experiments<-
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#'
#' @rdname Obj-IO
#' @export
#' @method reloadObj MultiAssayExperiment
reloadObj.MultiAssayExperiment <- function(object, path, verbose = TRUE, ...) {
  verboseMsg("Reload ", class(object)[1], " from ", path)
  experiments(object) <- reloadObj(
    experiments(object),
    path = file.path(path, "experiments"),
    verbose = verbose
  )
  return(object)
}

#' @rdname Obj-IO
#' @export
#' @method reloadObj SingleCellMultiExperiment
reloadObj.SingleCellMultiExperiment <- function(object, path, verbose = TRUE, ...) {
  verboseMsg("Reload ", class(object)[1], " from ", path)
  experiments(object) <- reloadObj(
    experiments(object),
    path = file.path(path, "experiments"),
    verbose = verbose
  )
  int_metadata(object)[[.path_key]] <- file_path_as_absolute(path)
  return(object)
}
