
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S3 Methods ###################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param path Path where `object` to be saved.
#' @param overwrite Whether or not to overwrite the existing `path`.
#' @param verbose `r .vb_param`
#'
#' @rdname Obj-IO
#' @export
#' @method writeObj default
writeObj.default <- function(
    object,
    path,
    overwrite = FALSE,
    verbose = TRUE,
    ...
) {
  .simple_writeObj(
    object = object,
    path = path,
    overwrite = overwrite,
    verbose = verbose,
    ...
  )
  return(reloadObj(object, path = path, verbose = verbose))
}

#' @rdname Obj-IO
#' @export
#' @method writeObj list
writeObj.list <- function(
    object,
    path,
    overwrite = FALSE,
    verbose = TRUE,
    ...
) {
  nms <- names(object)
  for (i in seq_along(object)) {
    if (!any(nchar(nms[i]))) {
      names(object)[i] <- as.character(i - 1L)
    }
  }
  .simple_writeObj(
    object = object,
    path = path,
    overwrite = overwrite,
    verbose = verbose,
    ...
  )
  return(reloadObj(object, path = path, verbose = verbose))
}

#' @importClassesFrom S4Vectors SimpleList
#'
#' @rdname Obj-IO
#' @export
#' @method writeObj SimpleList
writeObj.SimpleList <- function(
    object,
    path,
    overwrite = FALSE,
    verbose = TRUE,
    ...
) {
  old_func <- getS3method(f = "writeObj", class = "list")
  old_func(
    object,
    path = path,
    overwrite = overwrite,
    verbose = verbose,
    ...
  )
}

#' @param to.mem Whether to keep the matrix in RAM. If `FALSE`, will use
#' the on-disk format:
#' \itemize{
#' \item `r .doc_links("IterableMatrix")` for sparse matrix
#' \item `r .doc_links("HDF5Array")` for dense matrix
#' }
#'
#' @importFrom HDF5Array HDF5Array writeHDF5Array
#'
#' @rdname Obj-IO
#' @export
#' @method writeObj matrix
writeObj.matrix <- function(
    object,
    path,
    overwrite = FALSE,
    verbose = TRUE,
    to.mem = TRUE,
    ...
) {
  .simple_writeObj(
    object = object,
    path = path,
    overwrite = overwrite,
    verbose = verbose,
    ...
  )
  return(reloadObj(object, path = path, verbose = verbose, to.mem = to.mem))
}

#' @rdname Obj-IO
#' @export
#' @method writeObj CsparseMatrix
writeObj.CsparseMatrix <- function(
    object,
    path,
    overwrite = FALSE,
    verbose = TRUE,
    to.mem = TRUE,
    ...
) {
  .simple_writeObj(
    object = object,
    path = path,
    overwrite = overwrite,
    verbose = verbose,
    ...
  )
  return(reloadObj(object, path = path, verbose = verbose, to.mem = to.mem))
}

#' @rdname Obj-IO
#' @export
#' @method writeObj TsparseMatrix
writeObj.TsparseMatrix <- function(
    object,
    path,
    overwrite = FALSE,
    verbose = TRUE,
    to.mem = TRUE,
    ...
) {
  .simple_writeObj(
    object = object,
    path = path,
    overwrite = overwrite,
    verbose = verbose,
    ...
  )
  return(reloadObj(object, path = path, verbose = verbose, to.mem = to.mem))
}

#' @rdname Obj-IO
#' @export
#' @method writeObj IterableMatrix
writeObj.IterableMatrix <- function(
    object,
    path,
    overwrite = FALSE,
    verbose = TRUE,
    ...
) {
  if (!is_sparse(object)) {
    warning(
      "The ", class(object)[1], " is dense, which is not supported for ",
      "BPCells::write_matrix_dir() currently. Skip writeObj.",
      immediate. = TRUE, call. = FALSE
    )
    return(object)
  }
  writeObj.default(
    object = object,
    path = path,
    overwrite = overwrite,
    verbose = verbose,
    ...
  )
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.writeobj_ver <- "1.0"

.overwrite_obj_path <- function(path, overwrite = FALSE, verbose = TRUE) {
  if (!dir.exists(path)) {
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    return(file_path_as_absolute(path))
  }
  files <- list.dirs(path = path, full.names = FALSE, recursive = FALSE)
  if (length(files) == 0) {
    return(file_path_as_absolute(path))
  }
  if (is_false(x = overwrite)) {
    stop(path, " was already written. Please set 'overwrite = TRUE'.")
  }
  verboseMsg("Overwriting existing path '", path, "'")
  tmp.path <- tempfile(tmpdir = dirname(path))
  dir.create(path = tmp.path, showWarnings = FALSE, recursive = TRUE)
  return(file_path_as_absolute(tmp.path))
}

.overwrite_rename <- function(tmp.path, path) {
  path <- file_path_as_absolute(path)
  tmp.path <- file_path_as_absolute(tmp.path)
  if (!identical(tmp.path, path)) {
    unlink(path, recursive = TRUE, force = TRUE)
    file.rename(tmp.path, path)
  }
  return(invisible(NULL))
}

.simple_writeObj <- function(
    object,
    path,
    overwrite = FALSE,
    verbose = TRUE,
    ...
) {
  tmp.path <- .overwrite_obj_path(path, overwrite, verbose)
  path <- file_path_as_absolute(path)
  tmp.path <- file_path_as_absolute(tmp.path)
  if (!identical(tmp.path, path)) {
    on.exit(unlink(tmp.path, force = TRUE, recursive = TRUE))
  }
  dumpObj(object, path = tmp.path, verbose = verbose, ...)
  .overwrite_rename(tmp.path, path)
  return(invisible(NULL))
}
