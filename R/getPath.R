
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S3 methods ###################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @examples
#' ### For merged BPCells matrix
#' bpce <- load_example_sce()
#' getPath(counts(bpce))
#'
#' @export
#' @rdname getPath
#' @method getPath IterableMatrix
getPath.IterableMatrix <- function(object, ...) {
  return(.get_path_for_bpcells(object, ...))
}

#' @examples
#' ### For merged BPCells fragments
#' cbpce <- load_example_csce()
#' getPath(counts(cbpce))
#' getPath(fragments(cbpce))
#'
#' ### For merged BPCells object
#' cbpce2 <- load_example_csce("unsorted")
#' cbpce <- merge(cbpce, cbpce2)
#' getPath(counts(cbpce))
#' getPath(fragments(cbpce))
#'
#' @export
#' @rdname getPath
#' @method getPath IterableFragments
getPath.IterableFragments <- function(object, ...) {
  return(.get_path_for_bpcells(object, ...))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Get path for BPCells objects ################################################

#' Fetch all paths of backend files for BPCells object
#'
#' @importFrom hdf5r.Extra h5AbsLinkName
#' @noRd
.get_path_for_bpcells <- function(object, ...) {
  slots <- slotNames(object)
  use.slots <- c("matrix_list", "fragments_list", "path", "dir")
  empty.out <- data.frame(
    path = character(),
    type = character(),
    group = character()
  )
  out <- list(empty.out)
  if (!any(use.slots %in% slots)) {
    for (i in slots) {
      slot.obj <- slot(object, name = i)
      is.iter_frag <- inherits(slot.obj, "IterableFragments")
      is.iter_mat <- inherits(slot.obj, "IterableMatrix")
      if (is.iter_mat || is.iter_frag) {
        out[[i]] <- .get_path_for_bpcells(slot.obj, ...)
      }
    }
    names(out) <- NULL
    return(do.call(rbind, out))
  }
  if ("path" %in% slots) {
    path <- slot(object, name = "path") %>%
      file_path_as_absolute()
    group <- NA
    type <- "tsv"
    if (inherits(object, "10xMatrixH5")) {
      type <- "h5"
    }
    if ("group" %in% slots) {
      group <- h5AbsLinkName(name = slot(object, name = "group"))
      type <- "h5"
    }
    return(data.frame(path = path, type = type, group = group))
  }
  if ("dir" %in% slots) {
    path <- slot(object, name = "dir") %>%
      file_path_as_absolute()
    return(data.frame(path = path, type = "dir", group = NA))
  }
  iter_obj.list <- list()
  if ("fragments_list" %in% slots) {
    iter_obj.list <- c(
      iter_obj.list,
      slot(object, name = "fragments_list")
    )
  }
  if ("matrix_list" %in% slots) {
    iter_obj.list <- c(
      iter_obj.list,
      slot(object, name = "matrix_list")
    )
  }
  for (i in seq_along(iter_obj.list)) {
    out[[i]] <- .get_path_for_bpcells(iter_obj.list[[i]], ...)
  }
  names(x = out) <- NULL
  return(do.call(rbind, out, ...))
}

#' Check if a path has been occupied by BPCells object
#'
#' @param x The BPCells object to be check.
#' @param path Path to be check.
#' @param name Name of HDF5 link to be check.
#'
#' @noRd
#' @importFrom hdf5r.Extra h5AbsLinkName
.check_bpcells_file_occupy <- function(x, path, name = NULL) {
  path.use <- .get_path_for_bpcells(x)
  path <- normalizePath(path = path, mustWork = FALSE)
  name <- name %iff% h5AbsLinkName(name)
  path.use <- path.use[path.use$path %in% path, , drop = FALSE]
  if (any(name %in% path.use$group)) {
    stop(
      "\n  The destination has been occupied ",
      "by the iterable object to be written:",
      "\n  Path: ", path,
      "\n  Group: ", name
    )
  }
  if (nrow(path.use) > 0) {
    stop(
      "\n  The destination has been occupied ",
      "by the iterable object to be written:",
      "\n  Path: ", path
    )
  }
  return(invisible(NULL))
}
