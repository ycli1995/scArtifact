
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S3 methods ###################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## openFragments ###############################################################

#' @examples
#' frags_file <- system.file(
#'   "extdata", "pbmc_sorted_fragments.tsv.gz",
#'   package = "BPCellsExperiment"
#' )
#' tmp.file <- tempfile("pbmc_sorted_fragments", fileext = ".tsv.gz")
#' file.copy(frags_file, tmp.file)
#'
#' ## Open using file name
#' frags <- openFragments(tmp.file)
#'
#' @importFrom hdf5r is_hdf5
#' @importFrom hdf5r.Extra h5AbsLinkName
#' @importFrom dplyr filter select
#' @importFrom tibble add_row
#' @importFrom tools md5sum
#' @export
#' @rdname openFragments
#' @method openFragments character
openFragments.character <- function(x, group = NULL, ...) {
  group <- group %||% "fragments"
  group <- h5AbsLinkName(name = group)
  path <- file_path_as_absolute(x)
  is.dir <- dir.exists(paths = path)
  is.tsv <- grepl(pattern = "tsv(\\.gz)?$", x = path)
  is.h5 <- is_hdf5(name = path)
  if (is.h5) {
    return(open_fragments_hdf5(path = path, group = group, ...))
  }
  if (is.dir) {
    return(open_fragments_dir(dir = path, ...))
  }
  if (!is.tsv) {
    stop("'x' must be an H5 file, a directory or an indexed .tsv")
  }
  myhash <- md5sum(files = path) %>%
    paste(collapse = "_")
  if (myhash %in% .BPCells.envs$frags_filemap$hash) {
    # The query .tsv file has been cached into a directory, just open it.
    file.map <- .BPCells.envs$frags_filemap %>%
      filter(hash == myhash) %>%
      select(bpcells_dir)
    return(open_fragments_dir(dir = file.map$bpcells_dir, ...))
  }
  # Cache the query .tsv into an H5 file
  tmp.file <- basename(path) %>%
    gsub(pattern = "\\.tsv(\\.gz)?$", replacement = ".") %>%
    tempfile(tmpdir = tempdir()) %>%
    normalizePath(mustWork = FALSE)
  frags <- open_fragments_10x(path = path, ...) %>%
    write_fragments_dir(dir = tmp.file, overwrite = FALSE)
  .BPCells.envs$frags_filemap <- .BPCells.envs$frags_filemap %>%
    add_row(
      tsv = path,
      hash = myhash,
      bpcells_dir = file_path_as_absolute(x = tmp.file)
    )
  return(frags)
}

#' @examples
#' file.copy(paste0(frags_file, ".tbi"), paste0(tmp.file, ".tbi"))
#'
#' ## Open using Signac::Fragment
#' frags <- Signac::CreateFragmentObject(tmp.file, cells = cellNames(frags))
#' frags <- openFragments(frags)
#'
#' @importFrom rlang is_empty
#' @importClassesFrom Signac Fragment
#' @export
#' @rdname openFragments
#' @method openFragments Fragment
openFragments.Fragment <- function(x, group = NULL, ...) {
  file <- slot(x, name = "path")
  cells <- slot(x, name = "cells")
  frags <- openFragments(x = file, group = NULL, ...)
  if (!is_empty(cells)) {
    frags <- select_cells(frags, cell_selection = cells)
    cellNames(frags) <- names(cells)
  }
  return(frags)
}

#' @examples
#' ## Open an already existing Fragment
#' frags <- openFragments(frags)
#'
#' @export
#' @rdname openFragments
#' @method openFragments IterableFragments
openFragments.IterableFragments <- function(x, group = NULL, ...) {
  return(x)
}

#' @section Check directories for opened fragments:
#' When \code{openFragments} meets a `r .doc_links("Fragments")` object or
#' a `.tsv.gz` file, it will first check whether there is already a directory
#' created by `r .doc_links("write_fragments_dir")` and cached the fragments
#' data to open. If so, it will directly open the fragments. Otherwise, it will
#' calculate the \code{\link[tools]{md5sum}} of the `.tsv.gz`, then cache it
#' into a new directory. Use `showFragmentsFileMap()` to get the file mapping
#' table created by current R session.
#'
#' @export
#' @rdname openFragments
showFragmentsFileMap <- function() {
  .BPCells.envs$frags_filemap
}

#' @section Clear temporary fragments files:
#' When all needed fragments data in current R session have been stored in
#' proper destination paths, \code{clearFragmentsFileMap()} can be called to
#' delete all the temporary h5 files shown by \code{showFragmentsFileMap()}, so
#' that the disk storage can be saved.
#'
#' @export
#' @rdname openFragments
clearFragmentsFileMap <- function() {
  unlink(
    .BPCells.envs$frags_filemap$bpcells_dir,
    force = TRUE,
    recursive = TRUE
  )
  .init_frags_filemap()
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Fragments file map ##########################################################

.init_frags_filemap <- function() {
  .BPCells.envs$frags_filemap <- data.frame(
    tsv = character(),
    hash = character(),
    bpcells_dir = character()
  )
  return(invisible(x = NULL))
}
