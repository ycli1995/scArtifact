
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Constructor ##################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The ChromExperiment class
#'
#' The \code{ChromExperiment} class is designed to represent single-cell
#' chromatin accessibility (such as ATAC-seq) data using \pkg{BPCells}. It
#' inherits from the `r .doc_links("SingleCellExperiment")` class and is used in
#' the same manner. In addition, it supports storage of barcoded and aligned
#' fragments information via \code{\link{fragments}}, in which the
#' `r .doc_links("IterableFragments")` objects are used as the underlying data.
#'
#' @param ... Arguments passed to the constructor function of
#' `r .doc_links("SingleCellExperiment")` to fill the slots of the base class.
#' @param sep Separators to use for strings encoding genomic coordinates. The
#' first element is used to separate the chromosome from the coordinates, the
#' second element is used to separate the start from end coordinate. Only used
#' if `ranges` is `NULL`.
#' @param ranges A set of `r .doc_links("GRanges")` corresponding to the rows
#' (peaks) of the input matrix.
#' @param genome A `r .doc_links("Seqinfo")` object containing basic information
#' about the genome used. Alternatively, the name of a UCSC genome
#' can be provided and the sequence information will be downloaded from UCSC.
#' @param annotations A `r .doc_links("GRanges")` containing annotations for
#' the genome used
#' @param fragments Fragments data for the input matrix. Can be one of the
#' following:
#' \itemize{
#' \item Tabix-indexed fragment files (*.tsv.gz from 10x Genomics).
#' \item Directories created by `r .doc_links("write_fragments_dir")`.
#' \item HDF5 files created by `r .doc_links("write_fragments_hdf5")`.
#' \item One or a list of `r .doc_links("IterableFragments")` objects. Note that
#' if the `r .doc_links("IterableFragments")` is created by opening the 10x
#' `*.tsv.gz` using `r .doc_links("open_fragments_10x")` directly, the cell
#' names cannot be accessed or modified by `r .doc_links("cellNames")`.
#' }
#'
#' @details
#' In this class, rows should represent genomic coordinates (e.g., peaks) while
#' columns represent samples generated from single cells.
#'
#' The extra arguments in the constructor (e.g., `sep`, `ranges`, `genome`,
#' `annotations` and \code{\link{fragments}}) represent the main extensions
#' implemented in the \code{ChromExperiment} class. Readers are referred to
#' the specific documentation pages for more details.
#'
#' @return
#' A `ChromExperiment` object
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#'
#' @export
#' @docType class
ChromExperiment <- function(
    ...,
    sep = c("-", "-"),
    ranges = NULL,
    genome = NULL,
    annotations = NULL,
    fragments = NULL
) {
  sce <- SingleCellExperiment(...)
  .sce_to_csce(
    sce = sce,
    sep = sep,
    ranges = ranges,
    genome = genome,
    annotations = annotations,
    fragments = fragments
  )
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods ###################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## fragments ###################################################################

#' Get or set fragments data
#'
#' Methods to get or set fragments data for a \code{\link{ChromExperiment}}
#' object, to hold information needed for working with fragment files.
#'
#' @param x A \code{\link{ChromExperiment}} object.
#' @param ... Passed to \code{\link{openFragments}}
#'
#' @note Currently we only support storing fragments as a single
#' `r .doc_links("IterableFragments")` object. If the input object inherits from
#' `r .doc_links("Fragment")`, it will be converted into a
#' `r .doc_links("IterableFragments")` via \code{\link{openFragments}}.
#'
#' @returns
#' \itemize{
#' \item `fragments`: Retrieve the stored `r .doc_links("IterableFragments")`.
#' \item `fragments<-`: `x` with updated `fragments`
#' }
#'
#' @seealso
#' \itemize{
#' \item `r .doc_links("IterableFragments")`
#' \item `r .doc_links("Fragment")`
#' }
#'
#' @name fragments
NULL

#' @importFrom SingleCellExperiment int_metadata
#' @export
#' @rdname fragments
setMethod(
  f = "fragments",
  signature = "ChromExperiment",
  definition = function(x, ...) {
    return(int_metadata(x)[[.frag_key]])
  }
)

#' @param value \itemize{
#' \item A list of `r .doc_links("IterableFragments")` or
#' `r .doc_links("Fragment")` objects.
#' \item A single `r .doc_links("IterableFragments")`  or
#' `r .doc_links("Fragment")` object.
#' \item A string specifying the path to a tabix-indexed fragments file.
#' \item `NULL`. This will remove all `fragments` data.
#' }
#'
#' @importFrom SingleCellExperiment int_metadata<-
#' @export
#' @rdname fragments
setMethod(
  f = "fragments<-",
  signature = c("ChromExperiment", "IterableFragments"),
  definition = function(x, ..., value) {
    frag.cells <- cellNames(value)
    cells <- colnames(x)
    if (!any(frag.cells %in% cells)) {
      warning(
        "fragments<-: No cell for the input IterableFragments is found in 'x'",
        immediate. = TRUE, call. = FALSE
      )
      return(x)
    }
    value <- value %>%
      select_cells(fastIntersect(frag.cells, cells)) %>%
      select_chromosomes(seqlevels(rowRanges(x)))
    cells.notfound <- setdiff(cells, cellNames(value))
    if (length(cells.notfound) > 0) {
      warning(
        "fragments<-: ", length(cells.notfound), " cells were not found ",
        "in the input IterableFragments:\n  ",
        paste(head(cells.notfound), collapse = ", "), "...",
        immediate. = TRUE, call. = FALSE
      )
    }
    int_metadata(x)[[.frag_key]] <- value
    return(x)
  }
)

#' @importClassesFrom Signac Fragment
#' @export
#' @rdname fragments
setMethod(
  f = "fragments<-",
  signature = c("ChromExperiment", "Fragment"),
  definition = function(x, ..., value) {
    fragments(x) <- openFragments(value, ...)
    return(x)
  }
)

#' @export
#' @rdname fragments
setMethod(
  f = "fragments<-",
  signature = c("ChromExperiment", "character"),
  definition = function(x, ..., value) {
    frags <- list()
    for (i in seq_along(value)) {
      frags[[i]] <- openFragments(value[i], ...)
    }
    if (length(frags) == 0) {
      return(x)
    }
    if (length(frags) == 1) {
      fragments(x) <- frags[[1]]
      return(x)
    }
    fragments(x) <- Reduce(f = c, x = frags)
    return(x)
  }
)

#' @importClassesFrom S4Vectors list_OR_List
#' @export
#' @rdname fragments
setMethod(
  f = "fragments<-",
  signature = c("ChromExperiment", "list_OR_List"),
  definition = function(x, ..., value) {
    cells <- colnames(x)
    cells.loaded <- character()
    for (i in seq_along(value)) {
      check_inherits_for_func(
        value = value[[i]],
        classes = .support_frgas,
        func = "fragments<-",
        name = "value"
      )
      value[[i]] <- openFragments(value[[i]], ...)
      cell_selection <- cellNames(value[[i]]) %>%
        setdiff(cells.loaded) %>%
        fastIntersect(x = cells)
      value[[i]] <- select_cells(value[[i]], cell_selection)
      cells.loaded <- c(cells.loaded, cellNames(value[[i]]))
    }
    cells.notfound <- setdiff(cells, cells.loaded)
    if (length(cells.notfound) > 0) {
      warning(
        length(cells.notfound), " cells not found in loaded fragments:\n  ",
        paste(head(cells.notfound), collapse = ", "), "...",
        immediate. = TRUE
      )
    }
    if (length(value) == 0) {
      return(x)
    }
    if (length(value) == 1) {
      fragments(x) <- value[[1]]
      return(x)
    }
    fragments(x) <- Reduce(f = c, x = value)
    return(x)
  }
)

#' @importFrom SingleCellExperiment int_metadata<-
#' @export
#' @rdname fragments
setMethod(
  f = "fragments<-",
  signature = c("ChromExperiment", "NULL"),
  definition = function(x, ..., value) {
    int_metadata(x)[[.frag_key]] <- NULL
    return(x)
  }
)

## annotations #################################################################

#' Get or set genomic annotation
#'
#' Methods to get or set a `r .doc_links("GRanges")` specifying the genomic
#' annotation in a \code{\link{ChromExperiment}}.
#'
#' @param x A \code{\link{ChromExperiment}} object.
#' @param ... Arguments passed to other methods.
#'
#' @seealso
#' `r .doc_links("Annotation")`
#'
#' @name annotations
NULL

#' @importFrom SingleCellExperiment int_metadata
#' @export
#' @rdname annotations
setMethod(
  f = "annotations",
  signature = "ChromExperiment",
  definition = function(x, ...) {
    return(int_metadata(x)[[.annot_key]])
  }
)

#' @param value \itemize{
#' \item A \code{GRanges} object.
#' \item \code{NULL}. This will remove the existing \code{annotations}
#' }
#'
#' @importFrom SingleCellExperiment int_metadata<-
#' @importClassesFrom GenomicRanges GRanges
#' @export
#' @rdname annotations
setMethod(
  f = "annotations<-",
  signature = c("ChromExperiment", "GRanges"),
  definition = function(x, ..., value) {
    current.genome <- genome(x) %>%
      unique()
    annot.genome <- genome(value) %>%
      unique()
    if (!is.null(current.genome)) {
      if (!is.na(annot.genome) & (current.genome != annot.genome)) {
        stop("annotations genome does not match genome of the object")
      }
    }
    int_metadata(x)[[.annot_key]] <- value
    return(x)
  }
)

#' @importFrom SingleCellExperiment int_metadata<-
#' @export
#' @rdname annotations
setMethod(
  f = "annotations<-",
  signature = c("ChromExperiment", "NULL"),
  definition = function(x, ..., value) {
    int_metadata(x)[[.annot_key]] <- value
    return(x)
  }
)

## seqinfo #####################################################################

#' Access and modify sequence information for ChromExperiment
#'
#' Methods for accessing and modifying the `r .doc_links("Seqinfo")` object
#' stored in a \code{\link{ChromExperiment}} object.
#'
#' @param x A \code{\link{ChromExperiment}} object.
#' @param value A `r .doc_links("Seqinfo")` object or name of a UCSC genome to
#' store in the \code{\link{ChromExperiment}}
#'
#' @details
#' These methods are intended to be consistent with methods for
#' `r .doc_links("ChromatinAssay")` in the \pkg{Signac} package.
#'
#' @seealso
#' \itemize{
#' \item \code{\link[GenomeInfoDb]{seqinfo}} in the \pkg{GenomeInfoDb} package.
#' \item \code{\link[Signac]{seqinfo-methods}}
#' }
#'
#' @name seqinfo-methods
NULL

#' @importFrom SingleCellExperiment int_metadata
#' @importFrom GenomeInfoDb seqinfo
#' @export
#' @rdname seqinfo-methods
setMethod(
  f = "seqinfo",
  signature = "ChromExperiment",
  definition = function(x) {
    return(int_metadata(x)[[.sinfo_key]])
  }
)

#' @importFrom SingleCellExperiment int_metadata<-
#' @importFrom GenomeInfoDb seqinfo<- Seqinfo
#' @importClassesFrom GenomeInfoDb Seqinfo
#' @export
#' @rdname seqinfo-methods
setMethod(
  f = "seqinfo<-",
  signature = "ChromExperiment",
  definition = function(x, value) {
    if (inherits(value, what = "Seqinfo")) {
      int_metadata(x)[[.sinfo_key]] <- value
    } else if (is(value, "character")) {
      int_metadata(x)[[.sinfo_key]] <- Seqinfo(genome = value)
    } else if (is.null(value)) {
      int_metadata(x)[[.sinfo_key]] <- NULL
    } else {
      stop(
        "Unknown object supplied. Choose a Seqinfo object ",
        "or the name of a UCSC genome"
      )
    }
    return(x)
  }
)

## seqlevels ###################################################################

#' @importFrom GenomeInfoDb seqlevels
#' @export
#' @rdname seqinfo-methods
setMethod(
  f = "seqlevels",
  signature = "ChromExperiment",
  definition = function(x) {
    x <- seqinfo(x)
    if (is.null(x)) {
      return(NULL)
    }
    callGeneric()
  }
)

#' @importFrom GenomeInfoDb seqlevels<-
#' @export
#' @rdname seqinfo-methods
setMethod(
  f = "seqlevels<-",
  signature = "ChromExperiment",
  definition = function(x, value) {
    sinfo <- seqinfo(x)
    seqlevels(sinfo) <- value
    seqinfo(x) <- sinfo
    return(x)
  }
)

## seqnames ####################################################################

#' @importFrom GenomeInfoDb seqnames
#' @export
#' @rdname seqinfo-methods
setMethod(
  f = "seqnames",
  signature = "ChromExperiment",
  definition = function(x) {
    x <- seqinfo(x)
    if (is.null(x)) {
      return(NULL)
    }
    callGeneric()
  }
)

#' @importFrom GenomeInfoDb seqnames<-
#' @export
#' @rdname seqinfo-methods
setMethod(
  f = "seqnames<-",
  signature = "ChromExperiment",
  definition = function(x, value) {
    sinfo <- seqinfo(x)
    seqnames(sinfo) <- value
    seqinfo(x) <- sinfo
    return(x)
  }
)

## seqlengths ##################################################################

#' @importFrom GenomeInfoDb seqlengths
#' @export
#' @rdname seqinfo-methods
setMethod(
  f = "seqlengths",
  signature = "ChromExperiment",
  definition = function(x) {
    x <- seqinfo(x)
    if (is.null(x)) {
      return(NULL)
    }
    callGeneric()
  }
)

#' @importFrom GenomeInfoDb seqlengths<-
#' @export
#' @rdname seqinfo-methods
setMethod(
  f = "seqlengths<-",
  signature = "ChromExperiment",
  definition = function(x, value) {
    sinfo <- seqinfo(x)
    seqlengths(sinfo) <- value
    seqinfo(x) <- sinfo
    return(x)
  }
)

## genome ######################################################################

#' @importFrom GenomeInfoDb genome
#' @export
#' @rdname seqinfo-methods
setMethod(
  f = "genome",
  signature = "ChromExperiment",
  definition = function(x) {
    x <- seqinfo(x)
    if (is.null(x)) {
      return(NULL)
    }
    callGeneric()
  }
)

#' @importFrom GenomeInfoDb genome<-
#' @export
#' @rdname seqinfo-methods
setMethod(
  f = "genome<-",
  signature = "ChromExperiment",
  definition = function(x, value) {
    seqinfo(x) <- value
    return(x)
  }
)

## isCircular ##################################################################

#' @importFrom GenomeInfoDb isCircular
#' @export
#' @rdname seqinfo-methods
setMethod(
  f = "isCircular",
  signature = "ChromExperiment",
  definition = function(x) {
    x <- seqinfo(x)
    if (is.null(x)) {
      return(NULL)
    }
    callGeneric()
  }
)

#' @importFrom GenomeInfoDb isCircular<-
#' @export
#' @rdname seqinfo-methods
setMethod(
  f = "isCircular<-",
  signature = "ChromExperiment",
  definition = function(x, value) {
    sinfo <- seqinfo(x)
    isCircular(sinfo) <- value
    seqinfo(x) <- sinfo
    return(x)
  }
)

## Show ########################################################################

#' @param object A `ChromExperiment` object
#'
#' @importFrom SingleCellExperiment altExpNames mainExpName reducedDimNames
#' @importFrom S4Vectors coolcat
#' @importFrom methods selectMethod show
#' @export
#' @rdname ChromExperiment
setMethod(
  f = "show",
  signature = "ChromExperiment",
  definition = function(object) {
    .old_show_sce(object)
    cat("fragments: \n")
    if (!is.null(x = fragments(object))) {
      frag_paths <- getPath(object = fragments(object))
      for (i in seq_len(nrow(frag_paths))) {
        cat(" path:", frag_paths[i, 1], "\n")
        if (!is.na(x = frag_paths[i, 3])) {
          cat(" group:", frag_paths[i, 3], "\n")
        }
      }
    }
    .show_backend_path(object)
    return(invisible(x = NULL))
  }
)

## Subset ######################################################################

#' Subset ChromExperiment object
#'
#' Deal with the dimensions and subset of \code{\link{ChromExperiment}} object.
#' These methods are generally identical to those inherited operations, with
#' additional manipulation for `fragments(x)`.
#'
#' @param x A \code{\link{ChromExperiment}}.
#' @param ... Arguments passed to other methods.
#'
#' @name ChromExperiment-subset
NULL

#' @param value An object of a class specified in the S4 method signature.
#'
#' @examples
#' csce <- load_example_csce()
#' csce
#'
#' colnames(csce) <- paste0("Sorted_", colnames(csce))
#' csce
#' fragments(csce)
#'
#' @export
#' @rdname ChromExperiment-subset
setMethod(
  f = "dimnames<-",
  signature = c("ChromExperiment", "list"),
  definition = function(x, value) {
    old.colnames <- colnames(x)
    x <- callNextMethod()
    if (identical(old.colnames, value[[2L]])) {
      return(x)
    }
    if (is.null(fragments(x))) {
      return(x)
    }
    old.cells.frags <- cellNames(fragments(x))
    match.idx <- fmatch(old.cells.frags, old.colnames)
    cellNames(fragments(x)) <- value[[2L]][match.idx]
    return(x)
  }
)

#' @param i,j Indices specifying rows or columns to extract.
#' @param drop For matrices and arrays. If TRUE the result is coerced to the
#' lowest possible dimension.
#'
#' @examples
#' csce[1:10, 1:20]
#'
#' csce[, sample(colnames(csce), 20)]
#'
#' @export
#' @rdname ChromExperiment-subset
setMethod(
  f = "[",
  signature = c("ChromExperiment", "ANY", "ANY"),
  definition = function(x, i, j, ..., drop = TRUE) {
    old.frags <- fragments(x)
    x <- callNextMethod()
    if (ncol(x) == 0) {
      fragments(x) <- NULL
      return(x)
    }
    if (is.null(old.frags)) {
      return(x)
    }
    fragments(x) <- select_cells(
      fragments = old.frags,
      cell_selection = colnames(x)
    )
    return(x)
  }
)

## Validity ####################################################################

#' @importFrom SummarizedExperiment rowRanges
#' @importClassesFrom GenomicRanges GRanges
.valid_csce <- function(object) {
  msg <- NULL
  if (!inherits(x = rowRanges(x = object), what = "GRanges")) {
    msg <- msg %>%
      c(paste("'rowRanges' of", class(x = object), "must be GRanges"))
  }
  valid.frags <- is.null(fragments(object)) ||
    inherits(fragments(object), what = "IterableFragments")
  if (!valid.frags) {
    msg <- msg %>%
      c(paste("'fragments' of", class(object), "must be IterableFragments"))
  }
  if (length(x = msg) > 0) {
    return(msg)
  }
  return(TRUE)
}

#' @importFrom S4Vectors setValidity2
setValidity2(Class = "ChromExperiment", method = .valid_csce)

## Coerce ######################################################################

setAs(
  from = "SingleCellExperiment",
  to = "ChromExperiment",
  def = function(from) {
    ranges <- rowRanges(x = from)
    if (!inherits(x = ranges, what = "GenomicRanges")) {
      ranges <- NULL
    }
    fragments <- int_metadata(from)[[.frag_key]]
    genome <- int_metadata(from)[[.sinfo_key]]
    annotations <- int_metadata(from)[[.annot_key]]
    return(.sce_to_csce(
      sce = from,
      sep = c(":", "-"),
      ranges = ranges,
      genome = genome,
      annotations = annotations,
      fragments = fragments
    ))
  }
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Coerce SingleCellExperiment to ChromExperiment
#'
#' @importFrom SummarizedExperiment rowData rowData<- rowRanges<-
#' @importFrom GenomicRanges granges
#' @importFrom Signac StringToGRanges
#'
#' @noRd
.sce_to_csce <- function(
    sce,
    sep = c("-", "-"),
    ranges = NULL,
    genome = NULL,
    annotations = NULL,
    fragments = NULL,
    ...
) {
  # old <- S4Vectors:::disableValidity()
  # if (!isTRUE(x = old)) {
  #   S4Vectors:::disableValidity(disabled = TRUE)
  #   on.exit(S4Vectors:::disableValidity(disabled = old))
  # }
  ranges <- ranges %||% StringToGRanges(regions = rownames(sce), sep = sep)
  if (length(ranges) != nrow(sce)) {
    stop(
      "Length of 'ranges' does not match number of rows",
      " in SingleCellExperiment"
    )
  }
  if (!isDisjoint(ranges)) {
    warning(
      "Overlapping 'ranges' supplied. Ranges should be non-overlapping.",
      immediate. = TRUE, call. = FALSE
    )
  }
  row.data <- rowData(sce)
  rowRanges(sce) <- ranges
  rownames(sce) <- rownames(row.data)
  rowData(sce) <- row.data
  csce <- new(Class = "ChromExperiment", sce)
  csce <- .format_csce_rownames(csce = csce)
  fragments(csce) <- fragments
  annotations(csce) <- annotations
  genome(csce) <- genome
  return(csce)
}
