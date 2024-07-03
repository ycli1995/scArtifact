
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods ###################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## fetchData ###################################################################

#' @param vars Vector of all variables to fetch.
#' @param cells Cells to collect data for (default is all cells).
#' @param assay The `r .doc_links("assay")` to collect data from.
#' @param use.Exp Which `r .doc_links("altExp")` to collect data from. If the
#' selected `r .doc_links("altExp")` doesn't exist, use the main Experiment.
#'
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom SingleCellExperiment altExp altExpNames reducedDims
#' @importFrom SeuratObject EmptyDF
#' @importFrom easy.utils fastIntersect fetchColnames
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @export
#' @rdname fetchData
setMethod(
  f = "fetchData",
  signature = "SingleCellExperiment",
  definition = function(
    object,
    vars = NULL,
    cells = NULL,
    assay = NULL,
    use.Exp = NULL
  ) {
    cells <- cells %||% colnames(object)
    if (is.numeric(cells)) {
      cells <- colnames(object)[cells]
    }
    cells <- fetchColnames(object, query = cells)
    data.fetched <- data.frame(row.names = cells)
    if (length(vars) == 0) {
      return(data.fetched)
    }
    # Check altExp in object
    if (length(use.Exp) > 1) {
      warning(
        "Only use the first one of 'use.Exp': ", use.Exp[1],
        call. = FALSE, immediate. = TRUE
      )
      use.Exp <- use.Exp[1]
    }
    if (length(use.Exp) == 0) {
      assay <- assay %||% .select_default_assay(object)
    }
    if (any(use.Exp %in% altExpNames(object))) {
      alt.exp <- altExp(object, e = use.Exp)
      assay <- assay %||% .select_default_assay(alt.exp)
      assay.data <- assay(alt.exp, i = assay)
    } else {
      if (length(use.Exp) > 0) {
        warning(
          "Ignore the non-existing Exp '", use.Exp, "'",
          call. = FALSE, immediate. = TRUE
        )
        use.Exp <- NULL
      }
      assay <- assay %||% .select_default_assay(object)
      assay.data <- assay(object, i = assay)
    }

    raw.vars <- vars
    # Find vars in assays
    assay.vars <- fastIntersect(vars, rownames(assay.data))
    if (length(assay.vars) > 0) {
      assay.fetched <- assay.data[assay.vars, cells, drop = FALSE] %>%
        as.matrix()
      data.fetched <- cbind(data.fetched, t(assay.fetched))
      vars <- setdiff(vars, assay.vars)
    }
    # Find vars in colData
    colData.vars <- fastIntersect(vars, names(colData(object)))
    if (length(colData.vars) > 0) {
      colData.fetched <- colData(object)[cells, colData.vars, drop = FALSE] %>%
        as.data.frame(optional = TRUE)
      data.fetched <- cbind(data.fetched, colData.fetched)
      vars <- setdiff(vars, colData.vars)
    }
    # Find all vars in reducedDims(object)
    # (The same as Keys in SeuratObject)
    all.reducs <- reducedDims(object)
    for (i in names(all.reducs)) {
      use.reduc <- all.reducs[[i]]
      # First find vars already in colnames
      rd.vars <- fastIntersect(vars, colnames(use.reduc))
      if (length(rd.vars) > 0) {
        rd.fetched <- use.reduc[cells, rd.vars, drop = FALSE] %>%
          as.matrix()
        data.fetched <- cbind(data.fetched, rd.fetched)
      }
      vars <- setdiff(vars, rd.vars)
      # Then find vars in format 'Key_1', 'Key_2'...
      colnames(use.reduc) <- paste0(i, "_", 1:ncol(use.reduc))
      rd.vars <- fastIntersect(vars, colnames(use.reduc))
      if (length(rd.vars) > 0) {
        rd.fetched <- use.reduc[cells, rd.vars, drop = FALSE] %>%
          as.matrix()
        data.fetched <- cbind(data.fetched, rd.fetched)
      }
      vars <- setdiff(vars, rd.vars)
    }
    if (length(vars) > 0) {
      warning(
        "The following query variables were not found: \n  ",
        paste(vars, collapse = ", "),
        immediate. = TRUE
      )
    }
    final.vars <- fastIntersect(raw.vars, colnames(data.fetched))
    data.fetched <- data.fetched[, final.vars, drop = FALSE]
    return(data.fetched)
  }
)

#' @importFrom SummarizedExperiment rowData
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @export
#' @rdname fetchData
setMethod(
  f = "variableFeatures",
  signature = "SingleCellExperiment",
  definition = function(object, ...) {
    return(rownames(object)[rowData(object)[['variable']]])
  }
)

#' @param value An object of a class specified in the S4 method signature.
#'
#' @importFrom SummarizedExperiment rowData<-
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @export
#' @rdname fetchData
setMethod(
  f = "variableFeatures<-",
  signature = c("SingleCellExperiment", "character"),
  definition = function(object, ..., value) {
    rowData(object)[['variable']] <- rownames(object) %in% value
    return(object)
  }
)

#' @importFrom SummarizedExperiment rowData<-
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @export
#' @rdname fetchData
setMethod(
  f = "variableFeatures<-",
  signature = c("SingleCellExperiment", "NULL"),
  definition = function(object, ..., value) {
    rowData(object)[['variable']] <- rownames(object) %in% value
    return(object)
  }
)

## Backend path ################################################################

#' Get the backend path for SingleCellExperiment
#'
#' @param object A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' @param ... Additional arguments, for use in specific methods.
#'
#' @name SCE-backend
NULL

#' @importFrom BiocGenerics path
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @export
#' @rdname SCE-backend
setMethod(
  f = "path",
  signature = "SingleCellExperiment",
  definition = function(object, ...) {
    return(int_metadata(object)[[.path_key]])
  }
)

## show ########################################################################

#' @importFrom methods getMethod show
#' @export
#' @rdname SCE-backend
setMethod(
  f = "show",
  signature = "SingleCellExperiment",
  definition = function(object) {
    .old_show_sce(object)
    .show_backend_path(object)
  }
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.old_show_sce <- function(object) {
  old_show <- selectMethod(f = "show", signature = "SummarizedExperiment")
  old_show(object)
  coolcat(fmt = "reducedDimNames(%d): %s\n", vals = reducedDimNames(object))
  me <- mainExpName(object)
  if (is.null(me)) {
    me <- "NULL"
  }
  cat(sprintf(fmt = "mainExpName: %s\n", me))
  coolcat(fmt = "altExpNames(%d): %s\n", vals = altExpNames(object))
}

#' @importFrom BiocGenerics path
.show_backend_path <- function(object) {
  p <- path(object)
  cat("Backend path: ")
  cat(p[1], "\n")
  return(invisible(NULL))
}

#' Select default assay from SummarizedExperiment-like object
#'
#' Default order of priority is: logcounts, normcounts, counts, scaled
#'
#' @importFrom SummarizedExperiment assayNames
#' @noRd
.select_default_assay <- function(x) {
  default_assays <- c("logcounts", "normcounts", "counts", "scaled")
  all_assays <- assayNames(x)
  ret_assay <- intersect(default_assays, all_assays)
  if (length(ret_assay) > 0) {
    return(ret_assay[1])
  }
  stop(
    "\n  The SCE didn't contain any of the following assay(s): ",
    paste(default_assays, collapse = ", "), "\n",
    "\n  Cannot automatically select default assay."
  )
}
