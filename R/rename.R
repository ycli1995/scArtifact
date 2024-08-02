
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods ######################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Rename the dimensions
#'
#' Methods to set new dimension names to an object. These helpers are mostly
#' designed for `r .doc_links("MultiAssayExperiment")` and its sub-classes.
#'
#' @param x A matrix-like object or object inheriting from
#' `r .doc_links("MultiAssayExperiment")`
#' @param new.names Strings to set as the new dimension names. Can be a single
#' vector of characters or a list containing multiple vectors.
#' @param i For a `r .doc_links("MultiAssayExperiment")` object, should be
#' integers or characters specifying the experiment(s) to be renamed.
#' \itemize{
#' \item When `new.names` is a list, the length of `i` should be the same as
#' `new.names`. If `i` is `NA` or missing, all experiments will be used by
#' default.
#' \item When `new.names` is a character vector, `i` should also be a single
#' integer or character. If `i` is `NA` or missing, will change the primary
#' dimension names of `x`.
#' }
#' @param ... `r .dot_param`
#'
#' @name set-dimnames
NULL

## setRownames #################################################################

#' @export
#' @rdname set-dimnames
setMethod(
  f = "setRownames",
  definition = function(x, new.names, ...) {
    rownames(x) <- new.names
    return(x)
  }
)

#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @importClassesFrom S4Vectors list_OR_List
#' @export
#' @rdname set-dimnames
setMethod(
  f = "setRownames",
  signature = c("MultiAssayExperiment", "list_OR_List"),
  definition = function(x, new.names, i, ...) {
    if (missing(i)) {
      x <- .set_dimname_mae(
        x = x,
        i = names(x),
        new.names = new.names,
        type = "feature"
      )
      return(x)
    }
    x <- .set_dimname_mae(
      x = x,
      i = i,
      new.names = new.names,
      type = "feature"
    )
    return(x)
  }
)

#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @export
#' @rdname set-dimnames
setMethod(
  f = "setRownames",
  signature = c("MultiAssayExperiment", "character"),
  definition = function(x, new.names, i, ...) {
    if (missing(i)) {
      i <- NA
    }
    if (inherits(i, "Character_OR_Numeric") || is.na(i)) {
      x <- .set_dimname_mae(
        x = x,
        i = i,
        new.names = new.names,
        type = "feature"
      )
    }
    return(x)
  }
)

## setColnames #################################################################

#' @export
#' @rdname set-dimnames
setMethod(
  f = "setColnames",
  definition = function(x, new.names, ...) {
    colnames(x) <- new.names
    return(x)
  }
)

#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @importClassesFrom S4Vectors list_OR_List
#' @export
#' @rdname set-dimnames
setMethod(
  f = "setColnames",
  signature = c("MultiAssayExperiment", "list_OR_List"),
  definition = function(x, new.names, i, ...) {
    if (missing(i)) {
      x <- .set_dimname_mae(
        x = x,
        i = names(x),
        new.names = new.names,
        type = "sample"
      )
      return(x)
    }
    if (inherits(i, "Character_OR_Numeric")) {
      x <- .set_dimname_mae(
        x = x,
        i = i,
        new.names = new.names,
        type = "sample"
      )
    }
    return(x)
  }
)

#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @export
#' @rdname set-dimnames
setMethod(
  f = "setColnames",
  signature = c("MultiAssayExperiment", "character"),
  definition = function(x, new.names, i, ...) {
    if (missing(i)) {
      i <- NA
    }
    if (inherits(i, "Character_OR_Numeric") || is.na(i)) {
      x <- .set_dimname_mae(
        x = x,
        i = i,
        new.names = new.names,
        type = "sample"
      )
    }
    return(x)
  }
)

#' @importFrom BiocBaseUtils setSlots
#' @export
#' @rdname set-dimnames
setMethod(
  f = "setColnames",
  signature = c("SingleCellMultiExperiment", "character"),
  definition = function(x, new.names, i, ...) {
    if (missing(i)) {
      i <- NA
    }
    if (is.na(i)) {
      int_sce <- int_SCE(x)
      default.exp <- defaultExp(x)
      x <- .set_dimname_mae(
        x = as(x, Class = "MultiAssayExperiment"),
        i = i,
        new.names = new.names,
        type = "sample"
      )
      colnames(int_sce) <- rownames(colData(x))
      return(new(
        Class = "SingleCellMultiExperiment",
        x,
        int_SCE = int_sce,
        defaultExp = default.exp
      ))
    }
    callNextMethod(x = x, new.names = new.names, i = i)
  }
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom fastmatch fmatch
#' @importFrom BiocBaseUtils setSlots
#' @importFrom MultiAssayExperiment experiments listToMap mapToList
#' renamePrimary
.set_dimname_mae <- function(x, i, new.names, type = c("sample", "feature")) {
  type <- match.arg(type)
  if (any(is.na(i))) {
    # set primary names
    if (type == "sample") {
      stopifnot(is.character(new.names))
      x <- renamePrimary(x = x, value = new.names)
    }
    return(x)
  }
  stopifnot(length(i) == length(new.names))
  margin <- switch(EXPR = type, "feature" = 1L, "sample" = 2L)
  exps <- experiments(x)
  dmap <- mapToList(dfmap = sampleMap(x))
  for (j in seq_along(i)) {
    exp.idx <- i[j]

    if (type == "sample") {
      idx <- fmatch(
        dmap[[exp.idx]][['colname']],
        dimnames(exps[[exp.idx]])[[margin]]
      )
      dmap[[exp.idx]][['colname']] <- new.names[[j]][idx]
    }
    dimnames(exps[[exp.idx]])[[margin]] <- new.names[[j]]
  }
  dmap <- listToMap(listmap = dmap)
  x <- setSlots(x, sampleMap = dmap, ExperimentList = exps)
  return(x)
}

.check_set_rownames <- function(x, new.names) {
  if (is_true(all.equal(new.names, rownames(x), check.attributes = FALSE))) {
    return(x)
  }
  if (length(new.names) == 0) {
    return(x)
  }
  rownames(x) <- new.names
  return(x)
}

.check_set_colnames <- function(x, new.names) {
  if (is_true(all.equal(new.names, colnames(x), check.attributes = FALSE))) {
    return(x)
  }
  if (length(new.names) == 0) {
    return(x)
  }
  colnames(x) <- new.names
  return(x)
}
