
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S3 methods ###################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @export
#' @method formatObj default
formatObj.default <- function(object, ...) {
  return(object)
}

#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @rdname formatObj
#' @export
#' @method formatObj SummarizedExperiment
formatObj.SummarizedExperiment <- function(object, ...) {
  .format_sce_names(object, func = "assays", ns = "SummarizedExperiment")
}

#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @rdname formatObj
#' @export
#' @method formatObj SingleCellExperiment
formatObj.SingleCellExperiment <- function(object, ...) {
  object <- .format_sce_names(
    object,
    func = "assays",
    ns = "SummarizedExperiment"
  )
  fields <- c("colPairs", "rowPairs", "altExps")
  for (f in fields) {
    object <- .format_sce_names(object, func = f,  ns = "SingleCellExperiment")
  }
  object <- .format_sce_reducs(object)
  return(object)
}

#' @rdname formatObj
#' @export
#' @method formatObj ChromExperiment
formatObj.ChromExperiment <- function(object, ...) {
  object <- .format_csce_rownames(object)
  old_func <- getS3method(f = "formatObj", class = "SingleCellExperiment")
  object <- old_func(object, ...)
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.format_sce_names <- function(sce, func, ns) {
  pfx <- sub("s$", "", func)
  get.func <- getExportedValue(ns = ns, name = func)
  set.func <- getExportedValue(ns = ns, name = paste0(func, "<-"))
  old.data <- get.func(sce)
  old.names <- names(old.data)
  new.data <- list()
  for (i in seq_along(old.data)) {
    nm <- old.names[i]
    if (length(nm) == 0) {
      nm <- paste0(pfx, i)
    }
    new.data[[nm]] <- old.data[[i]]
  }
  sce <- set.func(sce, value = new.data)
  return(sce)
}

#' Format the structure of reducedDims(SCE)
#'
#' Make sure the \code{reducedDims(sce)} to be \code{LinearEmbeddingMatrix}
#'
#' @importFrom SingleCellExperiment LinearEmbeddingMatrix reducedDims
#' reducedDims<-
#' @importClassesFrom SingleCellExperiment LinearEmbeddingMatrix
#' @noRd
.format_sce_reducs <- function(sce) {
  old.reds <- reducedDims(sce, withDimnames = TRUE)
  red.names <- names(old.reds)
  new.reds <- list()
  for (i in seq_along(old.reds)) {
    nm <- red.names[i]
    if (length(nm) == 0) {
      nm <- paste0("reducedDim", i)
    }
    tmp.reds <- old.reds[[i]]
    if (!inherits(tmp.reds, "LinearEmbeddingMatrix")) {
      tmp.reds <- tmp.reds %>%
        as.matrix() %>%
        as("LinearEmbeddingMatrix")
    }
    new.reds[[nm]] <- tmp.reds
  }
  reducedDims(sce, withDimnames = TRUE) <- new.reds
  return(sce)
}

#' Format the row names of ChromExperiment
#'
#' Make sure the row names of \code{ChromExperiment} are formatted as
#' \code{"chr:start-end"}
#'
#' @importFrom SummarizedExperiment rowRanges
#'
#' @noRd
.format_csce_rownames <- function(csce) {
  rownames(csce) <- as.character(rowRanges(csce))
  return(csce)
}
