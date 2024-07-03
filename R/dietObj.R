
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S3 methods ###################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param assays Select \code{\link[SummarizedExperiment]{assays}} specified
#' here. If is \code{NA}, remove all assays.
#' @param verbose `r .vb_param`
#'
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @rdname dietObj
#' @export
#' @method dietObj SummarizedExperiment
dietObj.SummarizedExperiment <- function(
    object,
    assays = NULL,
    verbose = TRUE,
    ...
) {
  object <- formatObj(object)
  object <- .diet_sce_field(
    object = object,
    func = "assays",
    name.func = "assayNames",
    ns = "SummarizedExperiment",
    names = assays,
    verbose = verbose
  )
  return(object)
}

#' @param reducedDims Keep a subset of dimension reductions specified here. If
#' is \code{NA}, remove all \code{\link[SingleCellExperiment]{reducedDims}}.
#' @param colPairs Keep a subset of \code{\link[SingleCellExperiment]{colPairs}}
#' specified here. If is \code{NA}, remove all \code{colPairs}.
#' @param rowPairs Similar as \code{colPairs}.
#' @param altExps Keep a subset of \code{\link[SingleCellExperiment]{altExps}}
#' specified here. If is \code{NA}, remove all \code{altExps}.
#'
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @rdname dietObj
#' @export
#' @method dietObj SingleCellExperiment
dietObj.SingleCellExperiment <- function(
    object,
    assays = NULL,
    reducedDims = NULL,
    colPairs = NULL,
    rowPairs = NULL,
    altExps = NULL,
    verbose = TRUE,
    ...
) {
  old_func <- getS3method(f = "dietObj", class = "SummarizedExperiment")
  object <- old_func(object, assays = assays, verbose = verbose, ...)
  fields <- list(
    reducedDims = reducedDims,
    colPairs = colPairs,
    rowPairs = rowPairs,
    altExps = altExps
  )
  for (i in names(fields)) {
    name.func <- paste0(sub("s$", "", i), "Names")
    object <- .diet_sce_field(
      object = object,
      func = i,
      name.func = name.func,
      ns = "SingleCellExperiment",
      names = fields[[i]],
      verbose = verbose
    )
  }
  return(object)
}

#' @param fragments Logical scalar specifying whether to keep
#' \code{\link{fragments}}.
#' @param annotations Logical scalar specifying whether to keep
#' \code{\link{annotations}}.
#' @param seqinfo Logical scalar specifying whether to keep
#' \code{\link{seqinfo}}.
#'
#' @rdname dietObj
#' @export
#' @method dietObj ChromExperiment
dietObj.ChromExperiment <- function(
    object,
    assays = NULL,
    reducedDims = NULL,
    colPairs = NULL,
    rowPairs = NULL,
    altExps = NULL,
    fragments = NULL,
    annotations = NULL,
    seqinfo = NULL,
    verbose = TRUE,
    ...
) {
  old_func <- getS3method(f = "dietObj", class = "SingleCellExperiment")
  out <- old_func(
    object = object,
    assays = assays,
    reducedDims = reducedDims,
    colPairs = colPairs,
    rowPairs = rowPairs,
    altExps = altExps,
    verbose = verbose,
    ...
  )
  if (is_false(fragments)) {
    verboseMsg("Remove fragements")
    fragments(object) <- NULL
  }
  if (is_false(annotations)) {
    verboseMsg("Remove annotations")
    annotations(object) <- NULL
  }
  if (is_false(seqinfo)) {
    verboseMsg("Remove seqinfo")
    seqinfo(object) <- NULL
  }
  return(object)
}

#' @param experiments Select \code{\link[MultiAssayExperiments]{experiments}}
#' specified here. Cannot be \code{NA} (remove all assays).
#'
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#'
#' @rdname dietObj
#' @export
#' @method dietObj MultiAssayExperiment
dietObj.MultiAssayExperiment <- function(
    object,
    experiments = NULL,
    verbose = TRUE,
    ...
) {
  if (length(experiments) == 0) {
    return(object)
  }
  if (is.na(experiments)) {
    stop("Cannot remove all experiments")
  }
  old.data <- experiments(object)
  old.names <- names
  if (is_bare_numeric(experiments)) {
    experiments <- intersect(experiments, seq_along(old.data))
  }
  if (is_bare_character(experiments)) {
    experiments <- intersect(experiments, names(object))
  }
  if (length(x = experiments) == 0) {
    warning(
      "None of the following experiments is found: \n  ", old.names,
      call. = FALSE, immediate. = TRUE
    )
    return(object)
  }
  verboseMsg("Keep experiments: ", paste(experiments, collapse = ", "))
  experiments(object) <- old.data[experiments]
  return(object)
}

#' @rdname dietObj
#' @export
#' @method dietObj SingleCellMultiExperiment
dietObj.SingleCellMultiExperiment <- function(
    object,
    experiments = NULL,
    reducedDims = NULL,
    colPairs = NULL,
    verbose = TRUE,
    ...
) {
  int_sce <- formatObj(int_SCE(object))
  old_func <- getS3method(f = "dietObj", class = "ChromExperiment")
  object <- old_func(object, experiments = experiments, verbose = verbose, ...)
  fields <- list(reducedDims = reducedDims, colPairs = colPairs)
  for (i in names(fields)) {
    name.func <- paste0(sub("s$", "", i), "Names")
    int_sce <- .diet_sce_field(
      object = int_sce,
      func = i,
      name.func = name.func,
      ns = "SingleCellExperiment",
      names = fields[[i]],
      verbose = verbose
    )
  }
  int_SCE(object) <- int_sce
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom rlang is_bare_numeric
.diet_sce_field <- function(
    object,
    func,
    name.func,
    ns,
    names = NULL,
    verbose = TRUE
) {
  get.func <- getExportedValue(ns = ns, name = func)
  set.func <- getExportedValue(ns = ns, name = paste0(func, "<-"))
  getname.func <- getExportedValue(ns = ns, name = name.func)

  if (length(names) == 0) {
    verboseMsg("Skip ", func)
    return(object)
  }
  if (is.na(names)) {
    verboseMsg("Remove all ", func)
    return(set.func(object, value = list()))
  }
  old.data <- get.func(object)
  old.names <- names
  if (is_bare_numeric(names)) {
    names <- intersect(names, seq_along(old.data))
  }
  if (is_bare_character(names)) {
    curr.names <- getname.func(object)
    names <- intersect(names, curr.names)
  }
  if (length(x = names) == 0) {
    warning(
      "None of the following ", func, " is found: \n  ", old.names,
      call. = FALSE, immediate. = TRUE
    )
    return(object)
  }
  verboseMsg("Keep ", func, ": ", paste(names, collapse = ", "))
  return(set.func(object, value = old.data[names]))
}
