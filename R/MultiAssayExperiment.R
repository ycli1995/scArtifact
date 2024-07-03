
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions ####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Extracting an assay with primaries as column names
#'
#' Change the column names of extracted experiment with primaries stored in the
#' `r .doc_links("colData")` of `r .doc_links("MultiAssayExperiment")`.
#'
#' @param x A `r .doc_links("MultiAssayExperiment")` object.
#' @param i Which experiment to be fetched. Can be an integer or character.
#'
#' @return The extracted experiment object.
#'
#' @importFrom MultiAssayExperiment experiments mapToList sampleMap
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @export
getWithPrimaries <- function(x, i) {
  if (!inherits(x, "MultiAssayExperiment")) {
    stop("Provide a MultiAssayExperiment as input")
  }
  stopifnot(
    is.numeric(i) || is.character(i),
    identical(length(i), 1L),
    !is.na(i),
    !is.logical(i)
  )
  exp <- experiments(x)[[i]]
  sampMap <- mapToList(dfmap = sampleMap(x))[[i]]
  nameMap <- setNames(sampMap[['primary']], sampMap[['colname']])
  colnames(exp) <- nameMap[colnames(exp)]
  return(exp)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods ###################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Accessing single experiment data in MultiAssayExperiment
#'
#' Supplemental methods for accessing single experiment from
#' `r .doc_links("MultiAssayExperiment")` objects.
#'
#' @param x A `r .doc_links("MultiAssayExperiment")` object
#' @param e Which experiment to fetch. Can be a integer index or an experiment
#' name. If missing, will fetch the first experiment by default.
#' @param ... Arguments passed to other methods.
#'
#' @return
#' The extracted experiment object
#'
#' @seealso `r .doc_links("experiments")`
#'
#' @name experiment
NULL

#' @param withPrimaries Whether or not to use primary names to replace the
#' original column names
#' @param withColData Whether or not to also fetch the associated
#' `r .doc_links("colData")` of `x`.
#' @param mode String indicating how `r .doc_links("MultiAssayExperiment")`
#' metadata should be added or just replaced to the extracted experiment.
#' Passed to `r .doc_links("getWithColData")`.
#'
#' @importFrom MultiAssayExperiment experiments getWithColData mapToList
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @export
#' @rdname experiment
setMethod(
  f = "experiment",
  signature = c("MultiAssayExperiment", "character"),
  definition = function(
    x,
    e,
    withPrimaries = FALSE,
    withColData = FALSE,
    mode = c("append", "replace"),
    ...
  ) {
    mode <- match.arg(mode)
    if (!withColData && !withPrimaries) {
      return(experiments(x)[[e]])
    }
    if (!withColData && withPrimaries) {
      return(getWithPrimaries(x, i = e))
    }
    # MultiAssayExperiment::getWithColData will automatically set colnames to
    # primaries.
    return(tryCatch(
      expr = getWithColData(x = x, i = e, mode = mode),
      error = function(err) {
        message(
          "MultiAssayExperiment::getWithColData:\n ", err, "\n",
          "Force to set 'withColData = FALSE'"
        )
        if (withPrimaries) {
          return(getWithPrimaries(x = x, i = e))
        }
        return(experiments(x)[[e]])
      }
    ))
  }
)

#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @export
#' @rdname experiment
setMethod(
  f = "experiment",
  signature = c("MultiAssayExperiment", "numeric"),
  definition = function(
    x,
    e,
    withPrimaries = FALSE,
    withColData = FALSE,
    mode = c("append", "replace"),
    ...
  ) {
    mode <- match.arg(mode)
    exp.name <- names(x)[[e]]
    return(experiment(
      x = x,
      e = exp.name,
      withPrimaries = withPrimaries,
      withColData = withColData,
      mode = mode,
      ...
    ))
  }
)

#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @export
#' @rdname experiment
setMethod(
  f = "experiment",
  signature = c("MultiAssayExperiment", "missing"),
  definition = function(
    x,
    e,
    withPrimaries = FALSE,
    withColData = FALSE,
    mode = c("append", "replace"),
    ...
  ) {
    mode <- match.arg(mode)
    exp.name <- names(x)[[1]]
    return(experiment(
      x = x,
      e = exp.name,
      withPrimaries = withPrimaries,
      withColData = withColData,
      mode = mode,
      ...
    ))
  }
)
