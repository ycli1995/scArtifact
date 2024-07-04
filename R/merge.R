
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 Methods ###################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Merge single-cell data objects
#'
#' Methods to combine multiple single-cell data objects. The behaviors of
#' merging are the same as \code{\link[SeuratObject]{merge}} function from
#' \pkg{Seurat}.
#'
#' @param x A \code{SingleCellExperiment} or \code{SingleCellMultiExperiment}
#' object.
#' @param y A single object or a list of objects. The class of each object
#' should be the same as \code{x}.
#' @param ... Arguments passed to other methods.
#' @param assays Names of \code{\link[SummarizedExperiment]{assays}} to be
#' merged. If set to \code{NULL}, common assays will be automatically searched
#' from input objects. If set to \code{FALSE}, no assay will be merged.
#' @param reducedDims Names of \code{\link[SingleCellExperiment]{reducedDims}}
#' to be merged. If set to \code{NULL}, common reductions will be automatically
#' searched from input objects. If set to \code{FALSE}, no reduction will be
#' merged.
#' @param altExps Names of \code{\link[SingleCellExperiment]{altExps}} to be
#' merged. If set to \code{NULL}, common alternative experiments will be
#' automatically searched from input objects. If set to \code{FALSE}, no
#' \code{altExps} will be merged.
#' @param altExps.assays Names of \code{assays} to be merged in \code{altExps}.
#' @param experiments For \code{SingleCellMultiExperiment}, names of
#' \code{\link[MultiAssayExperiment]{experiments}} to be merged. If set to
#' \code{NULL}, all experiments will be used.
#' @param global.reducedDims For \code{SingleCellMultiExperiment}, names of the
#' global \code{reducedDims} to be merged.
#' @param label Add a new column to \code{\link[SummarizedExperiment]{colData}}
#' of each input to hold the batch information.
#' @param add.idx A character vector of \code{length(x = c(x, y))}; appends the
#' corresponding values to the start of each objects' cell names.
#' @param idx.collapse String to collapse \code{add.idx} and cell names of each
#' input object.
#' @param verbose Display progress.
#'
#' @return
#' A merged object
#'
#' @aliases merge
#' merge,SingleCellExperiment,SingleCellExperiment-method
#' merge,SingleCellExperiment,list_OR_List-method
#' merge,ChromExperiment,ChromExperiment-method
#' merge,ChromExperiment,list_OR_List-method
#'
#' @name SCE-merge
NULL

## SingleCellExperiment ########################################################

#' @examples
#' # BPCExperiment
#' bpce1 <- load_example_sce()
#' bpce2 <- load_example_sce(dataset = "unsorted")
#' bpce <- merge(bpce1, bpce2, add.idx = c("sorted", "unsorted"))
#' bpce
#'
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @export
#' @rdname SCE-merge
setMethod(
  f = "merge",
  signature = c("SingleCellExperiment", "SingleCellExperiment"),
  definition = function(
    x,
    y,
    assays = NULL,
    reducedDims = NULL,
    altExps = NULL,
    altExps.assays = NULL,
    label = NULL,
    add.idx = NULL,
    idx.collapse = "_",
    verbose = TRUE,
    ...
  ) {
    return(merge(
      x = x,
      y = list(y),
      assays = assays,
      reducedDims = reducedDims,
      altExps = altExps,
      altExps.assays = altExps.assays,
      label = label,
      add.idx = add.idx,
      idx.collapse = idx.collapse,
      verbose = verbose,
      ...
    ))
  }
)

#' @importClassesFrom S4Vectors list_OR_List
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @export
#' @rdname SCE-merge
setMethod(
  f = "merge",
  signature = c("SingleCellExperiment", "list_OR_List"),
  definition = function(
    x,
    y,
    assays = NULL,
    reducedDims = NULL,
    altExps = NULL,
    altExps.assays = NULL,
    label = NULL,
    add.idx = NULL,
    idx.collapse = "_",
    verbose = TRUE,
    ...
  ) {
    check_inherits_for_vec(x = y, class = "SingleCellExperiment", name = "y")
    SCEs <- c(list(x), as.list(x = y))
    return(.merge_SCEs(
      SCEs = SCEs,
      assays = assays,
      reducedDims = reducedDims,
      altExps = altExps,
      altExps.assays = altExps.assays,
      label = label,
      add.idx = add.idx,
      idx.collapse = idx.collapse,
      verbose = verbose,
      ...
    ))
  }
)

## ChromExperiment #############################################################

#' @examples
#' # ChromExperiment
#' cbpce1 <- load_example_csce()
#' cbpce2 <- load_example_csce(dataset = "unsorted")
#' cbpce <- merge(cbpce1, cbpce2, add.idx = c("sorted", "unsorted"))
#' cbpce
#'
#' @export
#' @rdname SCE-merge
setMethod(
  f = "merge",
  signature = c("ChromExperiment", "ChromExperiment"),
  definition = function(
    x,
    y,
    assays = NULL,
    reducedDims = NULL,
    altExps = NULL,
    altExps.assays = NULL,
    label = NULL,
    add.idx = NULL,
    idx.collapse = "_",
    verbose = TRUE,
    ...
  ) {
    return(merge(
      x = x,
      y = list(y),
      assays = assays,
      reducedDims = reducedDims,
      altExps = altExps,
      altExps.assays = altExps.assays,
      label = label,
      add.idx = add.idx,
      idx.collapse = idx.collapse,
      verbose = verbose,
      ...
    ))
  }
)

#' @importClassesFrom S4Vectors list_OR_List
#' @export
#' @rdname SCE-merge
setMethod(
  f = "merge",
  signature = c("ChromExperiment", "list_OR_List"),
  definition = function(
    x,
    y,
    assays = NULL,
    reducedDims = NULL,
    altExps = NULL,
    altExps.assays = NULL,
    label = NULL,
    add.idx = NULL,
    idx.collapse = "_",
    verbose = TRUE,
    ...
  ) {
    check_inherits_for_vec(x = y, class = "ChromExperiment", name = "y")
    SCEs <- c(list(x), as.list(x = y))
    return(.merge_ChromSCEs(
      SCEs = SCEs,
      assays = assays,
      reducedDims = reducedDims,
      altExps = altExps,
      altExps.assays = altExps.assays,
      label = label,
      add.idx = add.idx,
      idx.collapse = idx.collapse,
      verbose = verbose,
      ...
    ))
  }
)

## SingleCellMultiExperiment ###################################################

#' @examples
#' # SingleCellMultiExperiment
#' scme1 <- load_example_scme()
#' scme2 <- load_example_scme(dataset = "unsorted")
#' scme <- merge(scme1, scme2)
#' scme
#'
#' @export
#' @rdname SCE-merge
setMethod(
  f = "merge",
  signature = c("SingleCellMultiExperiment", "SingleCellMultiExperiment"),
  definition = function(
    x,
    y,
    experiments = NULL,
    global.reducedDims = NULL,
    assays = NULL,
    reducedDims = NULL,
    altExps = NULL,
    altExps.assays = NULL,
    label = NULL,
    add.idx = NULL,
    idx.collapse = "_",
    verbose = TRUE,
    ...
  ) {
    return(merge(
      x = x,
      y = list(y),
      experiments = experiments,
      global.reducedDims = global.reducedDims,
      assays = assays,
      reducedDims = reducedDims,
      altExps = altExps,
      altExps.assays = altExps.assays,
      label = label,
      add.idx = add.idx,
      idx.collapse = idx.collapse,
      verbose = verbose,
      ...
    ))
  }
)

#' @importClassesFrom S4Vectors list_OR_List
#' @export
#' @rdname SCE-merge
setMethod(
  f = "merge",
  signature = c("SingleCellMultiExperiment", "list_OR_List"),
  definition = function(
    x,
    y,
    experiments = NULL,
    global.reducedDims = NULL,
    assays = NULL,
    reducedDims = NULL,
    altExps = NULL,
    altExps.assays = NULL,
    label = NULL,
    add.idx = NULL,
    idx.collapse = "_",
    verbose = TRUE,
    ...
  ) {
    check_inherits_for_vec(
      x = y,
      class = "SingleCellMultiExperiment",
      name = "y"
    )
    SCMEs <- c(list(x), as.list(x = y))
    return(.merge_SCMEs(
      SCMEs = SCMEs,
      experiments = experiments,
      global.reducedDims = global.reducedDims,
      assays = assays,
      reducedDims = reducedDims,
      altExps = altExps,
      altExps.assays = altExps.assays,
      label = label,
      add.idx = add.idx,
      idx.collapse = idx.collapse,
      verbose = verbose,
      ...
    ))
  }
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Checking functions ##########################################################

.check_merge_by <- function(merge.by = NULL) {
  if (length(merge.by) == 0) {
    return(merge.by)
  }
  if (length(merge.by) > 1) {
    warning("Only take the first 'merge.by': ", merge.by[1], immediate. = TRUE)
    merge.by <- merge.by[1]
  }
  .merge_by <- c("same", "unique", "first", "only")
  if (merge.by %in% .merge_by) {
    return(merge.by)
  }
  .merge_by <- paste(.merge_by, collapse = ", ")
  stop("'merge.by' should only be one of the following: \n  ", .merge_by)
}

.check_same <- function(list) {
  vapply(list, identical, FUN.VALUE = logical(1L), y = list[[1]])
}

.check_null <- function(list) {
  vapply(list, is.null, FUN.VALUE = logical(1L))
}

## Combine List or list ########################################################
.concat_lists_helper <- function(lists, merge.by) {
  merge.by <- .check_merge_by(merge.by)
  new.list <- list()
  if (length(merge.by) == 0) {
    return(new.list)
  }
  all.nms <- lapply(lists, names)
  if (merge.by == "same") {
    # Keep elements that are the same in each of the lists.
    use.nms <- Reduce(f = intersect, x = all.nms)
    for (i in use.nms) {
      tmp <- lapply(lists, function(l) l[[i]])
      check.same <- .check_same(tmp)
      if (all(check.same)) {
        new.list[[i]] <- tmp[[1]]
      }
    }
    return(new.list)
  }
  if (merge.by == "unique") {
    # Keep elements for which there is only one possible value.
    use.nms <- all.nms %>%
      unlist() %>%
      unique()
    for (i in use.nms) {
      tmp <- lapply(
        X = lists,
        FUN = function(l) tryCatch(
          expr = l[[i]],
          error = function(e) return(NULL)
        )
      )
      check.null <- .check_null(list = tmp)
      tmp <- tmp[!check.null]
      check.same <- .check_same(list = tmp)
      if (all(check.same)) {
        new.list[[i]] <- tmp[[1]]
      }
    }
    return(new.list)
  }
  if (merge.by == "first") {
    # The first element seen at each from each position.
    use.nms <- unique(x = unlist(x = all.nms))
    for (i in seq_len(length.out = length(x = lists))) {
      use.nms <- setdiff(x = use.nms, y = names(x = new.list))
      if (length(x = use.nms) > 0 & any(use.nms %in% names(x = lists[[i]]))) {
        tmp.nms <- intersect(x = use.nms, y = names(x = lists[[i]]))
        new.list[tmp.nms] <- lists[[i]][tmp.nms]
      }
    }
    return(new.list)
  }
  if (merge.by == "only") {
    # Elements that show up in only one of the objects.
    use.nms <- names(x = which(x = table(unlist(x = all.nms)) == 1))
    for (i in seq_len(length.out = length(x = lists))) {
      if (any(use.nms %in% names(x = lists[[i]]))) {
        tmp.nms <- intersect(x = use.nms, y = names(x = lists[[i]]))
        new.list[[tmp.nms]] <- lists[[i]][[tmp.nms]]
      }
    }
    return(new.list)
  }
  return(new.list)
}

## Combine DataFrames or data.frames ###########################################

#' @importFrom easy.utils fastIntersect
.concatDFs_helper <- function(df.list, new.df, merge.by) {
  merge.by <- .check_merge_by(merge.by = merge.by)
  use.rows <- rownames(new.df)
  all.cols <- lapply(df.list, colnames)
  if (is.null(merge.by)) {
    return(new.df)
  }
  if (merge.by == "same") {
    # Keep columns that are the same in each of the objects.
    use.cols <- Reduce(intersect, all.cols)
    for (j in use.cols) {
      tmp <- lapply(df.list, function(df) {
        i <- fastIntersect(use.rows, rownames(df))
        df[i, j, drop = FALSE]
      })
      check.same <- .check_same(tmp)
      if (all(check.same)) {
        new.df[rownames(tmp[[1]]), j] <- tmp[[1]]
      }
    }
    return(new.df)
  }
  use.cols <- all.cols %>%
    unlist() %>%
    unique()
  if (merge.by == "unique") {
    # Keep columns for which there is only one possible value.
    use.cols <- all.cols %>%
      unlist() %>%
      unique()
    for (j in use.cols) {
      tmp <- lapply(df.list, function(df) {
        i <- fastIntersect(use.rows, rownames(df))
        tryCatch(
          expr = df[i, j, drop = FALSE],
          error = function(e) return(NULL)
        )
      })
      check.null <- .check_null(list = tmp)
      tmp <- tmp[!check.null]
      check.same <- .check_same(list = tmp)
      if (all(check.same)) {
        new.df[, j] <- NA
        new.df[rownames(tmp[[1]]), j] <- tmp[[1]][, j]
      }
    }
    return(new.df)
  }
  if (merge.by == "first") {
    # The first column seen at each from each position.
    if (length(use.cols) == 0) {
      return(new.df)
    }
    new.df[, use.cols] <- NA
    for (i in rev(seq_along(df.list))) {
      tmp.rows <- fastIntersect(use.rows, rownames(df.list[[i]]))
      tmp.cols <- intersect(use.cols, colnames(df.list[[i]]))
      new.df[tmp.rows, tmp.cols] <- df.list[[i]][tmp.rows, tmp.cols]
    }
    return(new.df)
  }
  if (merge.by == "only") {
    # Columns that show up in only one of the objects.
    all.cols <- unlist(all.cols)
    use.cols <- which(table(all.cols) == 1) %>%
      names()
    if (length(use.cols) == 0) {
      return(new.df)
    }
    new.df[, use.cols] <- NA
    for (i in seq_along(df.list)) {
      if (any(use.cols %in% colnames(df.list[[i]]))) {
        tmp.rows <- fastIntersect(use.rows, rownames(df.list[[i]]))
        tmp.cols <- intersect(use.cols, colnames(df.list[[i]]))
        new.df[tmp.rows, tmp.cols] <- df.list[[i]][tmp.rows, tmp.cols]
      }
    }
    new.df <- new.df[, use.cols, drop = FALSE]
    return(new.df)
  }
  return(new.df)
}

#' @importFrom easy.utils fastIntersect
.concat_DFs <- function(
    df.list,
    r.join = c("inner", "outer"),
    merge.by = NULL
) {
  r.join <- match.arg(r.join)
  all.rows <- lapply(df.list, rownames)
  all.cols <- lapply(df.list, colnames)
  if (r.join == "inner") {
    use.rows <- Reduce(fastIntersect, all.rows)
    new.df <- df.list[[1]][use.rows, character(), drop = FALSE]
  } else if (r.join == "outer") {
    all.rows <- all.rows %>%
      unlist()
    new.df <- lapply(df.list, function(df) df[, character(), drop = FALSE]) %>%
      Reduce(f = rbind)
    new.df <- new.df[!duplicated(all.rows), , drop = FALSE]
    rownames(new.df) <- unique(all.rows)
  }
  .concatDFs_helper(
    df.list = df.list,
    new.df = new.df,
    merge.by = merge.by
  )
}

#' @importFrom easy.utils fastIntersect
.rbind_DFs <- function(df.list, join = c("inner", "outer")) {
  join <- match.arg(join)
  all.cols <- lapply(X = df.list, FUN = colnames)
  if (join == "inner") {
    use.cols <- Reduce(f = fastIntersect, x = all.cols)
    for (i in seq_along(df.list)) {
      df.list[[i]] <- df.list[[i]][, use.cols, drop = FALSE]
    }
    new.df <- Reduce(rbind, df.list)
    return(new.df)
  }
  use.cols <- all.cols %>%
    unlist() %>%
    unique()
  for (i in seq_along(df.list)) {
    na.cols <- setdiff(use.cols, colnames(df.list[[i]]))
    if (length(na.cols) > 0) {
      df.list[[i]][, na.cols] <- NA
    }
  }
  new.df <- Reduce(rbind, df.list)
  return(new.df)
}
