
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom SummarizedExperiment colData
.check_duplicated_colnames_SCMEs <- function(
    SCMEs,
    collapse = "_",
    stop = FALSE
) {
  all.dimns <- SCMEs %>%
    lapply(function(x) rownames(colData(x))) %>%
    unlist()
  if (!any(duplicated(all.dimns))) {
    return(SCMEs)
  }
  if (stop) {
    stop("Duplicate colnames present across objects provided.")
  }
  warning(
    "Some colnames are duplicated across ", class(SCMEs[[1]]), " objects ",
    "provided. Renaming to enforce unique.",
    call. = FALSE, immediate. = TRUE
  )
  for (i in seq_along(SCMEs)) {
    new.names1 <- paste0(rownames(colData(SCMEs[[i]])), collapse, i)
    SCMEs[[i]] <- setColnames(x = SCMEs[[i]], new.names = new.names1)
    new.names2 <- lapply(
      X = colnames(SCMEs[[i]]),
      FUN = function(x) paste0(x, collapse, i)
    )
    SCMEs[[i]] <- setColnames(x = SCMEs[[i]], new.names = new.names2)
  }
  return(SCMEs)
}

#' @importFrom SummarizedExperiment colData
.prep_merge_SCMEs <- function(
    SCMEs,
    label = NULL,
    add.idx = NULL,
    idx.collapse = "_",
    verbose = TRUE
) {
  verboseMsg(">>>>>> Preparing SingleCellMultiExperiment(s) to be merged:")
  if (length(add.idx) > 0) {
    stopifnot(length(add.idx) == length(SCMEs))
    verboseMsg("Adding idx prefixes: ", paste(add.idx, collapse = ", "))
    for (i in seq_along(SCMEs)) {
      new.names1 <- paste0(
        add.idx[i],
        idx.collapse,
        rownames(colData(SCMEs[[i]]))
      )
      SCMEs[[i]] <- setColnames(x = SCMEs[[i]], new.names = new.names1)
      new.names2 <- lapply(
        X = colnames(SCMEs[[i]]),
        FUN = function(x) paste0(add.idx[i], idx.collapse, x)
      )
      SCMEs[[i]] <- setColnames(x = SCMEs[[i]], new.names = new.names2)
    }
  }
  if (!is.null(label)) {
    if (is.null(names(SCMEs))) {
      names(x = SCMEs) <- seq_along(SCMEs) %>%
        as.character()
    }
    verboseMsg(
      "Adding column '", label, "' to specify batch info: ",
      paste(names(x = SCMEs), collapse = ", ")
    )
    for (i in names(SCMEs)) {
      colData(SCMEs[[i]])[, label] <- i
    }
  }
  .check_duplicated_colnames_SCMEs(SCMEs = SCMEs)
}

#' @importFrom S4Vectors merge
#' @importFrom SummarizedExperiment assayNames assays assays<- colData
#' @importFrom SingleCellExperiment reducedDimNames reducedDims
#' @importFrom MultiAssayExperiment listToMap mapToList
.merge_SCMEs <- function(
    SCMEs,
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
  SCMEs <- .prep_merge_SCMEs(
    SCMEs = SCMEs,
    label = label,
    add.idx = add.idx,
    idx.collapse = idx.collapse,
    verbose = verbose
  )
  # colData
  verboseMsg(">>>>>> Merging colData...")
  cdata <- SCMEs %>%
    lapply(FUN = colData) %>%
    .rbind_DFs(join = "outer")

  # merge experiments
  all.exps <- SCMEs %>%
    lapply(FUN = names) %>%
    unlist() %>%
    table()
  e.use <- experiments %||% names(x = all.exps)
  e.use <- names(x = all.exps) %>%
    intersect(x = e.use)
  if (length(x = e.use) == 0) {
    stop("'experiments' not found: ", paste(experiments, collapse = "', '"))
  }
  new.exps <- list()
  new.smap <- list()
  for (e in e.use) {
    verboseMsg(">>>>>> Merging experiment: ", e)
    tmp.list <- lapply(X = SCMEs, FUN = experiment, e = e)
    if (length(x = tmp.list) == 1) {
      new.exps[[e]] <- tmp.list[[1]]
      if (any(assays %in% assayNames(x = new.exps[[e]]))) {
        a.use <- new.exps[[e]] %>%
          assayNames() %>%
          intersect(x = assays)
        assays(x = new.exps[[e]]) <- assays(x = new.exps[[e]])[a.use]
      }
      if (any(reducedDims %in% reducedDimNames(x = new.exps[[e]]))) {
        r.use <- new.exps[[e]] %>%
          reducedDimNames() %>%
          intersect(x = reducedDims)
        reducedDims(x = new.exps[[e]]) <- reducedDims(x = new.exps[[e]])[r.use]
      }
      if (any(altExps %in% altExpNames(x = new.exps[[e]]))) {
        alt.use <- new.exps[[e]] %>%
          altExpNames() %>%
          intersect(x = altExps)
        altExps(x = new.exps[[e]]) <- altExps(x = new.exps[[e]])[altExps]
      }
      next
    }
    new.exps[[e]] <- merge(
      x = tmp.list[[1]],
      y = tmp.list[2:length(x = tmp.list)],
      assays = assays,
      reducedDims = reducedDims,
      altExps = altExps,
      altExps.assays = altExps.assays,
      label = NULL,
      add.idx = NULL,
      idx.collapse = "_",
      verbose = verbose
    )
    new.smap[[e]] <- SCMEs %>%
      lapply(FUN = function(x) mapToList(dfmap = sampleMap(x = x))[[e]]) %>%
      Reduce(f = rbind)
  }
  new.smap <- listToMap(listmap = new.smap)
  cdata <- cdata[rownames(x = cdata) %in% new.smap[['primary']], , drop = FALSE]
  new.sce <- .merge_SCEs(
    SCEs = lapply(X = SCMEs, FUN = int_SCE),
    reducedDims = global.reducedDims,
    verbose = verbose
  )[, rownames(x = cdata)]
  SingleCellMultiExperiment(
    experiments = new.exps,
    colData = cdata,
    sampleMap = new.smap,
    reducedDims = reducedDims(x = new.sce)
  )
}
