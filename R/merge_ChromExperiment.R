
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom dplyr filter group_by mutate ungroup
.get_rows_to_merge2 <- function(all.ranges, reduced.ranges) {
  revmap <- reduced.ranges$revmap

  all.datasets <- unique(all.ranges$dataset)
  offsets <- c(0L, tabulate(all.ranges$dataset)[-length(all.datasets)]) %>%
    cumsum()

  revmap.idx <- unlist(revmap)
  map.df <- data.frame(
    dataset = all.ranges$dataset[revmap.idx],
    reduced.ranges.idx = rep(seq_along(revmap), lengths(revmap)),
    all.ranges.idx = revmap.idx
  ) %>%
    mutate(mat.idx = all.ranges.idx - offsets[dataset]) %>%
    mutate(new.rownames = as.character(reduced.ranges)[reduced.ranges.idx]) %>%
    group_by(dataset, reduced.ranges.idx) %>%
    mutate(n.merge.rows = n()) %>%
    ungroup()
  return(split(map.df, map.df$dataset))
}

#' @importFrom Seurat CreateCategoryMatrix
.add_overlapping_rows2 <- function(mat, map.df, verbose = TRUE) {
  map.rename <- map.df %>%
    filter(n.merge.rows < 2)
  mat.rename <- mat[map.rename$mat.idx, , drop = FALSE]
  rownames(mat.rename) <- map.rename$new.rownames

  verboseMsg("Merging overlapped peak data")
  map.to.megre <- map.df %>%
    filter(n.merge.rows > 1)
  mat.to.merge <- mat[map.to.megre$mat.idx, ]
  if (inherits(mat.to.merge, "IterableMatrix")) {
    mat.to.merge <- as(mat.to.merge, "dgCMatrix")
  }
  mat.to.merge <- t(mat.to.merge)
  cat.mat <- CreateCategoryMatrix(
    labels = map.to.megre$new.rownames,
    method = "aggregate",
    cells.name = colnames(mat.to.merge)
  )
  colnames(cat.mat) <- unique(map.to.megre$new.rownames)
  new.mat <- mat.to.merge %*% cat.mat
  new.mat <- t(new.mat)

  if (inherits(mat.rename, what = "IterableMatrix")) {
    new.mat <- new.mat %>%
      as("dgCMatrix") %>%
      as("IterableMatrix") %>%
      convert_matrix_type(type = matrix_type(mat))
  }
  mat <- rbind(mat.rename, new.mat)
  return(mat)
}

#' Find matrix indices corresponding to overlapping genomic ranges
#'
#' @param all.ranges Combined genomic ranges for all assays. This should be the
#' set of ranges that \code{reduce} was run on to get \code{reduced.ranges}, and
#' contain \code{dataset} info.
#' @param reduced.ranges A set of reduced genomic ranges containing the rev.map
#' information
#'
#' @noRd
.get_rows_to_merge <- function(all.ranges, reduced.ranges) {
  revmap <- as.vector(reduced.ranges$revmap)
  assay.list <- unique(all.ranges$dataset)

  # get indices of ranges that changed
  revmap.lengths <- lengths(revmap)
  which.changed <- which(revmap.lengths > 1)
  reduced.rownames <- reduced.ranges[which.changed] %>%
    as.character()

  # preallocate
  offsets <- list()
  results <- list()

  n.changed <- length(which.changed)
  mat.idx <- vector(mode = "numeric", length = n.changed * 2)
  new.rownames <- vector(mode = "character", length = n.changed * 2)
  for (i in seq_along(assay.list)) {
    all.ranges.idx <- which(all.ranges$dataset == i)
    offsets[[i]] <- min(all.ranges.idx) - 1L
    offsets[[i]][[2]] <- max(all.ranges.idx) + 1L
    results[['mat.idx']][[i]] <- mat.idx
    results[['new.rownames']][[i]] <- new.rownames
  }
  # find sets of ranges for each dataset
  counter <- vector(mode = "numeric", length = length(assay.list))
  for (i in seq_along(which.changed)) {
    idx <- which.changed[[i]]
    all.assay <- revmap[[idx]]
    for (j in seq_along(assay.list)) {
      selected <- (all.assay > offsets[[j]][1]) & (all.assay < offsets[[j]][2])
      this.assay <- all.assay[selected]
      mat.idx <- this.assay - offsets[[j]][1]
      mat.idx <- mat.idx[mat.idx < offsets[[j]][2] & mat.idx > 0]
      for (k in seq_along(mat.idx)) {
        counter[j] <- counter[j] + 1
        results$mat.idx[[j]][[counter[j]]] <- mat.idx[[k]]
        results$new.rownames[[j]][[counter[j]]] <- reduced.rownames[[i]]
      }
    }
  }
  # remove trailing extra values in each vector
  for (i in seq_along(assay.list)) {
    results$mat.idx[[i]] <- results$mat.idx[[i]][1:counter[i]]
    results$new.rownames[[i]] <- results$new.rownames[[i]][1:counter[i]]
  }
  return(results)
}

#' Merge rows of count matrices with overlapping genomic ranges
#'
#' @param mat The original matrix
#' @param mat.idx Rows to be merged
#' @param new.rownames Row names to be merged
#' @param verbose Display messages
#'
#' @noRd
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom Matrix rowSums
#' @importMethodsFrom Matrix t
.add_overlapping_rows <- function(mat, mat.idx, new.rownames, verbose = TRUE) {
  # transpose for faster access since matrix is column major
  is.bpcells <- inherits(mat, "IterableMatrix")
  mat <- t(mat)

  nrep <- rle(new.rownames)

  # allocate
  todelete <- c()
  n.new.rownames <- length(new.rownames)
  new.mat <- vector(mode = "list", length = n.new.rownames)
  new.mat.rownames <- vector(mode = "character", length = n.new.rownames)

  to.rename.n <- length(nrep$lengths)
  to.rename.idx <- vector(mode = "numeric", length = to.rename.n)
  to.rename.names <- vector(mode = "character", length = to.rename.n)

  if (verbose & to.rename.n > 1) {
    pb <- txtProgressBar(min = 1, max = to.rename.n, style = 3, file = stderr())
  }

  x <- 1  # row index for matrix
  y <- 1  # counter for list index
  idx.counter <- 0
  for (i in seq_along(nrep$lengths)) {
    n <- nrep$lengths[[i]]
    use.rowname <- nrep$values[[i]]
    use.mat.idx <- mat.idx[x:(x + n - 1)]
    if (n < 2) {
      idx.counter <- idx.counter + 1
      # no merge needed, just rename row in-place
      # store row indices and names to do the change in one step at the end
      to.rename.idx[idx.counter] <- use.mat.idx
      to.rename.names[idx.counter] <- use.rowname
    } else {
      # merge multiple rows and add to list
      new.mat[[y]] <- mat[, use.mat.idx] %>%
        rowSums()
      # mark merged row for deletion
      todelete <- c(todelete, use.mat.idx)
      # add row names
      new.mat.rownames[y] <- use.rowname
      y <- y + 1
    }
    if (verbose & to.rename.n > 1) {
      setTxtProgressBar(pb = pb, value = i)
    }
    x <- x + n
  }
  # remove extra elements in vectors
  to.rename.idx <- to.rename.idx[1:idx.counter]
  to.rename.names <- to.rename.names[1:idx.counter]
  new.mat <- new.mat[1:(y - 1)]
  new.mat.rownames <- new.mat.rownames[1:(y - 1)]

  # transpose back
  mat <- t(mat)

  # rename matrix rows that weren't merged
  rownames(mat)[to.rename.idx] <- to.rename.names

  if (y == 1) {
    # no rows were merged, can return mat
    return(mat)
  }
  dimnames <- list(new.mat.rownames, names(new.mat[[1]]))
  if (y == 2) {
    # only one element
    new.mat <- matrix(data = new.mat[[1]], nrow = 1, dimnames = dimnames)
    new.mat <- new.mat[, colnames(mat), drop = FALSE]
  } else {
    verboseMsg("\nBinding matrix rows")
    new.mat2 <- matrix(
      data = 0,
      nrow = length(new.mat),
      ncol = length(new.mat[[1]]),
      dimnames = dimnames
    )
    for (i in seq_along(new.mat)) {
      new.mat2[i, ] <- new.mat[[i]]
    }
    new.mat <- new.mat2
    # remove rows from old matrix that were merged
    tokeep <- nrow(mat) %>%
      seq_len() %>%
      setdiff(todelete)
    mat <- mat[tokeep, ]
  }
  # add new merged rows to counts
  if (inherits(x = mat, what = "IterableMatrix")) {
    new.mat <- new.mat %>%
      as(Class = "dgCMatrix") %>%
      as(Class = "IterableMatrix") %>%
      convert_matrix_type(type = matrix_type(x = mat))
  }
  mat <- rbind(mat, new.mat)
  return(mat)
}

#' @importFrom Signac StringToGRanges
#' @importFrom GenomicRanges isDisjoint reduce
.merge_peak_single_assay <- function(
    mats,
    all.ranges,
    reduced.ranges,
    all.nonoverlap,
    verbose = TRUE
) {
  # get the new rownames for the count matrix
  new.rownames <- as.character(reduced.ranges)
  if (all.nonoverlap) {
    # Directly merge, since no overlap range
    return(mergeMatrices(
      x = mats[[1]],
      y = mats[2:length(x = mats)]
    )[new.rownames, ])
  }
  # function to look up original
  # mergeinfo <- .get_rows_to_merge(
  #   all.ranges = all.ranges,
  #   reduced.ranges = reduced.ranges
  # )
  map.info <- .get_rows_to_merge2(
    all.ranges = all.ranges,
    reduced.ranges = reduced.ranges
  )
  for (i in seq_along(mats)) {
    # get rows to merge
    # mats[[i]] <- .add_overlapping_rows(
    #   mat = mats[[i]],
    #   mat.idx = mergeinfo$mat.idx[[i]],
    #   new.rownames = mergeinfo$new.rownames[[i]],
    #   verbose = verbose
    # )
    mats[[i]] <- .add_overlapping_rows2(
      mat = mats[[i]],
      map.df = map.info[[i]],
      verbose = verbose
    )
  }
  return(mergeMatrices(
    x = mats[[1]],
    y = mats[2:length(x = mats)]
  )[new.rownames, ])
}

.merge_seqinfo <- function(SCEs, verbose = TRUE) {
  # check genomes are all the same
  verboseMsg("-----> Check genomes...")
  genomes <- SCEs %>%
    lapply(FUN = function(x) unique(x = genome(x = x))) %>%
    unlist() %>%
    unique()
  if (length(x = genomes) > 1) {
    stop("Genomes do not match")
  }

  # merge seqinfo
  verboseMsg("-----> Merging seqinfo...")
  all.seqinfo <- lapply(X = SCEs, FUN = seqinfo)
  seqinfo.present <- !sapply(X = all.seqinfo, FUN = is.null)
  if (!any(seqinfo.present)) {
    return(NULL)
  }
  # need at least one non-NULL seqinfo, otherwise just set it as NULL
  all.seqinfo <- all.seqinfo[seqinfo.present]
  new.seqinfo <- all.seqinfo[[1]]
  if (length(x = all.seqinfo) == 1) {
    return(all.seqinfo[[1]])
  }
  new.seqinfo <- all.seqinfo[[1]]
  # iteratively merge seqinfo objects
  for (x in 2:length(x = all.seqinfo)) {
    new.seqinfo <- merge(x = new.seqinfo, y = all.seqinfo[[x]])
  }
  return(new.seqinfo)
}

.merge_annotations <- function(SCEs, verbose = TRUE) {
  verboseMsg("-----> Merging annotations...")
  all.annot <- lapply(X = SCEs, FUN = annotations)
  annot.present <- !sapply(X = all.annot, FUN = is.null)
  if (!any(annot.present)) {
    return(NULL)
  }
  all.annot <- all.annot[annot.present]
  if (length(x = all.annot) == 1) {
    return(all.annot[[1]])
  }
  new.annot <- all.annot[[1]]
  for (x in 2:length(x = all.annot)) {
    if (!identical(x = new.annot, y = all.annot[[x]])) {
      warning(
        "annotationss do not match, keeping annotation from ",
        "the first object only",
        call. = FALSE, immediate. = TRUE
      )
      return(new.annot)
    }
  }
  return(new.annot)
}

#' @importFrom SummarizedExperiment assay assayNames colData rowRanges
#' @importFrom SingleCellExperiment reducedDimNames
.merge_ChromSCEs <- function(
    SCEs,
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
  for (i in seq_along(SCEs)) {
    SCEs[[i]] <- .format_csce_rownames(SCEs[[i]])
  }
  SCEs <- .prep_merge_SCEs(
    SCEs = SCEs,
    label = label,
    add.idx = add.idx,
    idx.collapse = idx.collapse,
    verbose = verbose
  )
  verboseMsg("-----> Merging colData...")
  cdata <- SCEs %>%
    lapply(FUN = colData) %>%
    .rbind_DFs(join = "outer")

  verboseMsg("-----> Check overlapped rowRanges...")
  # check that all features are equal
  n.all.features <- SCEs %>%
    lapply(FUN = rownames) %>%
    unlist() %>%
    table()
  all.identical <- all(n.all.features == length(SCEs))
  all.nonoverlap <- SCEs %>%
    vapply(isDisjoint, FUN.VALUE = logical(1)) %>%
    all()
  if (!all.identical) {
    all.nonoverlap <- n.all.features[n.all.features < length(SCEs)] %>%
      names() %>%
      StringToGRanges(sep = c(":", "-")) %>%
      isDisjoint() %>%
      all(all.nonoverlap)
  }

  verboseMsg("-----> Get reduced rowRanges...")
  # Found overlaped ranges
  # First create a merged set of granges, preserving the assay of origin
  if (all.identical) {
    all.ranges <- rowRanges(SCEs[[1]])
  } else {
    all.ranges <- list()
    for (i in seq_along(SCEs)) {
      all.ranges[[i]] <- rowRanges(SCEs[[i]])
      all.ranges[[i]]$dataset <- i
    }
    all.ranges <- Reduce(f = c, x = all.ranges)
  }
  if (all.nonoverlap) {
    reduced.ranges <- all.ranges
  } else {
    # create reduced ranges, recording the indices of the merged ranges
    reduced.ranges <- reduce(all.ranges, with.revmap = TRUE)
  }

  assays <- .search_common_names(
    SCEs = SCEs,
    func = assayNames,
    names = assays,
    verbose = verbose
  )
  verboseMsg("-----> Merging peak assay...")
  new.assays <- list()
  for (i in assays) {
    verboseMsg("Merging '", i, "'")
    mats <- list()
    for (j in seq_along(SCEs)) {
      mats[[j]] <- assay(SCEs[[j]], i = i, withDimnames = TRUE)
    }
    new.assays[[i]] <- .merge_peak_single_assay(
      mats = mats,
      all.ranges = all.ranges,
      reduced.ranges = reduced.ranges,
      all.nonoverlap = all.nonoverlap,
      verbose = verbose
    )
  }
  reduced.ranges <- as.character(reduced.ranges) %>%
    StringToGRanges(sep = c(":", "-"))

  new.altExps <- .merge_sce_altExps(
    SCEs = SCEs,
    altExps = altExps,
    altExps.assays = altExps.assays,
    verbose = verbose
  )
  new.reducs <- .merge_sce_reducedDims(
    SCEs = SCEs,
    reducedDims = reducedDims,
    verbose = verbose
  )

  verboseMsg("-----> Merging fragments...")
  new.fragments <- SCEs %>%
    lapply(FUN = fragments) %>%
    Reduce(f = c)
  new.seqinfo <- .merge_seqinfo(SCEs = SCEs, verbose = verbose)
  new.annot <- .merge_annotations(SCEs = SCEs, verbose = verbose)

  new.SCE <- ChromExperiment(
    assays = new.assays,
    colData = cdata,
    reducedDims = new.reducs,
    altExps = new.altExps,
    ranges = reduced.ranges,
    genome = new.seqinfo,
    annotations = new.annot,
    fragments = new.fragments
  )
  return(new.SCE)
}
