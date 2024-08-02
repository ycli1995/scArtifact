
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S3 Methods ###################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname Obj-IO
#' @export
#' @method prepObj default
prepObj.default <- function(object, ...) {
  return(object)
}

#' @importFrom SingleCellExperiment factorData featureLoadings
#' LinearEmbeddingMatrix sampleFactors
#' @importFrom S4Vectors metadata
#' @importClassesFrom SingleCellExperiment LinearEmbeddingMatrix
#'
#' @rdname Obj-IO
#' @export
#' @method prepObj LinearEmbeddingMatrix
prepObj.LinearEmbeddingMatrix <- function(object, ...) {
  list(
    sampleFactors = sampleFactors(object) %>% as.matrix(),
    featureLoadings = featureLoadings(object) %>% as.matrix(),
    factorData = factorData(object) %>% as.data.frame(),
    metadata = metadata(object),
    samples = rownames(object),
    features = rownames(featureLoadings(object)),
    factors = colnames(object)
  )
}

#' @importClassesFrom GenomeInfoDb Seqinfo
#'
#' @rdname Obj-IO
#' @export
#' @method prepObj Seqinfo
prepObj.Seqinfo <- function(object, ...) {
  sinfo.df <- as.data.frame(object)
  sinfo.df <- sinfo.df[, colSums(is.na(sinfo.df)) == 0, drop = FALSE]
  return(sinfo.df)
}

#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#'
#' @rdname Obj-IO
#' @export
#' @method prepObj RangedSummarizedExperiment
prepObj.RangedSummarizedExperiment <- function(object, ...) {
  object <- formatObj(object)
  .prepObj_RSE(object, ...)
}

#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#'
#' @rdname Obj-IO
#' @export
#' @method prepObj RangedSummarizedExperiment
prepObj.SingleCellExperiment <- function(object, ...) {
  object <- formatObj(object)
  .prepObj_SCE(object, ...)
}

#' @rdname Obj-IO
#' @export
#' @method prepObj ChromExperiment
prepObj.ChromExperiment <- function(object, ...) {
  object <- formatObj(object)
  .prepObj_SCCE(object, ...)
}

#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#'
#' @rdname Obj-IO
#' @export
#' @method prepObj MultiAssayExperiment
prepObj.MultiAssayExperiment <- function(object, ...) {
  object <- formatObj(object)
  .prepObj_MAE(object)
}

#' @rdname Obj-IO
#' @export
#' @method prepObj SingleCellMultiExperiment
prepObj.SingleCellMultiExperiment <- function(object, ...) {
  object <- formatObj(object)
  .prepObj_SCME(object, ...)
}

#' @importFrom SeuratObject Cells Features LayerData Layers Misc
#' @importClassesFrom SeuratObject StdAssay
#'
#' @rdname Obj-IO
#' @export
#' @method prepObj StdAssay
prepObj.StdAssay <- function(object, ...) {
  layer.names <- Layers(object)
  layers <- list()
  for (i in layer.names) {
    layers[[i]] <- LayerData(object, layer = i, ...)
  }
  layers <- remove_dense_bpcells(layers)
  return(list(
    layers = layers,
    #cells = as.matrix(object@cells),
    #features = as.matrix(object@features),
    cells = Cells(object),
    features = Features(object),
    meta.data = object[[]],
    misc = Misc(object)
  ))
}

#' @importClassesFrom SeuratObject Graph
#'
#' @rdname Obj-IO
#' @export
#' @method prepObj Graph
prepObj.Graph <- function(object, ...) {
  as(object, "dgCMatrix")
}

#' @importFrom SeuratObject Cells Embeddings Features Loadings Misc
#' @importClassesFrom SeuratObject DimReduc
#'
#' @rdname Obj-IO
#' @export
#' @method prepObj DimReduc
prepObj.DimReduc <- function(object, ...) {
  out <- list(
    embeddings = Embeddings(object),
    loadings = Loadings(object),
    cells = Cells(object),
    features = Features(object),
    column_names = colnames(object),
    jackstraw = object@jackstraw,
    misc = Misc(object)
  )
  proj.loadings <- Loadings(object, projected = TRUE)
  if (!is_empty(proj.loadings)) {
    out$projected_loadings <- proj.loadings
  }
  return(out)
}

#' @importFrom SeuratObject Cells Distances Index
#' @importClassesFrom SeuratObject Neighbor
#'
#' @rdname Obj-IO
#' @export
#' @method prepObj Neighbor
prepObj.Neighbor <- function(object, ...) {
  if (length(Cells(object)) == 0) {
    object@cell.names <- rownames(Index(object))
  }
  list(
    idx = Index(object),
    dist = Distances(object),
    cells = Cells(object),
    alg = list(alg.idx = object@alg.idx, alg.info = object@alg.info)
  )
}

#' @importFrom SeuratObject AddMetaData Idents Misc
#' @importClassesFrom SeuratObject Seurat StdAssay
#'
#' @rdname Obj-IO
#' @export
#' @method prepObj Seurat
prepObj.Seurat <- function(object, ...) {
  object <- AddMetaData(object, Idents(object), col.name = "active.ident")
  all.assays <- list()
  assay.names <- Assays(object)
  for (i in assay.names) {
    all.assays[[i]] <- object@assays[[i]]
    if (!inherits(all.assays[[i]], "StdAssay")) {
      warning(
        "Assay '", i, "' does not inherit from 'StdAssay' of Seurat v5. ",
        "Converting into 'Assay5'. ",
        immediate. = TRUE, call. = FALSE
      )
      all.assays[[i]] <- as(all.assays[[i]], "Assay5")
    }
  }
  list(
    assays = all.assays,
    meta.data = object[[]],
    reductions = object@reductions,
    graphs = object@graphs,
    neighbors = object@neighbors,
    misc = Misc(object),
    commands = object@commands,
    tools = object@tools
  )
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Dense BPCells ###############################################################

remove_dense_bpcells <- function(rlist) {
  to.remove <- logical(length(rlist))
  for (i in seq_along(rlist)) {
    if (inherits(rlist[[i]], "IterableMatrix")) {
      to.remove[i] <- !is_sparse(rlist[[i]])
      if (to.remove[i]) {
        warning(
          "Skip writing [[", i, "]] element: ",
          "The ", class(rlist[[i]])[1], " is dense, which is not supported ",
          "for BPCells::write_matrix_dir() currently.",
          immediate. = TRUE, call. = FALSE
        )
      }
    }
  }
  rlist[!to.remove]
}

## Bioconductor S4 classes #####################################################

#' @importFrom SummarizedExperiment assays colData rowData rowRanges
#' @importFrom S4Vectors mcols<- metadata
#' @importClassesFrom GenomicRanges GRanges
.prepObj_RSE <- function(object, ...) {
  all.assays <- assays(object)
  all.assays <- remove_dense_bpcells(all.assays)
  for (i in seq_along(all.assays)) {
    if (inherits(all.assays[[i]], "dgCMatrix")) {
      all.assays[[i]] <- as(all.assays[[i]], "IterableMatrix")
    }
  }
  rr <- rowRanges(object)
  if (!inherits(rr, "GRanges")) {
    return(list(
      rowData = rowData(object),
      colData = colData(object),
      assays = all.assays,
      metadata = metadata(object)
    ))
  }
  mcols(rr) <- NULL
  return(list(
    rowRanges = rr,
    rowData = rowData(object),
    colData = colData(object),
    assays = all.assays,
    metadata = metadata(object)
  ))
}

#' @importFrom SingleCellExperiment altExps colPairs reducedDims rowPairs
.prepObj_SCE <- function(object, ...) {
  out <- .prepObj_RSE(object, ...)
  out[[.red_key]] <- reducedDims(object)
  out[[.colp_key]] <- colPairs(object, asSparse = FALSE)
  out[[.rowp_key]] <- rowPairs(object, asSparse = FALSE)
  for (i in seq_along(out[[.colp_key]])) {
    out[[.colp_key]][[i]] <- as(out[[.colp_key]][[i]], "CsparseMatrix")
  }
  for (i in seq_along(out[[.rowp_key]])) {
    out[[.rowp_key]][[i]] <- as(out[[.rowp_key]][[i]], "CsparseMatrix")
  }
  out[[.alt_key]] <- altExps(object)
  return(out)
}

.prepObj_SCCE <- function(object, ...) {
  out <- .prepObj_SCE(object, ...)
  out[[.frag_key]] <- fragments(object)
  out[[.annot_key]] <- annotations(object)
  out[[.sinfo_key]] <- seqinfo(object)
  return(out)
}

#' @importFrom SummarizedExperiment colData
#' @importFrom MultiAssayExperiment experiments sampleMap
#' @importFrom S4Vectors metadata
.prepObj_MAE <- function(object, ...) {
  all.exps <- experiments(object)
  for (i in seq_along(all.exps)) {
    if (inherits(all.exps[[i]], "dgCMatrix")) {
      all.exps[[i]] <- as(all.exps[[i]], "IterableMatrix")
    }
  }
  return(list(
    experiments = all.exps,
    colData = colData(object),
    metadata = metadata(object),
    sampleMap = sampleMap(object)
  ))
}

#' @importFrom SingleCellExperiment colPairs reducedDims
.prepObj_SCME <- function(object, ...) {
  out <- .prepObj_MAE(object, ...)
  out[[.red_key]] <- reducedDims(int_SCE(object))
  out[[.colp_key]] <- colPairs(int_SCE(object), asSparse = FALSE)
  for (i in seq_along(out[[.colp_key]])) {
    out[[.colp_key]][[i]] <- as(out[[.colp_key]][[i]], "CsparseMatrix")
  }
  return(out)
}
