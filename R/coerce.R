#' @importFrom methods coerce setAs as
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# setAs ########################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## SelfHits and sparse matrix ##################################################

#' @importFrom Matrix sparseMatrix
#' @importFrom S4Vectors from nLnode nRnode SelfHits to values values<-
#' @importClassesFrom Matrix dgCMatrix dgTMatrix
#' @importClassesFrom S4Vectors SelfHits
NULL

setAs(
  from = "SelfHits",
  to = "CsparseMatrix",
  def = function(from) {
    return(sparseMatrix(
      i = from(from) - 1L,
      j = to(from) - 1L,
      x = values(from)[, 1],
      dims = c(nLnode(from), nRnode(from)),
      repr = "C",
      index1 = FALSE
    ))
  }
)

setAs(
  from = "SelfHits",
  to = "TsparseMatrix",
  def = function(from) {
    return(sparseMatrix(
      i = from(from) - 1L,
      j = to(from) - 1L,
      x = values(from)[, 1],
      dims = c(nLnode(from), nRnode(from)),
      repr = "T",
      index1 = FALSE
    ))
  }
)

setAs(
  from = "TsparseMatrix",
  to = "SelfHits",
  def = function(from) {
    if (!identical(nrow(from), ncol(from))) {
      stop(
        "The input matrix must have identical row numbers and column numbers ",
        "to be coerced to 'SelfHits'."
      )
    }
    hits <- SelfHits(
      from = slot(from, name = "i") + 1L,
      to = slot(from, name = "j") + 1L,
      nnode = unique(dim(from))
    )
    values(hits) <- slot(from, name = "x")
    return(hits)
  }
)

setAs(
  from = "CsparseMatrix",
  to = "SelfHits",
  def = function(from) {
    from %>%
      as(Class = "TsparseMatrix") %>%
      as(Class = "SelfHits")
  }
)

## DimReduc and LinearEmbeddingMatrix ##########################################

#' @importFrom SeuratObject CreateDimReducObject Embeddings Key Loadings
#' Loadings<- Stdev
#' @importFrom SingleCellExperiment factorData featureLoadings
#' LinearEmbeddingMatrix sampleFactors
#' @importFrom S4Vectors DataFrame metadata
#' @importClassesFrom SingleCellExperiment LinearEmbeddingMatrix
#' @importClassesFrom SeuratObject DimReduc
NULL

setAs(
  from = "DimReduc",
  to = "LinearEmbeddingMatrix",
  def = function(from) {
    key <- Key(from)
    embeddings <- Embeddings(from)
    colnames(embeddings) <- paste0(key, seq_len(ncol(embeddings)))
    feature.loadings <- Loadings(from)
    if (is_empty(feature.loadings)) {
      feature.loadings <- matrix(0, nrow = 0, ncol = ncol(embeddings))
    }
    colnames(feature.loadings) <- colnames(embeddings)
    factor.data <- DataFrame(row.names = colnames(embeddings))
    stdev <- Stdev(from)
    if (length(stdev) == nrow(factor.data)) {
      factor.data$stdev <- stdev
    }
    metadata <- list()
    projected.loadings <- Loadings(from, projected = TRUE)
    if (!is_empty(projected.loadings)) {
      colnames(projected.loadings) <- colnames(embeddings)
      metadata$feature.loadings.projected <- projected.loadings
    }
    LinearEmbeddingMatrix(
      sampleFactors = embeddings,
      featureLoadings = feature.loadings,
      factorData = factor.data,
      metadata = metadata
    )
  }
)

setAs(
  from = "LinearEmbeddingMatrix",
  to = "DimReduc",
  def = function(from) {
    factor.data <- factorData(from)
    embeddings <- sampleFactors(from)
    loadings <- featureLoadings(from)

    stdev <- numeric()
    if (length(factor.data$stdev) > 0) {
      stdev <- factor.data$stdev
    }
    out <- CreateDimReducObject(
      embeddings = embeddings,
      loadings = loadings,
      stdev = stdev
    )
    projected.loadings <- metadata$feature.loadings.projected
    if (is_empty(projected.loadings)) {
      return(out)
    }
    if (identical(dimnames(projected.loadings), dimnames(loadings))) {
      Loadings(out, projected = TRUE) <- projected.loadings
      return(out)
    }
    warning(
      "The dimension names of projected loadings do not match the original ",
      "loadings, remove metadata$feature.loadings.projected",
      call. = FALSE, immediate. = TRUE
    )
    return(out)
  }
)

## matrix and LinearEmbeddingMatrix ##########################################

setAs(
  from = "matrix",
  to = "LinearEmbeddingMatrix",
  def = function(from) {
    loadings <- matrix(0, ncol = ncol(from))
    colnames(loadings) <- colnames(from)
    LinearEmbeddingMatrix(sampleFactors = from, featureLoadings = loadings)
  }
)


