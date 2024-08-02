#' @importFrom methods coerce setAs as
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# setAs ########################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## SelfHits and sparse matrix ##################################################

#' @importFrom Matrix sparseMatrix
#' @importFrom S4Vectors from nLnode nRnode SelfHits to values values<-
#' @importClassesFrom Matrix CsparseMatrix dgCMatrix dgTMatrix symmetricMatrix
#' triangularMatrix TsparseMatrix
#' @importClassesFrom SeuratObject Graph Neighbor
#' @importClassesFrom S4Vectors SelfHits
NULL

setAs(
  from = "SelfHits",
  to = "CsparseMatrix",
  def = function(from) {
    sparseMatrix(
      i = from(from),
      j = to(from),
      x = values(from)[, 1],
      dims = c(nLnode(from), nRnode(from)),
      repr = "C",
      index1 = TRUE
    )
  }
)

setAs(
  from = "SelfHits",
  to = "TsparseMatrix",
  def = function(from) {
    sparseMatrix(
      i = from(from),
      j = to(from),
      x = values(from)[, 1],
      dims = c(nLnode(from), nRnode(from)),
      repr = "T",
      index1 = TRUE
    )
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
    if (inherits(from, "symmetricMatrix")) {
      from <- as(from, "generalMatrix")
    }
    if (inherits(from, "triangularMatrix")) {
      from <- as(from, "generalMatrix")
    }
    hits <- SelfHits(
      from = from@i + 1L,
      to = from@j + 1L,
      nnode = unique(dim(from))
    )
    values(hits) <- from@x
    return(hits)
  }
)

setAs(
  from = "CsparseMatrix",
  to = "SelfHits",
  def = function(from) {
    from %>%
      as("TsparseMatrix") %>%
      as("SelfHits")
  }
)

setAs(
  from = "Neighbor",
  to = "CsparseMatrix",
  def = function(from) nnToGraph.Neighbor(from, repr = "C")
)

setAs(
  from = "Neighbor",
  to = "TsparseMatrix",
  def = function(from) nnToGraph.Neighbor(from, repr = "T")
)

setAs(
  from = "CsparseMatrix",
  to = "Neighbor",
  def = function(from) .graphToNeighbor(from)
)

setAs(
  from = "TsparseMatrix",
  to = "Neighbor",
  def = function(from) .graphToNeighbor(from)
)

setAs(
  from = "Neighbor",
  to = "SelfHits",
  def = function(from) {
    nnToGraph.Neighbor(from, repr = "T") %>%
      as("SelfHits")
  }
)

setAs(
  from = "SelfHits",
  to = "Neighbor",
  def = function(from) {
    from %>%
      as("TsparseMatrix") %>%
      .graphToNeighbor()
  }
)

## LinearEmbeddingMatrix and DimReduc ##########################################

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

setAs(
  from = "matrix",
  to = "LinearEmbeddingMatrix",
  def = function(from) {
    loadings <- matrix(0, ncol = ncol(from))
    colnames(loadings) <- colnames(from)
    LinearEmbeddingMatrix(sampleFactors = from, featureLoadings = loadings)
  }
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importClassesFrom SeuratObject Neighbor
.graphToNeighbor <- function(from) {
  cell.names <- colnames(from)
  if (length(cell.names) == 0) {
    cell.names <- character(0L)
  }
  nn <- graphToNN(from)
  new(
    Class = 'Neighbor',
    nn.idx = nn$idx,
    nn.dist = nn$dist,
    cell.names = cell.names
  )
}
