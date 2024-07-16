#' @importFrom methods coerce setAs as
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# setAs ########################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## SelfHits and sparse matrix ##################################################

#' @importFrom Matrix sparseMatrix
#' @importFrom S4Vectors from nLnode nRnode SelfHits to values values<-
#' @importClassesFrom Matrix dgCMatrix dgTMatrix
#' @importClassesFrom SeuratObject Graph Neighbor
#' @importClassesFrom S4Vectors SelfHits
NULL

setAs(
  from = "SelfHits",
  to = "dgCMatrix",
  def = function(from) {
    sparseMatrix(
      i = from(from) - 1L,
      j = to(from) - 1L,
      x = values(from)[, 1],
      dims = c(nLnode(from), nRnode(from)),
      repr = "C",
      index1 = FALSE
    )
  }
)

setAs(
  from = "SelfHits",
  to = "dgTMatrix",
  def = function(from) {
    sparseMatrix(
      i = from(from) - 1L,
      j = to(from) - 1L,
      x = values(from)[, 1],
      dims = c(nLnode(from), nRnode(from)),
      repr = "T",
      index1 = FALSE
    )
  }
)

setAs(
  from = "dgTMatrix",
  to = "SelfHits",
  def = function(from) {
    if (!identical(nrow(from), ncol(from))) {
      stop(
        "The input matrix must have identical row numbers and column numbers ",
        "to be coerced to 'SelfHits'."
      )
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
  from = "dgCMatrix",
  to = "SelfHits",
  def = function(from) {
    from %>%
      as("dgTMatrix") %>%
      as("SelfHits")
  }
)

setAs(
  from = "Neighbor",
  to = "dgCMatrix",
  def = function(from) nn2Sparse.Neighbor(from, repr = "C")
)

setAs(
  from = "Neighbor",
  to = "dgTMatrix",
  def = function(from) nn2Sparse.Neighbor(from, repr = "T")
)

setAs(
  from = "dgCMatrix",
  to = "Neighbor",
  def = function(from) .sparse2nn(from)
)

setAs(
  from = "dgTMatrix",
  to = "Neighbor",
  def = function(from) .sparse2nn(from)
)

setAs(
  from = "Neighbor",
  to = "SelfHits",
  def = function(from) {
    nn2Sparse.Neighbor(from, repr = "T") %>%
      as("SelfHits")
  }
)

setAs(
  from = "SelfHits",
  to = "Neighbor",
  def = function(from) {
    from %>%
      as("dgCMatrix") %>%
      .sparse2nn()
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

## matrix and LinearEmbeddingMatrix ############################################

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
# Helpers ######################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param repr One of "C" and "T", specifying the representation of the sparse
#' matrix result
#'
#' @importFrom Matrix sparseMatrix drop0
#'
#' @rdname sparse-NN
#' @export
#' @method nn2Sparse default
nn2Sparse.default <- function(object, repr = c("C", "T"), ...) {
  repr <- match.arg(repr)
  if (any(!c("idx", "dist") %in% names(repr))) {
    stop("'nn2Sparse()' needs a list containing 'idx' and 'dist'")
  }
  n.obs <- nrow(object[['idx']])
  n.neighbors <- ncol(object[['idx']])
  dist <- rcpp_get_sparse_dist(
    knn_index = object[['idx']],
    knn_dist = object[['dist']],
    n_obs = n.obs,
    n_neighbors = n.neighbors
  )
  dist <- sparseMatrix(
    i = dist$i,
    j = dist$j,
    x = dist$j,
    repr = repr,
    index1 = FALSE
  )
  dist
}

#' @importClassesFrom SeuratObject Neighbor
#'
#' @rdname sparse-NN
#' @export
#' @method nn2Sparse Neighbor
nn2Sparse.Neighbor <- function(object, repr = c("C", "T"), ...) {
  repr <- match.arg(repr)
  nn <- list(
    idx = slot(object, name = "nn.idx"),
    dist = slot(object, name = "nn.dist")
  )
  nn2Sparse(nn, repr = repr, ...)
}

#' @importFrom Matrix isSymmetric
#' @importClassesFrom Matrix dgCMatrix
#'
#' @rdname sparse-NN
#' @export
#' @method sparse2NN dgCMatrix
sparse2NN.dgCMatrix <- function(object, ...) {
  if (!isSymmetric(from)) {
    stop("Can only retrive NN indeces and distances from symmetric matrix")
  }
  indptr <- from@p
  indices <- from@i
  x <- from@x
  ncol <- indptr[2] - indptr[1]
  nrow <- ncol(object)
  nn.idx <- matrix(0L, nrow = nrow, ncol = ncol)
  nn.dist <- matrix(0, nrow = nrow, ncol = ncol)
  for (i in seq_len(nrow)) {
    idx <- ((i - 1) * ncol + 1):(i * ncol)
    ord <- order(x[idx])
    nn.dist[i, ] <- x[idx][ord]
    nn.idx[i, ] <- indices[idx][ord] + 1L
  }
  nn.dist <- cbind(rep(0, nrow), nn.dist)
  nn.idx <- cbind(1:nrow, nn.idx)
  return(list(idx = nn.idx, dist = nn.dist))
}

#' @importFrom Matrix isSymmetric
#' @importClassesFrom Matrix dgTMatrix dgCMatrix
#'
#' @rdname sparse-NN
#' @export
#' @method sparse2NN dgTMatrix
sparse2NN.dgTMatrix <- function(object, ...) {
  if (!isSymmetric(object)) {
    stop("Can only retrive NN indeces and distances from symmetric matrix")
  }
  object %>%
    as("dgCMatrix") %>%
    sparse2NN()
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importClassesFrom SeuratObject Neighbor
.sparse2nn <- function(from) {
  cell.names <- rownames(from)
  if (length(cell.names) == 0) {
    cell.names <- character(0L)
  }
  nn <- sparse2NN(from)
  new(
    Class = 'Neighbor',
    nn.idx = nn$idx,
    nn.dist = nn$dist,
    cell.names = cell.names
  )
}
