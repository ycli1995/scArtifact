#' @importFrom S4Arrays is_sparse
NULL

setMethod(
  f = "is_sparse",
  signature = "IterableMatrix",
  definition = function(x) TRUE
)

setMethod(
  f = "is_sparse",
  signature = "RenameDims",
  definition = function(x) is_sparse(x@matrix)
)

setMethod(
  f = "is_sparse",
  signature = "MatrixSubset",
  definition = function(x) is_sparse(x@matrix)
)

setMethod(
  f = "is_sparse",
  signature = "ConvertMatrixType",
  definition = function(x) is_sparse(x@matrix)
)

setMethod(
  f = "is_sparse",
  signature = "MatrixMultiply",
  definition = function(x) FALSE
)

setMethod(
  f = "is_sparse",
  signature = "MatrixMask",
  definition = function(x) is_sparse(x@matrix)
)

setMethod(
  f = "is_sparse",
  signature = "MatrixRankTransform",
  definition = function(x) is_sparse(x@matrix)
)

setMethod(
  f = "is_sparse",
  signature = "RowBindMatrices",
  definition = function(x) .is_sparse_bpcells_bind_matrices(x)
)

setMethod(
  f = "is_sparse",
  signature = "ColBindMatrices",
  definition = function(x) .is_sparse_bpcells_bind_matrices(x)
)

setMethod(
  f = "is_sparse",
  signature = "ColBindMatrices",
  definition = function(x) .is_sparse_bpcells_bind_matrices(x)
)

setMethod(
  f = "is_sparse",
  signature = "TransformedMatrix",
  definition = function(x) is_sparse(x@matrix)
)

setMethod(
  f = "is_sparse",
  signature = "SCTransformPearson",
  definition = function(x) FALSE
)

setMethod(
  f = "is_sparse",
  signature = "SCTransformPearsonTranspose",
  definition = function(x) FALSE
)

setMethod(
  f = "is_sparse",
  signature = "SCTransformPearsonSlow",
  definition = function(x) FALSE
)

setMethod(
  f = "is_sparse",
  signature = "SCTransformPearsonTransposeSlow",
  definition = function(x) FALSE
)

setMethod(
  f = "is_sparse",
  signature = "TransformScaleShift",
  definition = function(x) {
    !any(x@active_transforms[, "shift"])
  }
)

setMethod(
  f = "is_sparse",
  signature = "TransformLinearResidual",
  definition = function(x) FALSE
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.is_sparse_bpcells_bind_matrices <- function(x) {
  res <- vapply(x@matrix_list, FUN = is_sparse, FUN.VALUE = logical(1L))
  sum(res) / length(res) > 0.5
}
