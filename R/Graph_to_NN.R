
#' @importFrom easy.NN nnToGraph
#' @importClassesFrom SeuratObject Neighbor
#' @export
#' @method nnToGraph Neighbor
#' @concept nearest-neighbors
nnToGraph.Neighbor <- function(
    object,
    repr = c("C", "T"),
    use.weights = TRUE,
    self.loops = TRUE,
    ...
) {
  cell.names <- object@cell.names
  object <- list(
    idx = object@nn.idx,
    dist = object@nn.dist
  )
  rownames(object$idx) <- cell.names
  nnToGraph(
    object,
    repr = repr,
    use.weights = use.weights,
    self.loops = self.loops,
    ...
  )
}
