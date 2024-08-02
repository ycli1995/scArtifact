#' @importFrom methods new setClass setClassUnion
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Sparse matrix ################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importClassesFrom Matrix dgCMatrix dgRMatrix
setClassUnion(name = "CSC_OR_CSR_Matrix", members = c("dgCMatrix", "dgRMatrix"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ChromExperiment ##############################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @exportClass ChromExperiment
#' @rdname ChromExperiment
setClass(Class = "ChromExperiment", contains = "SingleCellExperiment")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SingleCellMultiExperiment ####################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @exportClass SingleCellMultiExperiment
#' @rdname SingleCellMultiExperiment
setClass(
  Class = "SingleCellMultiExperiment",
  contains = "MultiAssayExperiment",
  slots = list(
    int_SCE = "SingleCellExperiment",
    defaultExp = "character"
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Seurat extension classes #####################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# @name ChromatinAssay5-class
# @rdname ChromatinAssay5-class
# @importClassesFrom SeuratObject Assay5
# @importClassesFrom GenomicRanges GRanges
# @exportClass ChromatinAssay5
# @concept assay
# setClass(
#   Class = "ChromatinAssay5",
#   contains = "Assay5",
#   slots = list(
#     "ranges" = "GRanges",
#     "motifs" = "ANY",
#     "fragments" = "list",
#     "seqinfo" = "ANY",
#     "annotation" = "ANY",
#     "bias" = "ANY",
#     "positionEnrichment" = "list",
#     "links" = "GRanges"
#   )
# )

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Others #######################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setClassUnion(
  name = "Character_OR_Numeric",
  members = c("character", "numeric")
)
