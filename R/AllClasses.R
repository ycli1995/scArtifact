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
# Others #######################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setClassUnion(
  name = "Character_OR_Numeric",
  members = c("character", "numeric")
)
