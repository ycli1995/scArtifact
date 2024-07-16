
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Constant #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Field keys ##################################################################
.int_cd <- "int_colData"

.int_emd <- "int_elementMetadata"

.int_md <- "int_metadata"

.colp_key <- "colPairs"

.rowp_key <- "rowPairs"

.alt_key <- "altExps"

.red_key <- "reducedDims"

.frag_key <- "fragments"

.annot_key <- "annotations"

.sinfo_key <- "seqinfo"

.path_key <- "filepath"

## Supported classes ###########################################################

.support_frgas <- c("IterableFragments", "Fragment")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Global variables #############################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.BPCells.envs <- new.env(parent = emptyenv())

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# redefine functions in other packages #########################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.redefine_internal <- function() {
  environment(.matchReorderSub) <- asNamespace("MultiAssayExperiment")
  assignInNamespace(
    ".matchReorderSub",
    .matchReorderSub,
    ns = "MultiAssayExperiment"
  )
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package ######################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @import BPCells
#' @importFrom stats setNames
#' @importFrom utils head
#' @importFrom tools file_path_as_absolute
#' @importFrom rlang is_false is_true
#' @importFrom easy.utils fastIntersect verboseMsg
#' @importMethodsFrom Matrix t
#' @useDynLib scArtifact
#' @keywords internal
"_PACKAGE"

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Hooks ########################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.onLoad <- function(libname, pkgname) {
  .redefine_internal()
  .init_frags_filemap()
  return(invisible(x = NULL))
}
