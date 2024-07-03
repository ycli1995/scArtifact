
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Roxygen2 calls ###############################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.doc_links <- function(nm) {
  pkg <- .pkg_map[nm]
  paste0("\\code{\\link[", pkg, "]{", nm, "}}")
}

.dot_param <- "Arguments passed to other metheds."
.val_param <- "An object of a class specified in the S4 method signature."
.vb_param <- "Print progress."

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.pkg_map <- c(
  "metadata" = "S4Vectors",

  "GRanges" = "GenomicRanges:GRanges-class",

  "HDF5Array" = "HDF5Array",

  "Seqinfo" = "GenomeInfoDb:Seqinfo-class",

  "SummarizedExperiment" = "RangedSummarizedExperiment-class",
  "assay" = "SummarizedExperiment:SummarizedExperiment-class",
  "assays" = "SummarizedExperiment:SummarizedExperiment-class",
  "rowData" = "SummarizedExperiment:SummarizedExperiment-class",
  "colData" = "SummarizedExperiment:SummarizedExperiment-class",

  "SingleCellExperiment" = "SingleCellExperiment",
  "reducedDim" = "SingleCellExperiment",
  "reducedDims" = "SingleCellExperiment",
  "reducedDimNames" = "SingleCellExperiment",
  "colPair" = "SingleCellExperiment",
  "colPairs" = "SingleCellExperiment",
  "colPairNames" = "SingleCellExperiment",
  "rowPair" = "SingleCellExperiment",
  "rowPairs" = "SingleCellExperiment",
  "rowPairNames" = "SingleCellExperiment",
  "altExp" = "SingleCellExperiment",
  "altExps" = "SingleCellExperiment",
  "altExpNames" = "SingleCellExperiment",
  "int_colData" = "SingleCellExperiment:SCE-internals",
  "int_metadata" = "SingleCellExperiment:SCE-internals",
  "int_elementMetadata" = "SingleCellExperiment:SCE-internals",

  "MultiAssayExperiment" = "MultiAssayExperiment",
  "getWithColData" = "MultiAssayExperiment",
  "experiments" = "MultiAssayExperiment:MultiAssayExperiment-methods",
  "subsetBy" = "MultiAssayExperiment",

  "IterableFragments" = "BPCells:IterableFragments-methods",
  "IterableMatrix" = "BPCells:IterableMatrix-methods",
  "FragmentsDir" = "BPCells:IterableFragments-methods",
  "cellNames" = "BPCells",
  "write_fragments_dir" = "BPCells",
  "write_fragments_hdf5" = "BPCells",
  "open_fragments_10x" = "BPCells",
  "open_fragments_hdf5" = "BPCells",
  "open_fragments_dir" = "BPCells",

  "Fragment" = "Signac",
  "Annotation" = "Signac",
  "ChromatinAssay" = "Signac:ChromatinAssay-class"
)








