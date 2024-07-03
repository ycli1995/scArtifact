
#' Load in data from 10X
#'
#' Read the sparse data matrix provided by 10X genomics.
#'
#' @param path Directory containing the matrix.mtx, genes.tsv (or features.tsv),
#' and barcodes.tsv files provided by 10X.
#' @param gene.column Specify which column of genes.tsv or features.tsv to use
#' for gene names; default is 2
#' @param cell.column Specify which column of barcodes.tsv to use for cell
#' names; default is 1
#' @param use.BPCells Whether or not to use \pkg{BPCells} to hold the matrix.
#'
#' @importFrom utils read.delim read.table
#' @export
read10xMtx <- function(
    path,
    gene.column = 2,
    cell.column = 1,
    use.BPCells = FALSE
) {
  path <- file_path_as_absolute(path[1])
  has_dt <- requireNamespace("data.table", quietly = TRUE) &&
    requireNamespace("R.utils", quietly = TRUE)
  if (!dir.exists(paths = path)) {
    stop("Directory provided does not exist: ", path)
  }
  barcode.loc <- file.path(path, "barcodes.tsv")
  gene.loc <- file.path(path, "genes.tsv")
  features.loc <- file.path(path, "features.tsv.gz")
  matrix.loc <- file.path(path, "matrix.mtx")
  pre_ver_3 <- file.exists(gene.loc)
  if (!pre_ver_3) {
    barcode.loc <- paste0(barcode.loc, ".gz")
    matrix.loc <- paste0(matrix.loc, ".gz")
  }
  if (!file.exists(barcode.loc)) {
    stop("Barcode file missing: ", basename(barcode.loc))
  }
  if (!pre_ver_3 && !file.exists(features.loc)) {
    stop("Gene name or features file missing: ", basename(features.loc))
  }
  if (!file.exists(matrix.loc)) {
    stop("Expression matrix file missing: ", basename(matrix.loc))
  }
  ## Read cell barcodes
  if (has_dt) {
    cell.barcodes <- data.table::fread(
      input = barcode.loc,
      header = FALSE,
      data.table = FALSE
    )
  } else {
    cell.barcodes <- read.table(
      file = barcode.loc,
      header = FALSE,
      sep = "\t",
      row.names = NULL
    )
  }
  if (ncol(cell.barcodes) > 1) {
    rownames(cell.barcodes) <- cell.barcodes[, cell.column]
  } else {
    rownames(cell.barcodes) <- readLines(con = barcode.loc)
  }
  cell.barcodes <- data.frame(row.names = rownames(cell.barcodes))

  ## Read features
  features.file <- ifelse(test = pre_ver_3, yes = gene.loc, no = features.loc)
  if (has_dt) {
    features <- data.table::fread(
      input = features.file,
      header = FALSE,
      data.table = FALSE
    )
  } else {
    features <- read.delim(
      file = features.file,
      header = FALSE,
      stringsAsFactors = FALSE
    )
  }
  if (anyNA(features[, gene.column])) {
    warning(
      "Replacing NA features with ID from the opposite column requested",
      call. = FALSE, immediate. = TRUE
    )
    na.features <- which(is.na(features[, gene.column]))
    repl.column <- ifelse(test = gene.column == 2, yes = 1, no = 2)
    features[na.features, gene.column] <- features[na.features, repl.column]
  }
  rownames(features) <- make.unique(names = features[, gene.column])
  if (ncol(features) >= 1) {
    colnames(features)[1] <- "ID"
  }
  if (ncol(features) >= 2) {
    colnames(features)[2] <- "Name"
  }
  if (ncol(features) >= 3) {
    colnames(features)[3] <- "Type"
    features <- features[, 1:3]
  }
  sce.func <- ifelse(
    test = use.BPCells,
    yes = .read_10x_mtx_bpcells,
    no = .read_10x_mtx_sce
  )
  sce <- sce.func(matrix.loc, features, cell.barcodes)
  if (!"Type" %in% colnames(features)) {
    return(sce)
  }
  data_types <- factor(features$Type)
  all_types <- levels(data_types)
  if (length(all_types) == 1) {
    return(sce)
  }
  message(
    "10X data contains more than one type: ", paste(all_types, collapse = ", "),
    "\nReturn a list containing ", class(x = sce), " of each type."
  )
  expr_name <- "Gene Expression"
  if (expr_name %in% all_types) {
    all_types <- c(expr_name, all_types[-which(all_types == expr_name)])
  }
  sce <- lapply(sce, FUN = function(x) sce[data_types == x, , drop = FALSE])
  names(sce) <- all_types
  return(sce)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom BPCells import_matrix_market
.read_10x_mtx_bpcells <- function(mtx_path, rowData, colData) {
  mat <- import_matrix_market(
    mtx_path = mtx_path,
    row_names = rownames(rowData),
    col_names = rownames(colData)
  )
  return(SingleCellExperiment(
    assays = list(counts = mat),
    rowData = as(rowData, "DataFrame"),
    colData = as(colData, "DataFrame")
  ))
}

#' @importFrom Matrix readMM
#' @importFrom S4Vectors DataFrame
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom Matrix CsparseMatrix
.read_10x_mtx_sce <- function(mat_path, rowData, colData) {
  mat <- readMM(file = mat_path)
  rownames(mat) <- rownames(rowData)
  colnames(mat) <- rownames(colData)
  mat <- as(mat, "CsparseMatrix")
  gc(verbose = FALSE)
  return(SingleCellExperiment(
    assays = list(counts = mat),
    rowData = as(rowData, "DataFrame"),
    colData = as(colData, "DataFrame")
  ))
}
