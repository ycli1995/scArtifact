
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions ####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param path Path to the directory representing an object.
#' @param class String specifying the data class of the object.
#' @param extra Named list containing extra metadata to be written.
#'
#' @importFrom jsonlite toJSON
#' @export
#' @rdname ObjFile-IO
writeObjFile <- function(path, class, extra = list()) {
  extra$class <- class
  stopifnot(anyDuplicated(names(extra)) == 0L)
  toJSON(extra, auto_unbox = TRUE, pretty = 4) %>%
    write(file = file.path(path, .objfile))
}

#' @importFrom jsonlite fromJSON
#' @export
#' @rdname ObjFile-IO
readObjFile <- function(path) {
  fromJSON(file.path(path, .objfile))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.objfile <- "OBJECT.json"
