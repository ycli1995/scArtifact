
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Checking functions ###########################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom rlang inherits_any
check_inherits_for_vec <- function(x, class, name) {
  check.out <- vapply(x, inherits_any, FUN.VALUE = logical(1L), class = class)
  if (all(check.out)) {
    return(invisible(NULL))
  }
  idx <- seq_along(x)[!check.out]
  stop(
    "The following element(s) in '", name, "'",
    " doesn't inherit from ", paste(class, collapse = ", "), ":\n  ",
    paste(idx, collapse = ", ")
  )
}

#' @importFrom rlang inherits_any
check_inherits_for_func <- function(value, classes, func, name) {
  check.class <- inherits_any(value, class = classes)
  if (check.class) {
    return(invisible(NULL))
  }
  stop(
    " Unsupported class: ", class(value), ". \n",
    " '", name, "' in '", func, "' should be the following classes: \n",
    paste(classes, collapse = ", "),
    call. = FALSE
  )
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Printing functions ###########################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






