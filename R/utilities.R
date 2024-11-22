#' A function to handle validity of argumemts/check for deprecated arguments
#'
#' @importFrom lifecycle deprecate_warn
#' @importFrom methods is
#' @noRd
argumentChecks = function(function_name, user_vals) {
  
  rlang::local_options(lifecycle_verbosity = "warning")  
  
  # handle deprecated arguments
  handle_deprecated = function(old_arg, new_arg, user_vals) {
    if (old_arg %in% names(user_vals) && !(new_arg %in% names(user_vals))) {
      # warning(paste0("'", old_arg, "' was deprecated in 1.18.0. Please use '", new_arg, "' instead.\n"))
      deprecate_warn("1.14.4", paste0(function_name, "(", old_arg, ")"), paste0(function_name, "(", new_arg, ")"))
      assign(new_arg, user_vals[[old_arg]], envir = sys.frame(sys.parent(1)))
    }
  }

  handle_deprecated("BPPARAM", "cores", user_vals)
  handle_deprecated("Rs", "r", user_vals)
  
  # enforce mutually exclusive arguments
  check_exclusive = function(arg_set, user_vals) {
    provided_args = intersect(arg_set, names(user_vals)) 
    if (length(provided_args) > 1) {
      stop(paste("Please specify only one of", paste(shQuote(arg_set), collapse = ", "), "\n"))
    }
  }
  
  check_exclusive(c("cores", "BPPARAM"), user_vals)
  check_exclusive(c("Rs", "r"), user_vals)
  
  # validity checks for cores/BPPARAM
  if ("BPPARAM" %in% names(user_vals)) {
    # warning("'BPPARAM' was deprecated in 1.18.0. Please use 'cores' instead.\n")
    deprecate_warn("1.14.4", paste0(function_name, "(BPPARAM)"), paste0(function_name, "(cores)"))
    
    if (is(user_vals$BPPARAM, "MulticoreParam") || is(user_vals$BPPARAM, "SerialParam")) {
      assign("cores", user_vals$BPPARAM, envir = sys.frame(sys.parent(1)))
    } else {
      stop("'BBPARAM' must be a MulticoreParam or SerialParam object.")
    }
    
  } else if ("cores" %in% names(user_vals)) {
    if (is(user_vals$cores, "numeric")) {
      assign("cores", simpleSeg:::generateBPParam(cores = user_vals$cores))
    } else if (is(user_vals$cores, "MulticoreParam") || is(user_vals$cores, "SerialParam")) {
      assign("cores", user_vals$cores, sys.frame(sys.parent(1)))
    } else {
      stop("'cores'  must be either a numeric value, or a MulticoreParam or SerialParam object.\n")
    }
  }
}
