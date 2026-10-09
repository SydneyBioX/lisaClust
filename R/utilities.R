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
      # a number of cores is turned into a BiocParallel param by .bpparam()
    } else if (is(user_vals$cores, "MulticoreParam") || is(user_vals$cores, "SerialParam")) {
      assign("cores", user_vals$cores, sys.frame(sys.parent(1)))
    } else {
      stop("'cores'  must be either a numeric value, or a MulticoreParam or SerialParam object.\n")
    }
  }
}


# A BiocParallel param from a number of cores or a param.
#' @importFrom BiocParallel MulticoreParam SerialParam
.bpparam <- function(cores) {
  if (is(cores, "BiocParallelParam")) return(cores)
  if (is.null(cores) || cores <= 1) BiocParallel::SerialParam() else BiocParallel::MulticoreParam(workers = cores)
}

#' Put cells in a canonical data frame
#'
#' The columns imageID, cellType, x and y, and cellID and imageCellID when the data have none. imageID and
#' cellType become factors with levels in order of first appearance, and the cells are ordered by image.
#' Other columns are kept.
#'
#' @importFrom methods is
#' @importFrom S4Vectors as.data.frame
#' @noRd
.formatCells <- function(cells, imageID, cellType, spatialCoords) {
  if (is(cells, "SpatialExperiment")) {
    cd <- as.data.frame(SummarizedExperiment::colData(cells))
    sc <- data.frame(SpatialExperiment::spatialCoords(cells))
    # a SpatialExperiment's own coordinates are used unless spatialCoords names colData columns
    if (!all(spatialCoords %in% c(colnames(cd), colnames(sc)))) spatialCoords <- colnames(sc)[1:2]
    cd <- cd[, setdiff(colnames(cd), c("x", "y", colnames(sc))), drop = FALSE]
    cells <- cbind(cd, sc)
  } else if (is(cells, "SummarizedExperiment")) {
    cells <- as.data.frame(SummarizedExperiment::colData(cells))
  } else if (!is.data.frame(cells)) {
    stop("Data must be in the form of a SingleCellExperiment, SpatialExperiment, or data frame.")
  }
  cells <- as.data.frame(cells)
  for (col in c(imageID, cellType, spatialCoords)) {
    if (!col %in% colnames(cells)) stop(paste0("'", col, "' column not found in data"))
  }
  needed <- data.frame(
    imageID = cells[[imageID]], cellType = cells[[cellType]],
    x = cells[[spatialCoords[1]]], y = cells[[spatialCoords[2]]], stringsAsFactors = FALSE
  )
  cells <- cbind(cells[, setdiff(colnames(cells), c("imageID", "cellType", "x", "y")), drop = FALSE], needed)
  if (!is.factor(cells$imageID)) cells$imageID <- factor(cells$imageID, levels = unique(cells$imageID))
  if (!is.factor(cells$cellType)) cells$cellType <- factor(cells$cellType, levels = unique(cells$cellType))
  cells$.row <- seq_len(nrow(cells))
  cells <- cells[order(cells$imageID), , drop = FALSE]
  if (is.null(cells$cellID)) cells$cellID <- paste0("cell_", seq_len(nrow(cells)))
  if (is.null(cells$imageCellID)) {
    within <- stats::ave(seq_len(nrow(cells)), cells$imageID, FUN = seq_along)
    cells$imageCellID <- paste0(cells$imageID, "_", within)
  }
  rownames(cells) <- NULL
  cells
}

# The cells of each image: imageID, cellID, imageCellID, x, y and cellType, split by image.
#' @importFrom S4Vectors DataFrame split
.cellsByImage <- function(cells) {
  cells <- cells[, c("imageID", "cellID", "imageCellID", "x", "y", "cellType"), drop = FALSE]
  S4Vectors::split(S4Vectors::DataFrame(cells), cells$imageID)
}
