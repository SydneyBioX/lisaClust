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
      stop("'BPPARAM' must be a MulticoreParam or SerialParam object.")
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

# The task function and elements of the current .par_lapply() call, found by forked workers in their copy of
# this environment.
.parallel_env <- new.env(parent = emptyenv())
.run_task <- function(i) .parallel_env$task(.parallel_env$X[[i]])

# lapply(X, FUN, ...) over the workers of `cores` (a number or a BiocParallel param).
# - Forked workers start once, after the task exists, so they share the session's data without being sent it,
#   and only task numbers and results pass between processes. Forking a new process for every chunk
#   instead copies much of the session into each one.
# - BLAS runs one thread in each worker during the call. OpenBLAS otherwise starts one thread per core in every
#   worker (128 on a 128-core server), which slows the workers down many times. The setting is changed in this
#   session before the workers start, because changing it inside a worker first restarts every thread.
# - Tasks are handed out in about eight chunks per worker, so that quick and slow tasks balance out.
.par_lapply <- function(X, FUN, ..., cores = 1) {
  BPPARAM <- .bpparam(cores)
  workers <- BiocParallel::bpnworkers(BPPARAM)
  if (workers <= 1L || length(X) <= 1L) return(lapply(X, FUN, ...))
  task <- FUN
  if (...length()) {
    # The task's environment holds only FUN and the extra arguments, not X, so that socket workers aren't
    # sent all of X with every chunk.
    taskEnv <- list2env(list(FUN = FUN, args = list(...)), parent = asNamespace("lisaClust"))
    task <- function(x) do.call(FUN, c(list(x), args))
    environment(task) <- taskEnv
  }

  # Only when BLAS has several threads: restoring them after forking restarts every thread, and the new
  # threads spin for a moment, so a single-threaded BLAS is left alone.
  blasThreads <- RhpcBLASctl::blas_get_num_procs()
  if (blasThreads > 1L) {
    RhpcBLASctl::blas_set_num_threads(1L)
    on.exit(RhpcBLASctl::blas_set_num_threads(blasThreads), add = TRUE)
  }
  if (BiocParallel::bpisup(BPPARAM) || !is(BPPARAM, "MulticoreParam")) {
    # Workers that are already running, or are new R sessions, are sent the task with each chunk.
    threadVars <- Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS"), unset = NA)
    Sys.setenv(OPENBLAS_NUM_THREADS = 1L, OMP_NUM_THREADS = 1L)
    on.exit(.restoreEnv(threadVars), add = TRUE)
    return(BiocParallel::bplapply(X, task, BPPARAM = BPPARAM))
  }

  assign("task", task, envir = .parallel_env)
  assign("X", X, envir = .parallel_env)
  on.exit(rm(list = c("task", "X"), envir = .parallel_env), add = TRUE)
  tasks <- BiocParallel::bptasks(BPPARAM)
  if (tasks == 0L) {
    BiocParallel::bptasks(BPPARAM) <- min(length(X), 8L * workers)
    on.exit(BiocParallel::bptasks(BPPARAM) <- tasks, add = TRUE, after = FALSE)
  }
  BiocParallel::bpstart(BPPARAM)
  on.exit(BiocParallel::bpstop(BPPARAM), add = TRUE, after = FALSE)
  out <- BiocParallel::bplapply(seq_along(X), .run_task, BPPARAM = BPPARAM)
  names(out) <- names(X)
  out
}

# Put environment variables back as they were (NA: unset).
.restoreEnv <- function(values) {
  for (name in names(values)) {
    if (is.na(values[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, as.list(values[name]))
  }
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

# The given columns of a data frame, or of the colData of a SingleCellExperiment or SpatialExperiment.
.colDataFrame <- function(cells, columns) {
  if (is(cells, "SummarizedExperiment")) cells <- as.data.frame(SummarizedExperiment::colData(cells))
  if (!is.data.frame(cells)) stop("Data must be in the form of a SingleCellExperiment, SpatialExperiment, or data frame.")
  missing <- setdiff(columns, colnames(cells))
  if (length(missing)) stop("column(s) not found in data: ", paste0("'", missing, "'", collapse = ", "))
  as.data.frame(cells)[, columns, drop = FALSE]
}
