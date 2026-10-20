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
    
    if (is(user_vals$BPPARAM, "BiocParallelParam")) {
      assign("cores", user_vals$BPPARAM, envir = sys.frame(sys.parent(1)))
    } else {
      stop("'BPPARAM' must be a BiocParallel param object.")
    }
    
  } else if ("cores" %in% names(user_vals)) {
    if (is(user_vals$cores, "numeric")) {
      # a number of cores is turned into a BiocParallel param by .bpparam()
    } else if (is(user_vals$cores, "BiocParallelParam")) {
      assign("cores", user_vals$cores, sys.frame(sys.parent(1)))
    } else {
      stop("'cores' must be either a numeric value or a BiocParallel param object.\n")
    }
  }
}


# A BiocParallel param from a number of cores (forked workers, or socket workers on Windows) or a param.
#' @importFrom BiocParallel MulticoreParam SerialParam
.bpparam <- function(cores) {
  if (is(cores, "BiocParallelParam")) return(cores)
  if (is.null(cores) || cores <= 1) return(BiocParallel::SerialParam())
  if (.Platform$OS.type == "windows") BiocParallel::SnowParam(workers = cores) else BiocParallel::MulticoreParam(workers = cores)
}

# The task function and elements of the current .par_lapply() call, found by forked workers in their copy of
# this environment.
.parallel_env <- new.env(parent = emptyenv())
.run_task <- function(i) .parallel_env$task(.parallel_env$X[[i]])

## task(x) in a worker with one BLAS and OpenMP thread, for workers this session did not start itself.
.one_thread_task <- function(x, task) {
  if (isTRUE(RhpcBLASctl::blas_get_num_procs() > 1L)) RhpcBLASctl::blas_set_num_threads(1L)
  if (isTRUE(RhpcBLASctl::omp_get_max_threads() > 1L)) RhpcBLASctl::omp_set_num_threads(1L)
  task(x)
}

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
  if (isTRUE(blasThreads > 1L)) {
    RhpcBLASctl::blas_set_num_threads(1L)
    on.exit(RhpcBLASctl::blas_set_num_threads(blasThreads), add = TRUE)
  }
  # OpenMP likewise: forked workers can hang if this session has used several OpenMP threads.
  ompThreads <- RhpcBLASctl::omp_get_max_threads()
  if (isTRUE(ompThreads > 1L)) {
    RhpcBLASctl::omp_set_num_threads(1L)
    on.exit(RhpcBLASctl::omp_set_num_threads(ompThreads), add = TRUE)
  }
  # About eight chunks per worker, unless the param sets its own number of tasks.
  tasks <- BiocParallel::bptasks(BPPARAM)
  if (tasks == 0L && !BiocParallel::bpisup(BPPARAM)) {
    BiocParallel::bptasks(BPPARAM) <- min(length(X), 8L * workers)
    on.exit(BiocParallel::bptasks(BPPARAM) <- tasks, add = TRUE, after = FALSE)
  }
  if (BiocParallel::bpisup(BPPARAM) || !is(BPPARAM, "MulticoreParam")) {
    # Workers that are already running, or are new R sessions, are sent the task with each chunk. New
    # socket workers start with one BLAS and OpenMP thread; workers that were already running set
    # one thread themselves before their first task.
    threadVars <- Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS"), unset = NA)
    Sys.setenv(OPENBLAS_NUM_THREADS = 1L, OMP_NUM_THREADS = 1L)
    on.exit(.restoreEnv(threadVars), add = TRUE)
    return(BiocParallel::bplapply(X, .one_thread_task, task = task, BPPARAM = BPPARAM))
  }

  assign("task", task, envir = .parallel_env)
  assign("X", X, envir = .parallel_env)
  on.exit(rm(list = c("task", "X"), envir = .parallel_env), add = TRUE)
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
  cells <- cells[, intersect(c("imageID", "cellID", "imageCellID", "x", "y", "cellType", ".row"), colnames(cells)),
                 drop = FALSE]
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

# Soft cell type membership as triplets: the non-zero weights of a cells x labels matrix (base or sparse Matrix),
# checked against the number of cells. A vector gives each cell's weight for its own cell type.
.membership <- function(w, nCells, cellType = NULL) {
  if (is.numeric(w) && is.null(dim(w))) {
    if (length(w) != nCells)
      stop("`cellWeights` has ", length(w), " values but there are ", nCells, " cells (one per cell, in the order of ",
           "`cells`).", call. = FALSE)
    if (anyNA(w) || any(!is.finite(w)) || any(w < 0) || any(w > 1))
      stop("`cellWeights` must be between 0 and 1, with no missing values.", call. = FALSE)
    lab <- if (is.factor(cellType)) levels(droplevels(cellType)) else sort(unique(as.character(cellType)))
    j <- match(as.character(cellType), lab); keep <- w != 0 & !is.na(j)
    return(list(i = which(keep), j = j[keep], x = w[keep], labels = lab))
  }
  if (is.data.frame(w)) w <- as.matrix(w)
  if (is.null(colnames(w)) || anyDuplicated(colnames(w)) || any(colnames(w) == ""))
    stop("`cellWeights` needs unique column names (the labels).", call. = FALSE)
  if (nrow(w) != nCells)
    stop("`cellWeights` has ", nrow(w), " rows but there are ", nCells, " cells (one row per cell, in the ",
         "order of `cells`).", call. = FALSE)
  if (is.matrix(w)) {
    if (!is.numeric(w)) stop("`cellWeights` must be numeric.", call. = FALSE)
    nz <- which(w != 0, arr.ind = TRUE)
    trip <- list(i = unname(nz[, 1]), j = unname(nz[, 2]), x = w[nz])
  } else if (is(w, "sparseMatrix")) {
    s <- Matrix::summary(w)
    keep <- s$x != 0
    trip <- list(i = s$i[keep], j = s$j[keep], x = as.numeric(s$x[keep]))
  } else {
    stop("`cellWeights` must be a numeric matrix, data.frame or sparse Matrix.", call. = FALSE)
  }
  if (anyNA(trip$x) || any(!is.finite(trip$x)) || any(trip$x < 0))
    stop("`cellWeights` must be finite and non-negative, with no missing values.", call. = FALSE)
  c(trip, list(labels = colnames(w)))
}

# The membership of the cells `rows` (rows of the weight matrix), as compressed rows in that order: cell k's
# labels (1-based) are lab[ptr[k] + 1 .. ptr[k + 1]] with weights w; total is each label's summed weight.
.membershipOf <- function(memb, rows) {
  pos <- match(memb$i, rows)
  keep <- !is.na(pos)
  pos <- pos[keep]; j <- memb$j[keep]; x <- memb$x[keep]
  ord <- order(pos, j)
  total <- numeric(length(memb$labels))
  s <- rowsum(x, j)
  total[as.integer(rownames(s))] <- s[, 1]
  list(ptr = c(0L, cumsum(tabulate(pos, length(rows)))), lab = as.integer(j[ord]), w = x[ord], total = total,
       labels = memb$labels)
}
