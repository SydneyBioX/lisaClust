#' Use k-means clustering to cluster local indicators of spatial association. For other clustering use lisa.
#'
#' @param cells A SingleCellExperiment, SpatialExperiment or data frame that contains at least the
#' variables x and y, giving the  coordinates of each cell, imageID and cellType.
#' @param k The number of regions to cluster.
#' @param r A vector of the radii that the measures of association should be calculated.
#' @param imageID The column which contains image identifiers.
#' @param cellType The column which contains the cell types.
#' @param spatialCoords The columns which contain the x and y spatial coordinates.
#' @param regionName The output column for the lisaClust regions.
#' @param cores Number of cores to use for parallel processing, or a BiocParallel 
#' param object (such as MulticoreParam,
#' SnowParam or SerialParam).
#' @param window Should the window around the regions be 'square', 'convex' or 'concave'.
#' @param window.length A tuning parameter for controlling the level of concavity
#' when estimating concave windows.
#' @param whichParallel Should the function use parallization on the imageID or
#' the cellType.
#' @param sigma A numeric variable used for scaling when filting inhomogeneous L-curves.
#' @param lisaFunc Either "K" or "L" curve.
#' @param minLambda  Minimum value for density for scaling when fitting inhomogeneous L-curves.
#' @param BPPARAM \{DEPRECATED\} A BiocParallel param object.
#' @param Rs A vector of the radii that the measures of association should be calculated.
#' @param cellWeights Optional cell weights, a vector or a soft membership matrix in the order of `cells`; see
#'   [lisa()].
#'
#' @return `cells`, in the same order, with each cell's region in the column `regionName`.
#'
#' @examples
#' ## Generate toy data: two images, two cell types that sit in separate bands
#' set.seed(51773)
#' x <- round(c(
#'   runif(200), runif(200) + 1, runif(200) + 2, runif(200) + 3,
#'   runif(200) + 3, runif(200) + 2, runif(200) + 1, runif(200)
#' ), 4) * 100
#' y <- round(c(
#'   runif(200), runif(200) + 1, runif(200) + 2, runif(200) + 3,
#'   runif(200), runif(200) + 1, runif(200) + 2, runif(200) + 3
#' ), 4) * 100
#' cellType <- factor(paste("c", rep(rep(c(1:2), rep(200, 2)), 4), sep = ""))
#' imageID <- rep(c("s1", "s2"), c(800, 800))
#' cells <- data.frame(x, y, cellType, imageID)
#'
#' # Cluster the cells into two regions
#' cells <- lisaClust(cells, k = 2)
#' table(cells$region, cells$cellType)
#'
#' @export
#' @rdname lisaClust
#' @importFrom SummarizedExperiment colData
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom BiocParallel MulticoreParam
#' @importFrom stats kmeans
lisaClust <-
  function(cells,
           k = 2,
           r = NULL,
           imageID = "imageID",
           cellType = "cellType",
           spatialCoords = c("x", "y"),
           regionName = "region",
           cores = 1,
           window = "convex",
           window.length = NULL,
           whichParallel = "imageID",
           sigma = NULL,
           lisaFunc = "K",
           minLambda = 0.05,
           BPPARAM = NULL,
           Rs = r,
           cellWeights = NULL) {
    
    user_args = as.list(match.call())[-1]
    
    user_vals = lapply(names(user_args), function(arg) {
      if (arg %in% c("BPPARAM", "cores", "Rs")) {
        eval(user_args[[arg]])
      } 
    })
    
    names(user_vals) = names(user_args)
    argumentChecks("lisaClust", user_vals)
    
    
  
    if (!is(cells, "SummarizedExperiment") && !is(cells, "data.frame")) {
      stop("Data must be in the form of a SingleCellExperiment, SpatialExperiment, or data frame.")
    }
    
    cd <- .formatCells(cells, imageID, cellType, spatialCoords, allowNA = FALSE)
    # cd is ordered by image and lisa() keeps that order; .row records where each cell came from
    rowOf <- cd$.row
    cd <- cd[, c("imageID", "cellType", "x", "y", "cellID", "imageCellID")]
    if (!is.null(cellWeights) && NROW(cellWeights) != length(rowOf)) {
      stop("`cellWeights` has ", NROW(cellWeights), " ", if (is.null(dim(cellWeights))) "values" else "rows",
           " but there are ", length(rowOf), " cells (one per cell, in the order of `cells`).", call. = FALSE)
    }
    
    lisaCurves <- lisa(cd,
                       r = Rs,
                       cores = cores,
                       window = window,
                       window.length = window.length,
                       whichParallel = whichParallel,
                       sigma = sigma,
                       lisaFunc = lisaFunc,
                       minLambda = minLambda,
                       cellWeights = if (is.null(cellWeights)) NULL
                         else if (is.null(dim(cellWeights))) cellWeights[rowOf]
                         else cellWeights[rowOf, , drop = FALSE]
    )
    kM <- kmeans(lisaCurves, k)
    regions <- character(length(rowOf))
    regions[rowOf] <- paste("region", kM$cluster, sep = "_")
    
    if (is(cells, "SummarizedExperiment")) {
      SummarizedExperiment::colData(cells)[regionName] <- regions
    } else {
      cells[regionName] <- regions
    }
    
    cells
  }
