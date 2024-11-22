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
#' MulticoreParam or SerialParam object.
#' @param window Should the window around the regions be 'square', 'convex' or 'concave'.
#' @param window.length A tuning parameter for controlling the level of concavity
#' when estimating concave windows.
#' @param whichParallel Should the function use parallization on the imageID or
#' the cellType.
#' @param sigma A numeric variable used for scaling when filting inhomogeneous L-curves.
#' @param lisaFunc Either "K" or "L" curve.
#' @param minLambda  Minimum value for density for scaling when fitting inhomogeneous L-curves.
#' @param BPPARAM \{DEPRECATED\} A BiocParalell MulticoreParam or SerialParam object.
#' @param Rs A vector of the radii that the measures of association should be calculated.
#'
#' @return A matrix of LISA curves
#'
#' @examples
#' library(spicyR)
#' library(SingleCellExperiment)
#' # Read in data
#' isletFile <- system.file("extdata", "isletCells.txt.gz", package = "spicyR")
#' cells <- read.table(isletFile, header = TRUE)
#' cellExp <- SingleCellExperiment(
#'     assay = list(intensities = t(cells[, grepl(names(cells), pattern = "Intensity_")])),
#'     colData = cells[, !grepl(names(cells), pattern = "Intensity_")]
#' )
#'
#' # Cluster cell types
#' markers <- t(assay(cellExp, "intensities"))
#' kM <- kmeans(markers, 8)
#' colData(cellExp)$cluster <- paste("cluster", kM$cluster, sep = "")
#'
#' # Generate LISA
#' cellExp <- lisaClust(cellExp,
#'     k = 2,
#'     imageID = "ImageNumber",
#'     cellType = "cluster",
#'     spatialCoords = c("Location_Center_X", "Location_Center_Y")
#' )

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
           whichParallel = "imimageID",
           sigma = NULL,
           lisaFunc = "K",
           minLambda = 0.05,
           BPPARAM = BiocParallel::SerialParam(),
           Rs = r) {
    
    user_args = as.list(match.call())[-1]
    
    user_vals = lapply(names(user_args), function(arg) {
      if (arg %in% c("BPPARAM", "cores", "Rs")) {
        eval(user_args[[arg]])
      } 
    })
    
    names(user_vals) = names(user_args)
    argumentChecks("lisaClust", user_vals)
  
    if (is(cells, "SummarizedExperiment")) {
      cols = colnames(colData(cells))
    } else if (is(cells, "data.frame")) {
      cols = colnames(cells)
    } else {
      stop("Data must be in the form of a SingleCellExperiment, SpatialExperiment, or data frame.")
    }
    
    if (!(imageID %in% cols)) {
      stop(paste0("'", imageID, "' column not found in data"))
    }
    
    if (!(cellType %in% cols)) {
      stop(paste0("'", cellType, "' column not found in data"))
    }
    
    if (methods::is(cells, "SummarizedExperiment")) {
      cd <- spicyR:::.format_data(
        cells, imageID, cellType, spatialCoords, FALSE
      )
      
      lisaCurves <- lisa(cd,
                         r = r,
                         cores = cores,
                         window = window,
                         window.length = window.length,
                         whichParallel = whichParallel,
                         sigma = sigma,
                         lisaFunc = lisaFunc,
                         minLambda = minLambda
      )
      kM <- kmeans(lisaCurves, k)
      regions <- paste("region", kM$cluster, sep = "_")
      
      SummarizedExperiment::colData(cells)[regionName] <- regions
    } else if (is(cells, "data.frame")) {
      cd <- cells
      cd <- cd[, c(cellType, imageID, spatialCoords)]
      colnames(cd) <- c("cellType", "imageID", "x", "y")
      cd$cellID <- as.character(seq_len(nrow(cd)))
      cd$imageCellID <- as.character(seq_len(nrow(cd)))
      
      lisaCurves <- lisa(cd,
                         r = r,
                         cores = cores,
                         window = window,
                         window.length = window.length,
                         whichParallel = whichParallel,
                         sigma = sigma,
                         lisaFunc = lisaFunc,
                         minLambda = minLambda
      )
      
      kM <- kmeans(lisaCurves, k)
      regions <- paste("region", kM$cluster, sep = "_")
      
      cells[regionName] <- regions
    } else {
      stop(
        "Unsupported datatype for cells: please use",
        "SingleCellExperiment or SpatialExperiment"
      )
    }
    
    cells
  }
