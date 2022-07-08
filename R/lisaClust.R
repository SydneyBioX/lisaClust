#' Use k-means clustering to cluster local indicators of spatial association. For other clustering use lisa.
#'
#' @param cells A SegmentedCells or data frame that contains at least the 
#' variables x and y, giving the  coordinates of each cell, imageID and cellType.
#' @param k The number of regions to cluster.
#' @param Rs A vector of the radii that the measures of association should be calculated.
#' @param spatialCoords The columns which contain the x and y spatial coordinates.
#' @param cellType The column which contains the cell types.
#' @param imageID The column which contains image identifiers.
#' @param regionName The output column for the lisaClust regions.
#' @param BPPARAM A BiocParallelParam object.
#' @param window Should the window around the regions be 'square', 'convex' or 'concave'.
#' @param window.length A tuning parameter for controlling the level of concavity 
#' when estimating concave windows.
#' @param whichParallel Should the function use parallization on the imageID or 
#' the cellType.
#' @param sigma A numeric variable used for scaling when filting inhomogeneous L-curves.
#' @param lisaFunc Either "K" or "L" curve.
#' @param minLambda  Minimum value for density for scaling when fitting inhomogeneous L-curves.
#' @param fast A logical describing whether to use a fast approximation of the 
#' inhomogeneous local L-curves.
#'
#' @return A matrix of LISA curves
#'
#' @examples
#' library(spicyR)
#' # Read in data as a SegmentedCells objects
#' isletFile <- system.file("extdata","isletCells.txt.gz", package = "spicyR")
#' cells <- read.table(isletFile, header=TRUE)
#' cellExp <- SegmentedCells(cells, cellProfiler = TRUE)
#'
#' # Cluster cell types
#' markers <- cellMarks(cellExp)
#' kM <- kmeans(markers,8)
#' cellType(cellExp) <- paste('cluster',kM$cluster, sep = '')
#'
#' # Generate LISA
#' cellExp <- lisaClust(cellExp, k = 2)
#'
#' @export
#' @rdname lisaClust
#' @importFrom SummarizedExperiment colData
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom spicyR cellAnnotation<- cellAnnotation
#' @importFrom stats kmeans
#' @importFrom spicyR SegmentedCells
lisaClust <- 
    function(cells,
             k = 2,
             Rs = NULL,
             spatialCoords = c("x","y"),
             cellType = "cellType",
             imageID = "imageID",
             regionName = "region",
             BPPARAM = BiocParallel::SerialParam(),
             window = "convex",
             window.length = NULL,
             whichParallel = 'imageID',
             sigma = NULL,
             lisaFunc = "K",
             minLambda = 0.05,
             fast = TRUE) {
        
        
        
        
        if(is(cells, "SingleCellExperiment")){
            cd <- as.data.frame(SingleCellExperiment::colData(cells))
            cd <- cd[,c(cellType, imageID, spatialCoords)]
            colnames(cd) <- c("cellType", "imageID", "x", "y")
            cd$cellID <- as.character(seq_len(nrow(cd)))
            cd$imageCellID <- as.character(seq_len(nrow(cd)))
            cd <- spicyR::SegmentedCells(cd)
            lisaCurves <- lisa(cd,
                               Rs = Rs,
                               BPPARAM = BPPARAM,
                               window = window,
                               window.length = window.length,
                               whichParallel = whichParallel,
                               sigma = sigma,
                               lisaFunc = lisaFunc,
                               minLambda = minLambda,
                               fast = fast)
            
            kM <- kmeans(lisaCurves,k)
            regions <- paste('region',kM$cluster,sep = '_')
            SummarizedExperiment::colData(cells)[regionName] <- regions
        }
        
        if(is(cells, "SpatialExperiment")){
            cd <- cbind(as.data.frame(SingleCellExperiment::colData(cells)),as.data.frame(SpatialExperiment::spatialCoords(cells)))
            cd <- cd[,c(cellType, imageID, spatialCoords)]
            colnames(cd) <- c("cellType", "imageID", "x", "y")
            cd$cellID <- as.character(seq_len(nrow(cd)))
            cd$imageCellID <- as.character(seq_len(nrow(cd)))
            cd <- SegmentedCells(cd)
            lisaCurves <- lisa(cd,
                               Rs = Rs,
                               BPPARAM = BPPARAM,
                               window = window,
                               window.length = window.length,
                               whichParallel = whichParallel,
                               sigma = sigma,
                               lisaFunc = lisaFunc,
                               minLambda = minLambda,
                               fast = fast)
            
            kM <- kmeans(lisaCurves,k)
            SummarizedExperiment::colData(cells)[regionName] <- regions
        }
        
        if(is(cells, "SegmentedCells")){
        lisaCurves <- lisa(cells,
                           Rs = Rs,
                           BPPARAM = BPPARAM,
                           window = window,
                           window.length = window.length,
                           whichParallel = whichParallel,
                           sigma = sigma,
                           lisaFunc = lisaFunc,
                           minLambda = minLambda,
                           fast = fast)
        
        kM <- kmeans(lisaCurves,k)
        regions <- paste('region',kM$cluster,sep = '_')
        cellAnnotation(cells, regionName) <- regions
        }
        
     cells   
    
}
