################################################################################
#
# Generics for lisaClust
#
################################################################################


#' Region accessors for SegmentedCells
#'
#' Methods to access various components of the `SegmentedCells` object.
#'
#' @usage region(x, imageID = NULL, annot = FALSE)
#' @usage region(x, imageID = NULL) <- value
#'
#' @param x A `SegmentedCells` object.
#' @param imageID A vector of imageIDs to specifically extract.
#' @param annot Add cell annotation when selecting region information.
#' @param value The relevant information used to replace.
#'
#' @return DataFrame or a list of DataFrames
#' @name Region accessors
#'
#'
#' @examples
#' library(spicyR)
#' set.seed(51773)
#'
#' x <- round(c(runif(200),runif(200)+1,runif(200)+2,runif(200)+3,
#'              runif(200)+3,runif(200)+2,runif(200)+1,runif(200)),4)*100
#' y <- round(c(runif(200),runif(200)+1,runif(200)+2,runif(200)+3,
#'              runif(200),runif(200)+1,runif(200)+2,runif(200)+3),4)*100
#' cellType <- factor(paste('c',rep(rep(c(1:2),rep(200,2)),4),sep = ''))
#' imageID <- rep(c('s1', 's2'),c(800,800))
#' 
#' cells <- data.frame(x, y, cellType, imageID)
#'
#' cellExp <- SegmentedCells(cells, cellProfiler = TRUE)
#'
#'
#' # Generate LISA
#' lisaCurves <- lisa(cellExp)
#'
#' # Cluster the LISA curves
#' kM <- kmeans(lisaCurves,2)
#' region(cellExp) <- paste('region',kM$cluster,sep = '_')
#'
#' @aliases
#' region,SegmentedCells-method
#' region<-,SegmentedCells-method
#' region
#' region<-


### Get regions information

#' @export
#' @importFrom BiocGenerics do.call rbind
#' @importFrom spicyR cellSummary SegmentedCells
#' @importClassesFrom spicyR SegmentedCells
setGeneric("region", function(x, imageID = NULL, annot = FALSE)
    standardGeneric("region"))
setMethod("region", "SegmentedCells", function(x, imageID = NULL, annot = FALSE) {
    if (!is.null(imageID)) {
        x <- x[imageID,]
    }
    if (is.null(x$region))
        stop("There is no region information in your SegmentedCells yet")
    if (annot)
        return(data.frame(spicyR::cellSummary(x), region = BiocGenerics::do.call("rbind", x$region)))
    
    BiocGenerics::do.call("rbind", x$region)
})




#' @export
#' @importFrom S4Vectors DataFrame split
#' @importFrom methods new
#' @importFrom spicyR imageID SegmentedCells
#' @importClassesFrom spicyR SegmentedCells
setGeneric("region<-", function(x, imageID = NULL, value)
    standardGeneric("region<-"))
setReplaceMethod("region", "SegmentedCells", function(x, imageID = NULL, value) {
    if (is.null(imageID))
        imageID <- rownames(x)
    if (length(value) == length(imageID(x, imageID))) {
        if (is.null(x$region)) {
            x <- S4Vectors::DataFrame(x)
            by <- rep(rownames(x), unlist(lapply(x$cellSummary, nrow)))
            by <- factor(by, levels = unique(by))
            x$region <-
                S4Vectors::split(DataFrame(region = value), by)
            x <- new("SegmentedCells", x)
        }
        if (!is.null(x$region))
            by <- rep(rownames(x), unlist(lapply(x$cellSummary, nrow)))
        by <- factor(by, levels = unique(by))
        value <- S4Vectors::split(S4Vectors::DataFrame(region = value), by)
        x <- .putData(x, "region", value, imageID)
    }
    x
})

#' @importFrom methods slot slot<-
.putData <- function(object,variable, value, image = NULL){
    methods::slot(object, "listData")[[variable]] <- value
    object
}

