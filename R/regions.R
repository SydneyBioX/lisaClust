#' Name regions by the cell types they are enriched for
#'
#' The regions found by \code{lisaClust()} are numbered arbitrarily. \code{nameRegions()} gives each region the
#' name of the marker cell type it is most enriched for, as measured by \code{regionMap(type = "table")}.
#' Regions most enriched for the same marker get the same name, so several regions can be merged into one
#' domain.
#'
#' @param cells A SingleCellExperiment, SpatialExperiment or data frame with a region column.
#' @param markers A character vector of cell types. Its names, when given, are used as the names of the
#' regions; otherwise the cell types are.
#' @param region The column storing the regions.
#' @param cellType The column storing the cell types.
#' @param regionName The output column for the named regions; by default the regions are replaced.
#'
#' @return \code{cells} with the named regions in the column \code{regionName}.
#'
#' @examples
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
#' cells <- lisaClust(cells, k = 2)
#'
#' cells <- nameRegions(cells, c(first = "c1", second = "c2"), regionName = "domain")
#' table(cells$domain, cells$cellType)
#'
#' @export
#' @importFrom methods is
nameRegions <- function(cells, markers, region = "region", cellType = "cellType", regionName = region) {
  if (is.null(names(markers))) names(markers) <- markers
  enrichment <- .regionEnrichment(cells, cellType, region)
  absent <- setdiff(markers, rownames(enrichment))
  if (length(absent)) stop("cell type(s) not found in '", cellType, "': ", paste(absent, collapse = ", "))
  nameOf <- names(markers)[apply(enrichment[markers, , drop = FALSE], 2, which.max)]
  names(nameOf) <- colnames(enrichment)
  named <- unname(nameOf[as.character(.colDataFrame(cells, region)[[1]])])
  if (is(cells, "SummarizedExperiment")) {
    SummarizedExperiment::colData(cells)[regionName] <- named
  } else {
    cells[regionName] <- named
  }
  cells
}


#' Compare the share of each region between groups
#'
#' For each image, or each patient when \code{imageID} names a patient column, the share of its cells in each
#' region, drawn as box plots by \code{condition} with one panel per region. An image or patient with no cells
#' in a region has a share of zero.
#'
#' @param cells A SingleCellExperiment, SpatialExperiment or data frame with a region column.
#' @param condition The column storing the group of each image (or patient).
#' @param region The column storing the regions.
#' @param imageID The column identifying the units whose shares are compared: images, or patients to pool the
#' images of each patient.
#' @param regions The regions to show; all by default.
#'
#' @return A ggplot object.
#'
#' @examples
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
#' cells <- data.frame(x, y, cellType, imageID, group = rep(c("A", "B"), c(800, 800)))
#' cells <- lisaClust(cells, k = 2)
#'
#' regionBoxPlot(cells, condition = "group")
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_point position_jitter facet_wrap scale_y_continuous labs
#'   theme element_text
regionBoxPlot <- function(cells, condition, region = "region", imageID = "imageID", regions = NULL) {
  df <- .colDataFrame(cells, c(imageID, condition, region))
  units <- unique(df[, c(imageID, condition)])
  if (anyDuplicated(units[[imageID]])) stop("each '", imageID, "' must have a single '", condition, "'")
  shares <- prop.table(table(df[[imageID]], df[[region]]), 1)
  if (!is.null(regions)) shares <- shares[, regions, drop = FALSE]
  shares <- as.data.frame(shares, responseName = "share", stringsAsFactors = FALSE)
  colnames(shares)[1:2] <- c("unit", "region")
  shares$condition <- units[[condition]][match(shares$unit, units[[imageID]])]
  shares <- shares[!is.na(shares$condition), , drop = FALSE]

  ggplot2::ggplot(shares, ggplot2::aes(x = .data$condition, y = .data$share, colour = .data$condition)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.15, height = 0, seed = 51773), size = 1) +
    ggplot2::facet_wrap(~ region, nrow = 1) +
    ggplot2::scale_y_continuous(labels = function(v) paste0(100 * v, "%")) +
    ggplot2::labs(x = NULL, y = "share of cells") +
    ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
