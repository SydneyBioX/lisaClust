#' Generate local indicators of spatial association
#'
#' @param cells A SingleCellExperiment, SpatialExperiment or data frame that contains at least the
#' variables x and y, giving the  coordinates of each cell, imageID and cellType.
#' @param r A vector of the radii that the measures of association should be calculated.
#' @param imageID The column which contains image identifiers.
#' @param cellType The column which contains the cell types.
#' @param spatialCoords The columns which contain the x and y spatial coordinates.
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
#' @param BPPARAM \{DEPRECATED\} A BiocParallel MulticoreParam or SerialParam object. 
#' @param Rs \{DEPRECATED\} A vector of the radii that the measures of association should be calculated.
#' 
#' @return A matrix of LISA curves
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
#' # Generate LISA curves
#' lisaCurves <- lisa(cells, r = c(10, 20, 50))
#'
#' # Cluster the LISA curves
#' kM <- kmeans(lisaCurves, 2)
#'
#' @export
#' @rdname lisa
#' @importFrom methods is
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom S4Vectors DataFrame
#' @importFrom BiocGenerics do.call rbind
#' @importFrom dplyr bind_rows
lisa <- function(cells,
                 r = NULL,
                 imageID = "imageID",
                 cellType = "cellType",
                 spatialCoords = c("x", "y"),
                 cores = 1,
                 window = "convex",
                 window.length = NULL,
                 whichParallel = "imageID",
                 sigma = NULL,
                 lisaFunc = "K",
                 minLambda = 0.05,
                 BPPARAM = NULL,
                 Rs = r) {
  
  user_args = as.list(match.call())[-1]
  
  tryCatch({
    user_vals = lapply(names(user_args), function(arg) {
      if (arg %in% c("BPPARAM", "cores", "Rs")) {
        eval(user_args[[arg]])
      } 
    })
    
    names(user_vals) = names(user_args)
    
    if (length(user_vals) > 0 && any(!sapply(user_vals, is.null))) {
      argumentChecks("lisaClust", user_vals)
    }
  }, error = function(e) {
    if (grepl("object 'cd' not found", e$message)) {
      message("Skipping argument checks as `lisa()` is being called within `lisaClust()`")
    } else {
      stop(e)
    }
  })
  
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
  
  cells <- .formatCells(cells, imageID, cellType, spatialCoords)
  cellSummary <- .cellsByImage(cells)
  
  if (is.null(Rs)) {
    Rs <- c(20, 50, 100)
  }
  
  BPPARAM <- .bpparam(cores)
  
  message("Generating local indicators of spatial association.")
  
  curveList <-
    BiocParallel::bplapply(
      cellSummary,
      inhomLocalK,
      Rs = Rs,
      sigma = sigma,
      window = window,
      window.length = window.length,
      minLambda = minLambda,
      lisaFunc = lisaFunc,
      BPPARAM = BPPARAM
    )
  
  curvelist <- lapply(curveList, as.data.frame)
  curves <- as.matrix(dplyr::bind_rows(curvelist))
  rownames(curves) <- as.character(unlist(lapply(cellSummary, function(x) x$cellID)))
  
  curves[is.na(curves)] <- 0
  return(curves)
}




#' @importFrom spatstat.geom ppp
pppGenerate <- function(cells, window, window.length) {
  ow <- makeWindow(cells, window, window.length)
  pppCell <- spatstat.geom::ppp(
    cells$x,
    cells$y,
    window = ow,
    marks = cells$cellType
  )
  
  pppCell
}

#' @importFrom spatstat.geom owin convexhull ppp
#' @importFrom concaveman concaveman
makeWindow <-
  function(data,
           window = "square",
           window.length = NULL) {
    data <- data.frame(data)
    ow <-
      spatstat.geom::owin(xrange = range(data$x), yrange = range(data$y))
    
    if (window == "convex") {
      p <- spatstat.geom::ppp(data$x, data$y, ow)
      ow <- spatstat.geom::convexhull(p)
    }
    if (window == "concave") {
      message("Concave windows are temperamental. Try choosing values of window.length > and < 1 if you have problems.")
      if (is.null(window.length)) {
        window.length <- (max(data$x) - min(data$x)) / 20
      } else {
        window.length <- (max(data$x) - min(data$x)) / 20 * window.length
      }
      dist <- (max(data$x) - min(data$x)) / (length(data$x))
      # each cell and its 8 neighbours at distance dist, cell by cell
      bigDat <- cbind(
        rep(data$x, each = 9) + rep(c(0, 1, 0, -1, -1, 0, 1, -1, 1) * dist, times = nrow(data)),
        rep(data$y, each = 9) + rep(c(0, 1, 1, 1, -1, -1, -1, 0, 0) * dist, times = nrow(data))
      )
      ch <-
        concaveman::concaveman(bigDat,
                               length_threshold = window.length,
                               concavity = 1
        )
      poly <- as.data.frame(ch[nrow(ch):1, ])
      colnames(poly) <- c("x", "y")
      ow <-
        spatstat.geom::owin(
          xrange = range(poly$x),
          yrange = range(poly$y),
          poly = poly
        )
    }
    ow
  }


#' @importFrom spatstat.geom union.owin border inside.owin
#' @useDynLib lisaClust, .registration = TRUE
#' @importFrom Rcpp sourceCpp
borderEdge <- function(X, maxD) {
  W <- X$window
  bW <- spatstat.geom::union.owin(
    spatstat.geom::border(W, maxD, outside = FALSE),
    spatstat.geom::border(W, 2, outside = TRUE)
  )
  inB <- spatstat.geom::inside.owin(X$x, X$y, bW)
  e <- rep(1, X$n)
  if (any(inB)) {
    # area(intersect.owin(discs(X[inB], maxD), W)) / (pi maxD^2): each disc is the 128-gon
    # spatstat.geom::disc() builds, clipped to the window in C++.
    rings <- if (W$type == "rectangle") {
      list(list(x = W$xrange[c(1, 2, 2, 1)], y = W$yrange[c(1, 1, 2, 2)]))
    } else {
      W$bdry
    }
    e[inB] <- .discWindowArea(X$x[inB], X$y[inB], maxD, 128L, rings) / (pi * maxD^2)
  }
  e
}




#' Calculate the inhomogenous local K function.
#'
#' @param data The data.
#' @param Rs
#'   A vector of the radii that the measures of association should be
#'   calculated.
#' @param sigma
#'   A numeric variable used for scaling when filting inhomogeneous L-curves.
#' @param window
#'   Should the window around the regions be 'square', 'convex' or 'concave'.
#' @param window.length
#'   A tuning parameter for controlling the level of concavity.
#' @param minLambda
#'   Minimum value for density for scaling when fitting inhomogeneous L-curves.
#' @param lisaFunc Either "K" or "L" curve.
#'
#' @return A matrix of LISA curves
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
#' cells$cellID <- paste0("cell_", seq_len(nrow(cells)))
#'
#' # The curves of the first image
#' inhom <- inhomLocalK(cells[cells$imageID == "s1", ], Rs = c(10, 20, 50))
#'
#' @export
#' @rdname inhomLocalK
#' @importFrom spatstat.geom ppp area nearest.valid.pixel
#' @importFrom spatstat.explore density.ppp
#' @importFrom spatstat.random expand.owin
inhomLocalK <-
  function(data,
           Rs = c(20, 50, 100, 200),
           sigma = 10000,
           window = "convex",
           window.length = NULL,
           minLambda = 0.05,
           lisaFunc = "K") {
    ow <- makeWindow(data, window, window.length)
    ow <- spatstat.random::expand.owin(ow, distance = 0.01)
    X <-
      spatstat.geom::ppp(
        x = data$x,
        y = data$y,
        window = ow,
        marks = data$cellType
      )
    
    if (is.null(Rs)) {
      Rs <- c(20, 50, 100, 200)
    }
    if (is.null(sigma)) {
      sigma <- 100000
    }
    
    maxR <- min(ow$xrange[2] - ow$xrange[1], ow$yrange[2] - ow$yrange[1]) / 2.01
    Rs <- unique(pmin(c(0, sort(Rs)), maxR))
    
    den <- spatstat.explore::density.ppp(X, sigma = sigma)
    den <- den / mean(den)
    den$v <- pmax(den$v, minLambda)
    
    # inverse-density weight of each cell as a neighbour
    np <- spatstat.geom::nearest.valid.pixel(X$x, X$y, den)
    w <- den$v[cbind(np$row, np$col)]
    wt <- 1 / w * mean(w)
    rm(np)
    
    cellType <- data$cellType
    if (!is.factor(cellType)) cellType <- factor(cellType)
    lam <- as.numeric(table(cellType)) / spatstat.geom::area(X)
    edge <- vapply(Rs[-1], function(x) borderEdge(X, x), numeric(X$n))
    edge <- matrix(edge, nrow = X$n)
    labels <- as.character(Rs[-1])
    
    res <- .localCurves(
      as.numeric(data$x), as.numeric(data$y), as.integer(cellType), nlevels(cellType), as.numeric(Rs),
      as.numeric(labels), as.numeric(wt), lam, edge, lisaFunc == "L", includeSelf = TRUE
    )
    
    # one column per radius and neighbouring type, as radius_type, for the radii and types that occur
    keep <- which(res$type)
    curves <- do.call("cbind", c(list(matrix(numeric(0), nrow = X$n, ncol = 0)), lapply(which(res$bin), function(k) {
      m <- matrix(res$value[, keep, k], nrow = X$n)
      colnames(m) <- paste(labels[k], levels(cellType)[keep], sep = "_")
      m
    })))
    curves[!res$cell, ] <- NA
    rownames(curves) <- data$cellID
    curves
  }


#' Plot heatmap of cell type enrichment for lisaClust regions
#'
#' @param cells SingleCellExperiment, SpatialExperiment or data.frame
#' @param type Make a "bubble" or "heatmap" plot, or return the relative frequencies as a "table".
#' @param region The column storing the regions
#' @param cellType The column storing the cell types
#' @param limit limits to the lower and upper relative frequencies
#' @param ... Any arguments to be passed to the pheatmap package
#'
#' @return A bubble plot or heatmap, or with \code{type = "table"} a matrix of the relative frequencies with
#' one row per cell type and one column per region: how much more often the cell type is found in the region
#' than if cell types were spread evenly over the regions.
#'
#'
#' @examples
#'
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
#'
#' cells <- data.frame(x, y, cellType, imageID)
#'
#' cells <- lisaClust(cells, k = 2)
#'
#' regionMap(cells)
#' regionMap(cells, type = "table")
#'
#' @export
#' @importFrom SummarizedExperiment colData
#' @importFrom pheatmap pheatmap
#' @importFrom ggplot2 ggplot aes geom_point scale_colour_gradient2 theme_minimal labs
#' @importFrom dplyr mutate .data
#' @import SpatialExperiment SingleCellExperiment
regionMap <- function(cells, type = "bubble", cellType = "cellType", region = "region", limit = c(0.33, 3), ...) {
  tab <- .regionEnrichment(cells, cellType, region)
  if (type == "table") return(unclass(tab))
  
  ph <- pheatmap::pheatmap(pmax(pmin(tab, limit[2]), limit[1]), cluster_cols = FALSE, silent = TRUE, ...)
  
  if (type == "bubble") {
    p1 <- tab |>
      as.data.frame() |>
      dplyr::mutate(cellType = factor(.data$Var1, levels = levels(.data$Var1)[ph$tree_row$order]), region = .data$Var2,
                    Freq2 = pmax(pmin(.data$Freq, limit[2]), limit[1])) |>
      ggplot2::ggplot(ggplot2::aes(x = .data$region, y = .data$cellType, colour = .data$Freq2, size = .data$Freq2)) +
      ggplot2::geom_point() +
      ggplot2::scale_colour_gradient2(low = "#4575B4", mid = "grey90", high = "#D73027", midpoint = 1, guide = "legend") +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Region", y = "Cell-type", colour = "Relative\nFrequency", size = "Relative\nFrequency")
    
    return(p1)
  }
  
  pheatmap::pheatmap(pmax(pmin(tab, limit[2]), limit[1]), cluster_cols = FALSE, ...)
}


# The relative frequency of each cell type (rows) in each region (columns): observed count over the count
# expected if cell types were spread evenly over the regions.
.regionEnrichment <- function(cells, cellType, region) {
  df <- .colDataFrame(cells, c(cellType, region))
  tab <- table(df[, cellType], df[, region])
  # cell types or regions without cells (unused factor levels) have no enrichment
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  tab / rowSums(tab) %*% t(colSums(tab)) * sum(tab)
}
