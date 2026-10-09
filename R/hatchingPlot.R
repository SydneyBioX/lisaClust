#' hatchingPlot
#'
#' The hatchingPlot() function is used to create hatching patterns for representating
#' spatial regions and cell-types.
#'
#' @param cells A data.frame or SingleCellExperiment.
#' @param useImages A vector of images to plot.
#' @param region The region column to plot.
#' @param imageID The imageIDs column if using data.frame or SingleCellExperiment.
#' @param cellType The cellType column if using data.frame or SingleCellExperiment.
#' @param spatialCoords The spatial coordinates columns if using data.frame or SingleCellExperiment.
#' @param window Should the window around the regions be 'square', 'convex' or 'concave'.
#' @param line.spacing A integer indicating the spacing between hatching lines.
#' @param hatching.colour Colour for the hatching.
#' @param nbp Not used: regions are outlined by the Voronoi tiles of their cells.
#' @param window.length A tuning parameter for controlling the level of concavity
#' when estimating concave windows.
#'
#' @details Each region is outlined by the union of the Voronoi tiles of its cells, clipped to the window, and
#' its hatching is clipped to that outline. Up to 12 regions have their own hatching: none, /, \\, -, |, x, +,
#' dots, circles, / with dots, \\ with dots and - with dots.
#'
#' @return A ggplot object
#'
#' @examples
#' ## Generate toy data
#' set.seed(51773)
#' x <- round(c(
#'     runif(200), runif(200) + 1, runif(200) + 2, runif(200) + 3,
#'     runif(200) + 3, runif(200) + 2, runif(200) + 1, runif(200)
#' ), 4) * 100
#' y <- round(c(
#'     runif(200), runif(200) + 1, runif(200) + 2, runif(200) + 3,
#'     runif(200), runif(200) + 1, runif(200) + 2, runif(200) + 3
#' ), 4) * 100
#' cellType <- factor(paste("c", rep(rep(c(1:2), rep(200, 2)), 4), sep = ""))
#' imageID <- rep(c("s1", "s2"), c(800, 800))
#' cells <- data.frame(x, y, cellType, imageID)
#' cells <- SingleCellExperiment::SingleCellExperiment(colData = cells)
#'
#' ## Generate regions
#' cells <- lisaClust(cells, k = 2)
#'
#' ## Plot regions
#' hatchingPlot(cells)
#'
#' @export
#' @rdname hatchingPlot
#' @importFrom ggplot2 ggplot aes geom_point theme_minimal facet_wrap labs
#' @importFrom dplyr .data
hatchingPlot <-
  function(cells,
           useImages = NULL,
           region = "region",
           imageID = "imageID",
           cellType = "cellType",
           spatialCoords = c("x", "y"),
           window = "concave",
           line.spacing = 21,
           hatching.colour = 1,
           nbp = NULL,
           window.length = NULL) {
    df <- .formatCells(cells, imageID, cellType, spatialCoords)
    if (!region %in% colnames(df)) stop(paste0("'", region, "' column not found in data"))
    df <- data.frame(df[, c("imageID", "cellID", "imageCellID", "x", "y", "cellType")], region = df[[region]])
    
    if (is.null(useImages)) useImages <- df$imageID[1]
    
    if (any(!useImages %in% df$imageID)) {
      stop("Some of the useImages are not in your data.")
    }
    
    df <- df[df$imageID %in% useImages, ]
    p <-
      ggplot(df, aes(
        x = .data$x,
        y = .data$y,
        colour = cellType
      )) +
      geom_point() +
      facet_wrap(~imageID) +
      geom_hatching(
        aes(region = region),
        show.legend = TRUE,
        window = window,
        line.spacing = line.spacing,
        hatching.colour = hatching.colour,
        window.length = window.length
      )
    q <- p + theme_minimal() + scale_region() + labs(x = "x", y = "y")
    q
  }


################################################################################
##
## geomHatching
##
################################################################################



#' Hatching geom
#'
#' The hatching geom is used to create hatching patterns for representation of spatial regions.
#'
#' @param mapping Set of aesthetic mappings created by aes() or aes_(). If specified
#' and inherit.aes = TRUE (the default), it is combined with the default mapping
#' at the top level of the plot. You must supply mapping if there is no plot mapping.
#' @param data The data to be displayed in this layer. There are three options:
#'
#' If NULL, the default, the data is inherited from the plot data as specified
#' in the call to ggplot(). A data.frame, or other object, will override the plot
#' data. All objects will be fortified to produce a data frame. See fortify() for
#'  which variables will be created. A function will be called with a single argument,
#'  the plot data. The return value must be a data.frame, and will be used as the
#'  layer data. A function can be created from a formula (e.g. ~ head(.x, 10)).
#' @param stat The statistical transformation to use on the data for this layer as a string.
#' @param position adjustment, either as a string, or the result of a call to a
#' position adjustment function.
#' @param show.legend logical. Should this layer be included in the legends? NA,
#' the default, includes if any aesthetics are mapped. FALSE never includes, and
#' TRUE always includes. It can also be a named logical vector to finely select
#' the aesthetics to display.
#' @param inherit.aes If FALSE, overrides the default aesthetics, rather than
#' combining with them. This is most useful for helper functions that define both
#' data and aesthetics and shouldn't inherit behaviour from the default plot
#' specification, e.g. borders().
#' @param na.rm If FALSE, the default, missing values are removed with a warning.
#' If TRUE, missing values are silently removed.
#' @param line.spacing A integer indicating the spacing between hatching lines.
#' @param hatching.colour A colour for the hatching.
#' @param window Should the window around the regions be 'square', 'convex' or 'concave'.
#' @param window.length A tuning parameter for controlling the level of concavity
#' when estimating concave windows.
#' @param nbp Not used: regions are outlined by the Voronoi tiles of their cells.
#' @param line.width A numeric controlling the width of the hatching lines
#' @param ... Other arguments passed on to layer(). These are often aesthetics,
#' used to set an aesthetic to a fixed value, like colour = "red" or size = 3.
#' They may also be parameters to the paired geom/stat.
#'
#'
#' @return A ggplot geom
#'
#' @examples
#' ## Generate toy data
#' set.seed(51773)
#' library(ggplot2)
#' x <- round(c(
#'     runif(200), runif(200) + 1, runif(200) + 2, runif(200) + 3,
#'     runif(200) + 3, runif(200) + 2, runif(200) + 1, runif(200)
#' ), 4) * 100
#' y <- round(c(
#'     runif(200), runif(200) + 1, runif(200) + 2, runif(200) + 3,
#'     runif(200), runif(200) + 1, runif(200) + 2, runif(200) + 3
#' ), 4) * 100
#' cellType <- factor(paste("c", rep(rep(c(1:2), rep(200, 2)), 4), sep = ""))
#' imageID <- rep(c("s1", "s2"), c(800, 800))
#' cells <- data.frame(x, y, cellType, imageID)
#' ## Generate regions
#' cells <- lisaClust(cells, k = 2)
#'
#' # Plot the regions with geom_hatching()
#' ggplot(
#'     cells, aes(x = x, y = y, colour = cellType, region = region)
#' ) +
#'     geom_point() +
#'     facet_wrap(~imageID) +
#'     geom_hatching()
#'
#' @export
#' @rdname hatchingPlot
#' @importFrom methods is
#' @importFrom BiocParallel bplapply
#' @importFrom ggplot2 layer
geom_hatching <-
  function(mapping = NULL,
           data = NULL,
           stat = "identity",
           position = "identity",
           na.rm = FALSE,
           show.legend = NA,
           inherit.aes = TRUE,
           line.spacing = 21,
           hatching.colour = 1,
           window = "concave",
           window.length = NULL,
           nbp = NULL,
           line.width = 1,
           ...) {
    ggplot2::layer(
      geom = GeomHatching,
      mapping = mapping,
      data = data,
      stat = stat,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(
        na.rm = na.rm,
        line.spacing = line.spacing,
        hatching.colour = hatching.colour,
        window = window,
        window.length = window.length,
        line.width = line.width,
        ...
      )
    )
  }

#' Scale constructor for regions
#'
#' Region scale constructor.
#'
#' @param aesthetics The names of the aesthetics that this scale works with
#' @param ... Arguments passed on to discrete_scale
#' @param guide A function used to create a guide or its name. See guides() for more info.
#'
#' @return a ggplot guide
#'
#' @examples
#'
#' ## Generate toy data
#' set.seed(51773)
#' x <- round(c(
#'     runif(200), runif(200) + 1, runif(200) + 2, runif(200) + 3,
#'     runif(200) + 3, runif(200) + 2, runif(200) + 1, runif(200)
#' ), 4) * 100
#' y <- round(c(
#'     runif(200), runif(200) + 1, runif(200) + 2, runif(200) + 3,
#'     runif(200), runif(200) + 1, runif(200) + 2, runif(200) + 3
#' ), 4) * 100
#' cellType <- factor(paste("c", rep(rep(c(1:2), rep(200, 2)), 4), sep = ""))
#' imageID <- rep(c("s1", "s2"), c(800, 800))
#' cells <- data.frame(x, y, cellType, imageID)
#' cells <- SingleCellExperiment::SingleCellExperiment(colData = cells)
#'
#' ## Generate regions
#' cells <- lisaClust(cells, k = 2)
#'
#' # Plot the regions with hatchingPlot()
#' hatchingPlot(cells) +
#'     scale_region_manual(
#'         values = c(1, 4), labels = c("Region A", "Region B"),
#'         name = "Regions"
#'     )
#'
#' @export
#' @rdname scale_region
#' @importFrom ggplot2 discrete_scale
scale_region <-
  function(aesthetics = "region",
           ...,
           guide = "legend") {
    discrete_scale(
      "region",
      "region_d",
      palette = function(n) 
        seq_len(n),
      ...
    )
  }

#' @export
#' @param values a set of aesthetic values to map data values to. If this is a
#' named vector, then the values will be matched based on the names. If unnamed,
#' values will be matched in order (usually alphabetical) with the limits of the scale.
#' Any data values that don't match will be given na.value.
#' @rdname scale_region
#' @importFrom ggplot2 discrete_scale
scale_region_manual <- function(..., values) {
  force(values)
  pal <- function(n) {
    if (n > length(values)) {
      stop(
        "Insufficient values in manual scale. ",
        n,
        " needed but only ",
        length(values),
        " provided.",
        call. = FALSE
      )
    }
    if (any(!values %in% seq_len(nHatchings))) {
      stop("values must be between 1 and ", nHatchings)
    }
    
    values
  }
  discrete_scale("region", "manual", pal, ...)
}


ggname <- getFromNamespace("ggname", "ggplot2")

#' @importFrom grid grob polylineGrob pointsGrob rectGrob gpar gTree gList unit
draw_key_region <- function(data, params, size) {
  # the region's hatching in a unit square, three repeats across, inside a frame
  hatching.colour <- params$hatching.colour
  if (is.null(hatching.colour)) hatching.colour <- 1
  type <- as.integer(data$region)
  square <- list(list(x = c(0, 1, 1, 0), y = c(0, 0, 1, 1)))
  kids <- hatchGrobs(type, square, 1 / 3, gpar(col = hatching.colour, lwd = 1, fill = NA), size = 0.8)
  kids[[length(kids) + 1]] <- rectGrob(gp = gpar(col = hatching.colour, lwd = 1, fill = NA))
  gTree(children = do.call(gList, kids), name = "region_key")
}



hatchingLevels <- function(data, hatching = NULL) {
  if (!is.factor(data$region)) {
    data$region <- factor(data$region)
  }
  regionLevels <- levels(data$region)
  if (!any(hatching %in% seq_len(nHatchings)) & !is.null(hatching)) {
    stop("hatching must equal the number of regions and be <= ", nHatchings, ".")
  }
  if (all(regionLevels %in% names(hatching))) {
    hatching <- hatching[regionLevels]
  }
  if (is.null(hatching)) {
    hatching <- seq_len(length(regionLevels))
    names(hatching) <- regionLevels
  }
  if (length(hatching) == length(regionLevels)) {
    names(hatching) <- regionLevels
  } else {
    stop("hatching must be the same length as the number of regions.")
  }
  return(hatching)
}




# The hatching grobs of the most recently drawn panels, by a hash of their inputs.
.hatchingCache <- local({
  store <- list()
  size <- 20
  list(
    get = function(key) store[[key]],
    set = function(key, value) {
      store[[key]] <<- value
      if (length(store) > size) store <<- store[seq(length(store) - size + 1, length(store))]
      invisible(value)
    },
    clear = function() store <<- list()
  )
})

GeomHatching <-
  ggplot2::ggproto(
    "GeomHatching",
    ggplot2::GeomPoint,
    extra_params = c(
      "na.rm",
      "line.spacing", "window", "window.length", "line.width", "hatching.colour"
    ),
    draw_panel = function(data,
                          panel_params,
                          coord,
                          na.rm = FALSE,
                          line.spacing = 21,
                          window = "convex",
                          window.length = NULL,
                          nbp = NULL,
                          line.width = 1,
                          hatching.colour = 1) {
      region <- data$region
      if (is.factor(region)) region <- as.numeric(region)
      if (is.character(region)) region <- as.numeric(as.factor(region))
      if (max(region) > nHatchings) {
        warning("Can not plot more than ", nHatchings, " regions. Adding regions greater than ", nHatchings,
                " to region 1.")
        region[region > nHatchings] <- 1
      }
      
      # The region outlines depend only on the cells, their regions and the window, so a plot that is printed
      # again (for example when its window is resized) reuses them.
      key <- rlang::hash(list(data$x, data$y, region, window, window.length))
      polys <- .hatchingCache$get(key)
      if (is.null(polys)) {
        polys <- regionPolygons(data$x, data$y, region, makeWindow(data, window, window.length))
        .hatchingCache$set(key, polys)
      }
      
      # outlines in the panel's coordinates
      polys <- lapply(polys, function(rings) lapply(rings, function(r) {
        xy <- coord$transform(data.frame(x = r$x, y = r$y), panel_params)
        list(x = xy$x, y = xy$y)
      }))
      grid::gTree(polys = polys, spacing = 1 / line.spacing, col = hatching.colour, lwd = line.width,
                  name = "geom_hatching", cl = "lisaHatching")
    },
    draw_key = draw_key_region,
    required_aes = c("x", "y", "region"),
    non_missing_aes = c(
      "x",
      "y", "region"
    ),
    default_aes = ggplot2::aes(
      region = 0,
      size = 0.05,
      angle = 0,
      alpha = 1
    )
  )





################################################################################
##
## Region outlines and their hatching
##
################################################################################


# The outline of each region in data coordinates: the union of the Voronoi tiles of its cells, clipped to the
# window. A list with one element per region code (1 to 7), each a list of rings list(x, y).
#' @importFrom deldir deldir
#' @importFrom polyclip polyclip
regionPolygons <- function(x, y, region, window) {
  win <- if (window$type == "rectangle") {
    list(list(x = window$xrange[c(1, 2, 2, 1)], y = window$yrange[c(1, 1, 2, 2)]))
  } else {
    lapply(window$bdry, function(b) list(x = b$x, y = b$y))
  }
  keep <- !duplicated(cbind(x, y))
  x <- x[keep]; y <- y[keep]; region <- region[keep]
  out <- vector("list", nHatchings)
  if (length(unique(region)) == 1 || length(x) < 3) {
    out[[region[1]]] <- win
    return(out)
  }
  bb <- c(range(x), range(y))
  pad <- 0.05 * max(diff(bb[1:2]), diff(bb[3:4]))
  rw <- c(bb[1] - pad, bb[2] + pad, bb[3] - pad, bb[4] + pad)
  dd <- suppressWarnings(deldir::deldir(x, y, rw = rw, round = FALSE))
  e <- dd$dirsgs
  # A Voronoi tile is convex, so its corners in order of angle around its cell trace it: the ends of the edges
  # it shares with its neighbours, and the corners of the frame closest to its cell.
  gen <- c(e$ind1, e$ind1, e$ind2, e$ind2)
  vx <- c(e$x1, e$x2, e$x1, e$x2)
  vy <- c(e$y1, e$y2, e$y1, e$y2)
  cx <- rw[c(1, 2, 2, 1)]
  cy <- rw[c(3, 3, 4, 4)]
  gen <- c(gen, .nearestLabels(x, y, seq_along(x), cx, cy))
  vx <- c(vx, cx)
  vy <- c(vy, cy)
  d <- data.frame(gen = gen, vx = vx, vy = vy)
  d <- d[!duplicated(data.frame(d$gen, round(d$vx, 9), round(d$vy, 9))), ]
  d <- d[order(d$gen, atan2(d$vy - y[d$gen], d$vx - x[d$gen])), ]
  tiles <- lapply(split(d, d$gen), function(t) list(x = t$vx, y = t$vy))
  tileRegion <- region[as.integer(names(tiles))]
  for (r in unique(tileRegion)) {
    tl <- tiles[tileRegion == r]
    u <- if (length(tl) == 1) tl else polyclip::polyclip(tl[1], tl[-1], op = "union", fillA = "nonzero", fillB = "nonzero")
    out[[r]] <- polyclip::polyclip(u, win, op = "intersection")
  }
  out
}


# The number of hatchings: none, /, \, -, |, x, +, dots, circles, / with dots, \ with dots, - with dots.
nHatchings <- 12

# Line segments for hatching type `type`, as list(x0, y0, x1, y1) in units of the spacing, within one tile;
# repeated every tile, they join into continuous lines. NULL for a hatching without lines.
hatchSegments <- function(type) {
  type <- c(`10` = 2, `11` = 3, `12` = 4)[as.character(type)] %|NA|% type
  switch(as.character(type),
    `2` = list(x0 = c(-1, 0, 1), y0 = c(0, 0, 0), x1 = c(0, 1, 2), y1 = c(1, 1, 1)),
    `3` = list(x0 = c(-1, 0, 1), y0 = c(1, 1, 1), x1 = c(0, 1, 2), y1 = c(0, 0, 0)),
    `4` = list(x0 = 0, y0 = 0.5, x1 = 1, y1 = 0.5),
    `5` = list(x0 = 0.5, y0 = 0, x1 = 0.5, y1 = 1),
    `6` = list(x0 = c(-1, 0, 1, -1, 0, 1), y0 = c(0, 0, 0, 1, 1, 1), x1 = c(0, 1, 2, 0, 1, 2), y1 = c(1, 1, 1, 0, 0, 0)),
    `7` = list(x0 = c(0, 0.5), y0 = c(0.5, 0), x1 = c(1, 0.5), y1 = c(0.5, 1)),
    NULL
  )
}

`%|NA|%` <- function(a, b) if (is.na(a)) b else unname(a)

# The marks of hatching type `type`: "dot", "circle" or NULL. Dots or circles alone sit on a staggered grid, two
# per tile at (0.5, 0) and (0, 0.5); with lines, one per tile at (0.5, 0) in units of the spacing, which lies
# between the lines of /, \ and -.
hatchMarks <- function(type) {
  if (type %in% c(8, 10, 11, 12)) "dot" else if (type == 9) "circle" else NULL
}

# The hatching segments of type `type` over the unit square: the tile's segments repeated every w.
hatchLines <- function(type, w) {
  s <- hatchSegments(type)
  if (is.null(s)) return(list())
  k <- seq(-1, ceiling(1 / w) + 1)
  off <- expand.grid(i = k, j = k)
  lapply(seq_len(nrow(off) * length(s$x0)), function(m) {
    o <- off[(m - 1) %/% length(s$x0) + 1, ]
    q <- (m - 1) %% length(s$x0) + 1
    list(x = (c(s$x0[q], s$x1[q]) + o$i) * w, y = (c(s$y0[q], s$y1[q]) + o$j) * w)
  })
}

# Whether points are strictly inside the rings, by the even-odd rule (holes are rings inside rings).
#' @importFrom polyclip pointinpolygon
insideRings <- function(x, y, rings) {
  odd <- logical(length(x))
  for (r in rings) odd <- xor(odd, polyclip::pointinpolygon(list(x = x, y = y), r) == 1)
  odd
}

# The grobs of hatching type `type` clipped to `rings`: its lines, clipped by polyclip, and its marks inside.
hatchGrobs <- function(type, rings, w, gp, size = 1) {
  kids <- list()
  if (type <= 1) return(kids)
  lines <- hatchLines(type, w)
  if (length(lines)) {
    pieces <- polyclip::polyclip(lines, rings, op = "intersection", closed = FALSE)
    if (length(pieces)) {
      kids[[length(kids) + 1]] <- polylineGrob(
        unlist(lapply(pieces, `[[`, "x")), unlist(lapply(pieces, `[[`, "y")),
        id = rep(seq_along(pieces), lengths(lapply(pieces, `[[`, "x"))), gp = gp
      )
    }
  }
  marks <- hatchMarks(type)
  if (!is.null(marks)) {
    k <- seq(-1, ceiling(1 / w) + 1)
    g <- expand.grid(i = k, j = k)
    px <- (g$i + 0.5) * w
    py <- g$j * w
    if (type %in% c(8, 9)) {
      px <- c(px, g$i * w)
      py <- c(py, (g$j + 0.5) * w)
    }
    keep <- insideRings(px, py, rings)
    if (any(keep)) {
      lwd <- if (is.null(gp$lwd)) 1 else gp$lwd
      kids[[length(kids) + 1]] <- pointsGrob(
        px[keep], py[keep], pch = if (marks == "dot") 16 else 1,
        size = unit(size * (if (marks == "dot") 3 else 4.5) * lwd, "pt"), gp = gp
      )
    }
  }
  kids
}

#' @importFrom grid makeContent gList pathGrob polylineGrob gpar
#' @exportS3Method grid::makeContent
makeContent.lisaHatching <- function(x) {
  # Each region's hatching, clipped to its outline, and the outline itself.
  gp <- gpar(col = x$col, lwd = x$lwd, fill = NA, lineend = "butt")
  kids <- list()
  for (type in which(lengths(x$polys) > 0)) {
    rings <- x$polys[[type]]
    kids <- c(kids, hatchGrobs(type, rings, x$spacing, gp))
    kids[[length(kids) + 1]] <- pathGrob(
      unlist(lapply(rings, `[[`, "x")), unlist(lapply(rings, `[[`, "y")),
      id = rep(seq_along(rings), lengths(lapply(rings, `[[`, "x"))), rule = "evenodd", gp = gp
    )
  }
  grid::setChildren(x, do.call(gList, kids))
}
