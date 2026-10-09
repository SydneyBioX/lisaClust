## The C++ core against plain-R references written from the definitions.

sim <- function(n = 300, seed = 1) {
  set.seed(seed)
  data.frame(x = round(runif(n, 0, 200), 1), y = round(runif(n, 0, 150), 1),
             cellType = factor(sample(c("a", "b", "c"), n, TRUE, prob = c(0.5, 0.3, 0.2))))
}

# The local indicators from a distance matrix: weighted counts per distance band, accumulated over the bands
# with any pair, and the expectation with the cell's edge share where it has a neighbour of that type in the band.
refCurves <- function(d, Rs, wt, lam, edge, L = FALSE, self = TRUE) {
  n <- nrow(d); K <- nlevels(d$cellType); nb <- length(Rs) - 1
  # with self, each cell is its own neighbour at distance 0
  D <- as.matrix(dist(d[, c("x", "y")])); diag(D) <- if (self) 0 else Inf
  band <- matrix(findInterval(D, Rs, left.open = TRUE, rightmost.closed = FALSE), n)
  band[D == 0] <- 1; band[D > max(Rs)] <- NA
  type <- as.integer(d$cellType)
  S <- array(0, c(n, K, nb)); has <- array(FALSE, c(n, K, nb))
  for (i in seq_len(n)) for (j in which(!is.na(band[i, ]))) {
    S[i, type[j], band[i, j]] <- S[i, type[j], band[i, j]] + wt[j]; has[i, type[j], band[i, j]] <- TRUE
  }
  cellPresent <- apply(has, 1, any); typePresent <- apply(has, 2, any); binPresent <- apply(has, 3, any)
  out <- array(NA_real_, c(n, K, nb))
  for (i in which(cellPresent)) for (J in which(typePresent)) {
    cum <- 0
    for (k in which(binPresent)) {
      cum <- cum + S[i, J, k]
      e <- if (has[i, J, k]) edge[i, k] else 1
      E <- Rs[k + 1]^2 * pi * e * lam[J]
      out[i, J, k] <- if (L) sqrt(cum) - sqrt(E) else (cum - E) / sqrt(E)
    }
  }
  list(value = out, cell = cellPresent, type = typePresent, bin = binPresent)
}

test_that("the local curves match the plain-R reference", {
  for (seed in 1:3) for (L in c(FALSE, TRUE)) for (self in c(TRUE, FALSE)) {
    d <- sim(seed = seed)
    Rs <- c(0, 10, 25, 40)
    wt <- runif(nrow(d), 0.5, 1.5); lam <- as.numeric(table(d$cellType)) / (200 * 150)
    edge <- matrix(runif(nrow(d) * 3, 0.5, 1), nrow(d))
    got <- lisaClust:::.localCurves(d$x, d$y, as.integer(d$cellType), 3L, Rs, Rs[-1], wt, lam, edge, L, self)
    ref <- refCurves(d, Rs, wt, lam, edge, L, self)
    expect_equal(got$cell, ref$cell)
    expect_equal(got$type, ref$type)
    expect_equal(got$bin, ref$bin)
    expect_equal(got$value, ref$value, tolerance = 1e-12)
  }
})

test_that("pairs at exactly a radius fall in the band below it, and coincident cells are neighbours", {
  d <- data.frame(x = c(0, 3, 0, 0), y = c(0, 4, 0, 10), cellType = factor(c("a", "b", "b", "a")))
  got <- lisaClust:::.localCurves(d$x, d$y, as.integer(d$cellType), 2L, c(0, 5, 10), c(5, 10), rep(1, 4),
                                  c(1, 1), matrix(1, 4, 2), FALSE)
  ref <- refCurves(d, c(0, 5, 10), rep(1, 4), c(1, 1), matrix(1, 4, 2))
  expect_equal(got$value, ref$value)
})

test_that("disc areas match spatstat's polygon intersection", {
  d <- sim(200)
  W <- spatstat.geom::convexhull(spatstat.geom::ppp(d$x, d$y, spatstat.geom::owin(c(0, 200), c(0, 150))))
  X <- spatstat.geom::ppp(d$x, d$y, window = W)
  for (r in c(5, 30)) {
    ref <- vapply(seq_len(X$n), function(i) {
      spatstat.geom::area(spatstat.geom::intersect.owin(spatstat.geom::disc(r, c(X$x[i], X$y[i]), npoly = 128), W))
    }, numeric(1))
    got <- lisaClust:::.discWindowArea(X$x, X$y, r, 128L, W$bdry)
    expect_equal(got, ref, tolerance = 1e-6)
  }
})

test_that("nearest labels match a brute-force search", {
  set.seed(3)
  tr <- data.frame(x = runif(500), y = runif(500), l = sample(1:4, 500, TRUE))
  q <- data.frame(x = runif(2000, -0.1, 1.1), y = runif(2000, -0.1, 1.1))
  got <- lisaClust:::.nearestLabels(tr$x, tr$y, tr$l, q$x, q$y)
  ref <- apply(q, 1, function(p) tr$l[which.min((tr$x - p[1])^2 + (tr$y - p[2])^2)])
  expect_equal(got, ref)
  # equidistant cells vote, and a tied vote goes to the first label
  expect_equal(lisaClust:::.nearestLabels(c(0, 2), c(0, 0), c(2L, 1L), 1, 0), 1L)
  expect_equal(lisaClust:::.nearestLabels(c(0, 2, 1), c(0, 0, 1), c(2L, 1L, 2L), 1, 0), 2L)
})

ringsArea <- function(rings) sum(vapply(rings, function(r) {
  n <- length(r$x); 0.5 * sum(r$x * r$y[c(2:n, 1)] - r$x[c(2:n, 1)] * r$y)
}, numeric(1)))

test_that("region outlines partition the window, and each cell is inside its own region", {
  set.seed(4)
  n <- 600
  x <- runif(n, 0, 100); y <- runif(n, 0, 80)
  region <- ifelse(x + 20 * sin(y / 8) < 50, 1, ifelse(y < 40, 2, 3))
  for (window in c("square", "convex", "concave")) {
    W <- suppressMessages(lisaClust:::makeWindow(data.frame(x = x, y = y), window, NULL))
    polys <- lisaClust:::regionPolygons(x, y, region, W)
    expect_equal(sum(vapply(polys[lengths(polys) > 0], ringsArea, numeric(1))), spatstat.geom::area(W), tolerance = 1e-6)
    # cells on the window's edge (the corners of a convex hull or bounding box) lie on the outline itself
    interior <- spatstat.geom::bdist.points(spatstat.geom::ppp(x, y, window = W, check = FALSE)) > 1e-9
    for (r in 1:3) {
      Wr <- spatstat.geom::owin(poly = polys[[r]], check = FALSE)
      k <- region == r & interior
      expect_true(all(spatstat.geom::inside.owin(x[k], y[k], Wr)))
    }
  }
})

test_that("hatching lines clipped to a region stay inside it", {
  ring <- list(list(x = c(0.1, 0.9, 0.9, 0.5, 0.1), y = c(0.1, 0.1, 0.9, 0.5, 0.9)))
  for (type in 2:7) {
    pieces <- polyclip::polyclip(lisaClust:::hatchLines(type, 1 / 20), ring, op = "intersection", closed = FALSE)
    expect_true(length(pieces) > 0)
    W <- spatstat.geom::owin(poly = ring, check = FALSE)
    mids <- t(vapply(pieces, function(p) c(mean(range(p$x)), mean(range(p$y))), numeric(2)))
    expect_true(all(spatstat.geom::inside.owin(mids[, 1], mids[, 2], spatstat.geom::grow.rectangle(spatstat.geom::Frame(W), 1e-9)) ))
    expect_true(all(unlist(lapply(pieces, `[[`, "x")) >= 0.1 - 1e-9 & unlist(lapply(pieces, `[[`, "x")) <= 0.9 + 1e-9))
  }
})

test_that("all 12 hatchings draw inside their region, with legend keys", {
  ring <- list(list(x = c(0.1, 0.9, 0.9, 0.5, 0.1), y = c(0.1, 0.1, 0.9, 0.5, 0.9)))
  W <- spatstat.geom::owin(poly = ring, check = FALSE)
  gp <- grid::gpar(col = "black", lwd = 1)
  expect_length(lisaClust:::hatchGrobs(1, ring, 1 / 20, gp), 0)
  for (type in 2:12) {
    kids <- lisaClust:::hatchGrobs(type, ring, 1 / 20, gp)
    expect_true(length(kids) > 0)
    for (k in kids) {
      if (inherits(k, "points")) {
        expect_true(all(spatstat.geom::inside.owin(as.numeric(k$x), as.numeric(k$y), W)))
      }
    }
    expect_equal(any(vapply(kids, inherits, TRUE, "points")), type >= 8)
    expect_equal(any(vapply(kids, inherits, TRUE, "polyline")), !type %in% c(8, 9))
    key <- lisaClust:::draw_key_region(data.frame(region = type), list(hatching.colour = 1), NULL)
    expect_s3_class(key, "gTree")
  }
  # a hole is left empty
  holed <- list(list(x = c(0, 1, 1, 0), y = c(0, 0, 1, 1)), list(x = c(0.3, 0.3, 0.7, 0.7), y = c(0.3, 0.7, 0.7, 0.3)))
  pts <- Filter(function(k) inherits(k, "points"), lisaClust:::hatchGrobs(8, holed, 1 / 20, gp))[[1]]
  px <- as.numeric(pts$x); py <- as.numeric(pts$y)
  expect_false(any(px > 0.3 & px < 0.7 & py > 0.3 & py < 0.7))
})

test_that("12 regions each get a hatching; a 13th falls back to the first with a warning", {
  set.seed(7)
  d <- data.frame(x = runif(1300, 0, 130), y = runif(1300, 0, 50))
  d$region <- sprintf("region_%02d", pmin(floor(d$x / 10) + 1, 13))
  p12 <- ggplot2::ggplot(d[d$x < 120, ], ggplot2::aes(x, y, region = region)) + geom_hatching(window = "square") +
    scale_region()
  pdf(NULL); on.exit(dev.off())
  expect_silent(print(p12))
  p13 <- ggplot2::ggplot(d, ggplot2::aes(x, y, region = region)) + geom_hatching(window = "square") + scale_region()
  expect_warning(print(p13), "more than 12 regions")
})

test_that("a hatching plot drawn again reuses the region outlines", {
  set.seed(6)
  d <- data.frame(x = runif(400, 0, 100), y = runif(400, 0, 100), imageID = "a")
  d$region <- ifelse(d$x < 50, "region_1", "region_2")
  p <- ggplot2::ggplot(d, ggplot2::aes(x, y, region = region)) + geom_hatching(window = "square")
  lisaClust:::.hatchingCache$clear()
  g1 <- ggplot2::ggplotGrob(p)
  g2 <- ggplot2::ggplotGrob(p)
  expect_length(environment(lisaClust:::.hatchingCache$get)$store, 1)
  panel <- function(g) g$grobs[[grep("panel", g$layout$name)[1]]]
  h1 <- panel(g1)$children[[grep("hatching", names(panel(g1)$children))]]
  h2 <- panel(g2)$children[[grep("hatching", names(panel(g2)$children))]]
  expect_identical(h1$polys, h2$polys)
  # it draws
  pdf(NULL); on.exit(dev.off())
  expect_silent(grid::grid.draw(g1))
})

test_that("regions go back to the right cells when the images are not in order", {
  d <- sim(600)
  d$imageID <- rep(c("i1", "i2", "i3"), 200)  # interleaved
  d$id <- seq_len(nrow(d))
  sorted <- d[order(d$imageID), ]
  set.seed(1)
  a <- lisaClust(d, k = 3, r = c(10, 20))
  set.seed(1)
  b <- lisaClust(sorted, k = 3, r = c(10, 20))
  expect_equal(a[, names(d)], d)
  expect_equal(a$region[order(a$imageID)], b$region)
})
