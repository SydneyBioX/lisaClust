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

test_that("hatching lines cross each outline an even number of times, with their midpoints inside", {
  # a region on a grid, like those regionPoly() makes: outlines with many vertices on grid lines
  set.seed(4)
  mat <- matrix(runif(400) < 0.5, 20, 20)
  W <- spatstat.geom::as.polygonal(spatstat.geom::owin(c(0, 20), c(0, 20), mask = mat))
  edges <- lapply(W$bdry, function(b) {
    n <- length(b$x)
    list(x1 = b$x, y1 = b$y, x2 = b$x[c(seq_len(n)[-1], 1)], y2 = b$y[c(seq_len(n)[-1], 1)])
  })
  vx <- unlist(lapply(W$bdry, `[[`, "x")); vy <- unlist(lapply(W$bdry, `[[`, "y"))
  lines <- list(rbind(c(-1, 5), c(21, 5)),                 # along grid lines: through vertices and along edges
                rbind(c(vx[3], -1), c(vx[3], 21)),
                rbind(c(vx[7] - 30, vy[7] - 30), c(vx[7] + 30, vy[7] + 30)),  # diagonal through a vertex
                rbind(c(-1, 3.3), c(21, 7.7)))              # general position
  for (h in lines) {
    cr <- lisaClust:::hatchCrossings(edges, h[1, ], h[2, ])
    expect_true(nrow(cr) %% 2 == 0)
    if (nrow(cr) >= 2) {
      mid <- (cr[c(TRUE, FALSE), , drop = FALSE] + cr[c(FALSE, TRUE), , drop = FALSE]) / 2
      long <- sqrt(rowSums((cr[c(TRUE, FALSE), , drop = FALSE] - cr[c(FALSE, TRUE), , drop = FALSE])^2)) > 1e-9
      expect_true(all(spatstat.geom::inside.owin(mid[long, 1], mid[long, 2], W)))
    }
  }
})

test_that("hatching crossings agree with the line intersection used before, in general position", {
  # the scalar version lisaClust used before 1.21.1
  lineIntersection <- function(P1, P2, P3, P4) {
    P1 <- round(as.vector(P1), 10); P2 <- round(as.vector(P2), 10); P3 <- round(as.vector(P3), 10); P4 <- round(as.vector(P4), 10)
    dx1 <- P1[1] - P2[1]; dx2 <- P3[1] - P4[1]; dy1 <- P1[2] - P2[2]; dy2 <- P3[2] - P4[2]
    D <- det(rbind(c(dx1, dy1), c(dx2, dy2)))
    if (is.na(D) | D == 0) return(c(Inf, Inf))
    D1 <- det(rbind(P1, P2)); D2 <- det(rbind(P3, P4))
    X <- round(det(rbind(c(D1, dx1), c(D2, dx2))) / D, 10); Y <- round(det(rbind(c(D1, dy1), c(D2, dy2))) / D, 10)
    l1 <- -((X - P1[1]) * dx1 + (Y - P1[2]) * dy1) / (dx1^2 + dy1^2)
    l2 <- -((X - P3[1]) * dx2 + (Y - P3[2]) * dy2) / (dx2^2 + dy2^2)
    if (!((l1 >= 0) & (l1 <= 1) & (l2 >= 0) & (l2 <= 1))) return(c(NA, NA))
    c(X, Y)
  }
  set.seed(5)
  th <- sort(runif(30, 0, 2 * pi)); rr <- runif(30, 5, 10)
  ring <- list(x1 = rr * cos(th), y1 = rr * sin(th))
  ring$x2 <- ring$x1[c(2:30, 1)]; ring$y2 <- ring$y1[c(2:30, 1)]
  h <- rbind(c(-12, -3.1), c(12, 4.7))
  got <- lisaClust:::hatchCrossings(list(ring), h[1, ], h[2, ])
  old <- t(vapply(1:30, function(i) lineIntersection(c(ring$x1[i], ring$y1[i]), c(ring$x2[i], ring$y2[i]), h[1, ], h[2, ]),
                  numeric(2)))
  old <- old[!is.na(old[, 1]), , drop = FALSE]
  old <- old[order(old[, 1]), , drop = FALSE]
  expect_equal(unname(got), unname(old), tolerance = 1e-9)
})

test_that("a hatching plot drawn again comes from the cache", {
  set.seed(6)
  d <- data.frame(x = runif(400, 0, 100), y = runif(400, 0, 100), imageID = "a")
  d$cellType <- ifelse(d$x < 50, "l", "r")
  d$region <- ifelse(d$x < 50, "region_1", "region_2")
  p <- ggplot2::ggplot(d, ggplot2::aes(x, y, region = region)) + geom_hatching(window = "square", nbp = 50)
  lisaClust:::.hatchingCache$clear()
  t1 <- system.time(g1 <- ggplot2::ggplotGrob(p))[3]
  t2 <- system.time(g2 <- ggplot2::ggplotGrob(p))[3]
  panel <- function(g) g$grobs[[grep("panel", g$layout$name)[1]]]
  h1 <- panel(g1)$children[[grep("hatching", names(panel(g1)$children))]]
  h2 <- panel(g2)$children[[grep("hatching", names(panel(g2)$children))]]
  expect_identical(h1, h2)
  expect_length(environment(lisaClust:::.hatchingCache$get)$store, 1)
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
