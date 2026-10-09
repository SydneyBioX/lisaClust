## The C++ core against plain-R references written from the definitions.

sim <- function(n = 300, seed = 1) {
  set.seed(seed)
  data.frame(x = round(runif(n, 0, 200), 1), y = round(runif(n, 0, 150), 1),
             cellType = factor(sample(c("a", "b", "c"), n, TRUE, prob = c(0.5, 0.3, 0.2))))
}

# The local indicators from a distance matrix: weighted counts per distance band, accumulated over the bands
# with any pair, and the expectation with the cell's edge share where it has a neighbour of that type in the band.
refCurves <- function(d, Rs, wt, lam, edge, L = FALSE) {
  n <- nrow(d); K <- nlevels(d$cellType); nb <- length(Rs) - 1
  D <- as.matrix(dist(d[, c("x", "y")])); diag(D) <- Inf
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
  for (seed in 1:3) for (L in c(FALSE, TRUE)) {
    d <- sim(seed = seed)
    Rs <- c(0, 10, 25, 40)
    wt <- runif(nrow(d), 0.5, 1.5); lam <- as.numeric(table(d$cellType)) / (200 * 150)
    edge <- matrix(runif(nrow(d) * 3, 0.5, 1), nrow(d))
    got <- lisaClust:::.localCurves(d$x, d$y, as.integer(d$cellType), 3L, Rs, Rs[-1], wt, lam, edge, L)
    ref <- refCurves(d, Rs, wt, lam, edge, L)
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

test_that("hatching intersections match the scalar line intersection", {
  # the scalar version lisaClust used before
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
  set.seed(4)
  seg <- matrix(round(runif(400, 0, 100), sample(0:3, 400, TRUE)), ncol = 4)
  seg[1:20, 3] <- seg[1:20, 1]  # vertical edges
  seg[21:40, 4] <- seg[21:40, 2]  # horizontal edges
  for (h in list(rbind(c(10, 0), c(60, 100)), rbind(c(0, 37), c(100, 37)), rbind(c(42, 0), c(42, 100)))) {
    got <- lisaClust:::segmentIntersections(seg[, 1], seg[, 2], seg[, 3], seg[, 4], h[1, ], h[2, ])
    ref <- t(apply(seg, 1, function(s) lineIntersection(s[1:2], s[3:4], h[1, ], h[2, ])))
    expect_identical(unname(got), unname(ref))
  }
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
