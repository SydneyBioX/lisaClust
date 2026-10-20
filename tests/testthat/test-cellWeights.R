## cellWeights: 0/1 weights give the curves and regions of the cell types
set.seed(51773)
x <- round(c(runif(200), runif(200) + 1, runif(200) + 2, runif(200) + 3, runif(200) + 3, runif(200) + 2, runif(200) + 1, runif(200)), 4) * 100
y <- round(c(runif(200), runif(200) + 1, runif(200) + 2, runif(200) + 3, runif(200), runif(200) + 1, runif(200) + 2, runif(200) + 3), 4) * 100
cellType <- factor(paste("c", rep(rep(c(1:2), rep(200, 2)), 4), sep = ""))
cells <- data.frame(x, y, cellType, imageID = rep(c("s1", "s2"), c(800, 800)), cellID = paste0("cell", 1:1600))
onehot <- cbind(c1 = as.numeric(cells$cellType == "c1"), c2 = as.numeric(cells$cellType == "c2"))

test_that("lisa() with 0/1 weights equals the cell types", {
  a <- lisa(cells, r = c(10, 20, 50))
  b <- lisa(cells, r = c(10, 20, 50), cellWeights = onehot)
  expect_identical(a, b[rownames(a), colnames(a)])
})

test_that("lisa() with weights does not depend on the order of the cells", {
  w <- onehot * 0.8 + 0.1
  o <- sample(nrow(cells))
  a <- lisa(cells, r = c(10, 20), cellWeights = w)
  b <- lisa(cells[o, ], r = c(10, 20), cellWeights = w[o, ])
  expect_equal(a, b[rownames(a), colnames(a)])
})

test_that("lisaClust() with 0/1 weights gives the same regions", {
  set.seed(1); a <- lisaClust(cells, k = 2)
  set.seed(1); b <- lisaClust(cells, k = 2, cellWeights = onehot)
  expect_identical(a$region, b$region)
})

test_that("a vector of weights equals the matrix with each cell's weight on its own type", {
  q <- runif(nrow(cells), 0.5, 1)
  a <- lisa(cells, r = c(10, 20), cellWeights = q)
  b <- lisa(cells, r = c(10, 20), cellWeights = onehot * q)
  expect_equal(a, b[rownames(a), colnames(a)])
})
