## Inputs as users have them: numeric image IDs, images interleaved, cell IDs repeated in each image, a cell type
## column with another name next to a cellType column, and SingleCellExperiment or SpatialExperiment objects.

simImages <- function(n = 150, images = c(10, 2, 1), seed = 1) {
  set.seed(seed)
  d <- do.call(rbind, lapply(images, function(i) data.frame(
    imageID = i, x = round(runif(n, 0, 300), 1), y = round(runif(n, 0, 300), 1),
    label = sample(c("a", "b", "c"), n, TRUE), cellID = seq_len(n)
  )))
  d$cellType <- sample(c("p", "q"), nrow(d), TRUE)
  d[sample(nrow(d)), ]
}

test_that("lisaClust() gives every cell a region, in its own order, when cell IDs repeat across images", {
  d <- simImages()
  set.seed(1)
  out <- suppressMessages(lisaClust(d, k = 3, r = c(10, 20), cellType = "label"))
  expect_equal(out[, names(d)], d)
  expect_false(any(out$region == ""))
  # the same regions as with unique cell IDs
  u <- d
  u$cellID <- paste0("cell", seq_len(nrow(u)))
  set.seed(1)
  ref <- suppressMessages(lisaClust(u, k = 3, r = c(10, 20), cellType = "label"))
  expect_identical(out$region, ref$region)
})

test_that("lisaClust() keeps the cells of a SpatialExperiment in order (numeric image IDs)", {
  d <- simImages()
  spe <- SpatialExperiment::SpatialExperiment(
    colData = d[, c("imageID", "label", "cellID")], spatialCoords = as.matrix(d[, c("x", "y")])
  )
  colnames(spe) <- paste0("c", seq_len(ncol(spe)))
  out <- suppressMessages(lisaClust(spe, k = 2, r = c(10, 20), cellType = "label"))
  expect_identical(colnames(out), colnames(spe))
  expect_identical(out$imageID, spe$imageID)
  expect_false(any(out$region == ""))
})

test_that("lisa() runs on a SingleCellExperiment and returns the cells in their own order", {
  d <- simImages()
  d$cellID <- paste0("cell", seq_len(nrow(d)))
  sce <- SingleCellExperiment::SingleCellExperiment(colData = d)
  a <- suppressMessages(lisa(sce, r = c(10, 20), cellType = "label"))
  expect_identical(rownames(a), d$cellID)
  # the curves of each cell do not depend on the order of the cells
  sorted <- d[order(d$imageID), ]
  b <- suppressMessages(lisa(sorted, r = c(10, 20), cellType = "label"))
  expect_identical(a, b[rownames(a), ])
})

test_that("lisa() names cells without a cellID column cell_1, cell_2, ... in their own order", {
  d <- simImages(n = 60)
  d$cellID <- NULL
  a <- suppressMessages(lisa(d, r = c(10, 20)))
  expect_identical(rownames(a), paste0("cell_", seq_len(nrow(d))))
})

test_that("lisa() and lisaClust() stop on missing image IDs, cell types or coordinates", {
  d <- simImages(n = 60)
  d$label[3] <- NA
  expect_error(suppressMessages(lisa(d, cellType = "label")), "'label' has missing values")
  expect_error(suppressMessages(lisaClust(d, cellType = "label")), "'label' has missing values")
  d <- simImages(n = 60)
  d$imageID[3] <- NA
  expect_error(suppressMessages(lisa(d)), "'imageID' has missing values")
  d <- simImages(n = 60)
  d$y[3] <- NA
  expect_error(suppressMessages(lisaClust(d)), "'y' has missing values")
})

test_that("image IDs that are factors with unused levels are fine", {
  d <- simImages(n = 100, images = c("a", "b"))
  d$imageID <- factor(d$imageID, levels = c("z", "b", "a"))
  a <- suppressMessages(lisa(d, r = c(10, 20)))
  d$imageID <- as.character(d$imageID)
  expect_identical(a, suppressMessages(lisa(d, r = c(10, 20))))
})

test_that("lisaClust() checks the length of cellWeights", {
  d <- simImages(n = 60)
  expect_error(suppressMessages(lisaClust(d, cellWeights = rep(1, nrow(d) + 1))), "cellWeights")
})

test_that("hatchingPlot() uses the named cell type column and the right regions for a subset of images", {
  d <- simImages(n = 100)
  set.seed(1)
  d <- suppressMessages(lisaClust(d, k = 2, r = c(10, 20), cellType = "label"))
  sce <- SingleCellExperiment::SingleCellExperiment(colData = d)
  p <- expect_no_warning(suppressMessages(hatchingPlot(sce, useImages = "2", cellType = "label")))
  pd <- p$data
  expect_identical(as.character(pd$cellType), d$label[d$imageID == 2])
  expect_identical(as.character(pd$region), d$region[d$imageID == 2])
  pdf(NULL); on.exit(dev.off())
  expect_no_warning(suppressMessages(ggplot2::ggplotGrob(p)))
})

test_that("hatchingPlot() treats numeric regions as categories", {
  d <- simImages(n = 100)
  set.seed(1)
  d <- suppressMessages(lisaClust(d, k = 3, r = c(10, 20)))
  d$regionNumber <- as.integer(factor(d$region))
  p <- suppressMessages(hatchingPlot(d, useImages = c(2, 1), region = "regionNumber"))
  expect_true(is.factor(p$data$region))
  pdf(NULL); on.exit(dev.off())
  expect_no_error(suppressMessages(ggplot2::ggplotGrob(p)))
})
