## regionMap(type = "table"), nameRegions() and regionBoxPlot()

toy <- function() {
  data.frame(imageID = rep(c("a", "b", "c"), each = 4),
             group = rep(c("A", "A", "B"), each = 4),
             cellType = c("t", "t", "i", "i", "t", "t", "t", "i", "i", "i", "i", "t"),
             region = c("r1", "r1", "r2", "r2", "r1", "r1", "r3", "r2", "r2", "r2", "r2", "r3"),
             x = 1:12, y = 1:12)
}

test_that("regionMap(type = 'table') is observed over expected counts", {
  cells <- toy()
  tab <- table(cells$cellType, cells$region)
  expected <- outer(rowSums(tab), colSums(tab)) / sum(tab)
  enrichment <- regionMap(cells, type = "table")
  expect_true(is.matrix(enrichment))
  expect_equal(as.vector(enrichment), as.vector(tab / expected))
  expect_equal(dimnames(enrichment), dimnames(unclass(tab)))
})

test_that("nameRegions names each region by its most enriched marker", {
  cells <- nameRegions(toy(), c(tumour = "t", immune = "i"), regionName = "domain")
  # r1 and r3 hold only t cells, r2 only i cells
  expect_equal(cells$domain, unname(c(r1 = "tumour", r2 = "immune", r3 = "tumour")[cells$region]))
  expect_equal(cells$region, toy()$region)
  expect_equal(nameRegions(toy(), "i")$region, rep("i", 12))
  expect_error(nameRegions(toy(), "b"), "not found")
})

test_that("nameRegions writes to the colData of a SingleCellExperiment", {
  cells <- toy()
  sce <- SingleCellExperiment::SingleCellExperiment(colData = S4Vectors::DataFrame(cells))
  sce <- nameRegions(sce, c(tumour = "t", immune = "i"))
  expect_equal(sce$region, nameRegions(cells, c(tumour = "t", immune = "i"))$region)
})

test_that("regionBoxPlot gives one share per unit and region, with zeros", {
  p <- regionBoxPlot(toy(), condition = "group")
  expect_s3_class(p, "ggplot")
  d <- p$data
  expect_equal(nrow(d), 9)
  expect_equal(d$share[d$unit == "a" & d$region == "r3"], 0)
  expect_equal(d$share[d$unit == "c" & d$region == "r2"], 0.75)
  expect_equal(unique(d$condition[d$unit == "c"]), "B")
  expect_equal(nrow(regionBoxPlot(toy(), "group", regions = "r1")$data), 3)
  bad <- toy(); bad$group[1] <- "B"
  expect_error(regionBoxPlot(bad, "group"), "single")
})
