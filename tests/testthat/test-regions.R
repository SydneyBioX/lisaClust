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

test_that("regionBoxPlot(subject =) pools each patient's images, one point per patient", {
  cells <- toy(); cells$patient <- c(rep("p1", 8), rep("p2", 4))   # images a and b are patient p1's
  p <- regionBoxPlot(cells, condition = "group", subject = "patient")
  expect_identical(p$data, regionBoxPlot(cells, condition = "group", imageID = "patient")$data)
  expect_setequal(unique(p$data$unit), c("p1", "p2"))
  expect_equal(p$data$share[p$data$unit == "p1" & p$data$region == "r1"], 4 / 8)
})

test_that("regionMap() turns the region labels upright only when there are many regions", {
  few <- regionMap(toy())
  expect_null(few$theme$axis.text.x$angle)
  set.seed(1)
  many <- data.frame(cellType = sample(c("t", "i", "s"), 400, TRUE, prob = c(5, 3, 1)),
                     region = paste0("niche ", sample(1:8, 400, TRUE)))
  expect_equal(regionMap(many)$theme$axis.text.x$angle, 90)
})
