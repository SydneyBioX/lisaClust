test_that("the C++ concave hull gives concaveman's polygons exactly", {
  # polygons from concaveman::concaveman() 1.2.0, for the points makeWindow() builds around 150 cells
  ref <- readRDS(test_path("concaveman_reference.rds"))
  for (case in names(ref)) {
    r <- ref[[case]]
    expect_identical(.concaveHull(r$points[, 1], r$points[, 2], r$concavity, r$length_threshold), r$polygon,
                     label = case)
  }
})

test_that("makeWindow() makes a concave window", {
  ref <- readRDS(test_path("concaveman_reference.rds"))
  cells <- data.frame(x = ref$uniform$points[seq(1, nrow(ref$uniform$points), 9), 1],
                      y = ref$uniform$points[seq(1, nrow(ref$uniform$points), 9), 2])
  ow <- suppressMessages(makeWindow(cells, "concave"))
  expect_s3_class(ow, "owin")
  expect_true(all(spatstat.geom::inside.owin(cells$x, cells$y, ow)))
})
