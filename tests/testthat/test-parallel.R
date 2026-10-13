test_that(".par_lapply() gives the serial results and leaves BLAS threads as they were", {
  skip_on_os("windows")
  X <- setNames(as.list(1:20), letters[1:20])
  blasThreads <- RhpcBLASctl::blas_get_num_procs()
  expect_identical(.par_lapply(X, function(x, y) x * y, y = 3, cores = 2), lapply(X, function(x, y) x * y, y = 3))
  expect_identical(.par_lapply(X, sqrt, cores = 2), lapply(X, sqrt))
  expect_identical(.par_lapply(X, sqrt, cores = BiocParallel::MulticoreParam(2)), lapply(X, sqrt))
  expect_identical(RhpcBLASctl::blas_get_num_procs(), blasThreads)
  expect_error(.par_lapply(1:4, function(x) if (x == 3) stop("bad task") else x, cores = 2), "bad task")
})

test_that("lisa() gives the same curves on two cores as on one", {
  skip_on_os("windows")
  set.seed(51773)
  x <- round(c(runif(200), runif(200) + 1, runif(200) + 2, runif(200) + 3), 4) * 100
  y <- round(c(runif(200), runif(200) + 1, runif(200) + 2, runif(200) + 3), 4) * 100
  cells <- data.frame(
    x = c(x, rev(x)), y = c(y, y), cellType = factor(paste0("c", rep(rep(1:2, c(200, 200)), 4))),
    imageID = rep(c("s1", "s2", "s3", "s4"), each = 400)
  )
  serial <- suppressMessages(lisa(cells, r = c(10, 20, 50)))
  expect_identical(suppressMessages(lisa(cells, r = c(10, 20, 50), cores = 2)), serial)
})
