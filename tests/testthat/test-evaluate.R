library(testthat)

context("Function: findBestNrCluster")

test_that("findBestNrCluster functionallity tests (normal)", {
  dat <- matrix(data = c(1, 2, 3, 4,
                         2.1, 3.9, 6.2, 7.5,
                         0, 1, 3, 2,
                         2, 2, 11, 32,
                         2, 3, 41, 3,
                         2, 0, 1, 0
  ),
  nrow = 6, ncol = 4,
  byrow = TRUE,
  dimnames = list(c(1:6), c(1:4)))

  sim <- MultIS::getSimilarityMatrix(readouts = dat,
                                     self = 1,
                                     upper = TRUE,
                                     method = "rsquared",
                                     parallel = FALSE)

  bnc <- findBestNrCluster(data = dat, sim = sim)
  expect_equal(bnc, 2)
})

test_that("findBestNrCluster functionallity tests (report)", {
  dat <- matrix(data = c(1, 2, 3, 4,
                         2.1, 3.9, 6.2, 7.5,
                         0, 1, 3, 2,
                         2, 2, 11, 32,
                         2, 3, 41, 3,
                         2, 0, 1, 0
  ),
  nrow = 6, ncol = 4,
  byrow = TRUE,
  dimnames = list(c(1:6), c(1:4)))

  sim <- MultIS::getSimilarityMatrix(readouts = dat,
                                     self = 1,
                                     upper = TRUE,
                                     method = "rsquared",
                                     parallel = FALSE)

  expect_output(findBestNrCluster(data = dat, sim = sim, report = TRUE))
})

test_that("findBestNrCluster functionallity tests (returnAll)", {
  dat <- matrix(data = c(1, 2, 3, 4,
                         2.1, 3.9, 6.2, 7.5,
                         0, 1, 3, 2,
                         2, 2, 11, 32,
                         2, 3, 41, 3,
                         2, 0, 1, 0
  ),
  nrow = 6, ncol = 4,
  byrow = TRUE,
  dimnames = list(c(1:6), c(1:4)))

  sim <- MultIS::getSimilarityMatrix(readouts = dat,
                                     self = 1,
                                     upper = TRUE,
                                     method = "rsquared",
                                     parallel = FALSE)

  ev <- findBestNrCluster(data = dat, sim = sim, returnAll = TRUE)

  expect_equal(length(ev), 4)
})

context("Function: evaluateClustering")

test_that("evaluateClustering functionallity tests (silhouette)", {
  dat <- matrix(data = c(1, 2, 3, 4,
                         2.1, 3.9, 6.2, 7.5,
                         0, 1, 3, 2,
                         2, 2, 11, 32,
                         2, 3, 41, 3,
                         2, 0, 1, 0
  ),
  nrow = 6, ncol = 4,
  byrow = TRUE,
  dimnames = list(c(1:6), c(1:4)))

  sim <- MultIS::getSimilarityMatrix(readouts = dat,
                                     self = 1,
                                     upper = TRUE,
                                     method = "rsquared",
                                     parallel = FALSE)

  rec <- reconstruct(readouts = dat, targetCommunities = 3, sim = sim)

  ev <- evaluateClustering(readouts = dat, clustering = rec, sim = sim, method = "silhouette")
  expect_equal(ev, 0.38689220640418287012)
})

test_that("evaluateClustering functionallity tests (sdindex)", {
  dat <- matrix(data = c(1, 2, 3, 4,
                         2.1, 3.9, 6.2, 7.5,
                         0, 1, 3, 2,
                         2, 2, 11, 32,
                         2, 3, 41, 3,
                         2, 0, 1, 0
  ),
  nrow = 6, ncol = 4,
  byrow = TRUE,
  dimnames = list(c(1:6), c(1:4)))

  sim <- MultIS::getSimilarityMatrix(readouts = dat,
                                     self = 1,
                                     upper = TRUE,
                                     method = "rsquared",
                                     parallel = FALSE)

  rec <- reconstruct(readouts = dat, targetCommunities = 3, sim = sim)

  ev <- evaluateClustering(readouts = dat, clustering = rec, sim = sim, method = "sdindex")
  expect_equal(ev, 0.38461041245274385503)
})

test_that("evaluateClustering functionallity tests (ptbiserial)", {
  dat <- matrix(data = c(1, 2, 3, 4,
                         2.1, 3.9, 6.2, 7.5,
                         0, 1, 3, 2,
                         2, 2, 11, 32,
                         2, 3, 41, 3,
                         2, 0, 1, 0
  ),
  nrow = 6, ncol = 4,
  byrow = TRUE,
  dimnames = list(c(1:6), c(1:4)))

  sim <- MultIS::getSimilarityMatrix(readouts = dat,
                                     self = 1,
                                     upper = TRUE,
                                     method = "rsquared",
                                     parallel = FALSE)

  rec <- reconstruct(readouts = dat, targetCommunities = 3, sim = sim)

  ev <- evaluateClustering(readouts = dat, clustering = rec, sim = sim, method = "ptbiserial")
  expect_equal(ev, 0.83990304016878913895)
})

test_that("evaluateClustering functionallity tests (failure cases)", {
  expect_error(evaluateClustering(readouts = NULL, clustering = NULL, sim = NULL, method = "silhouette"))
  expect_error(evaluateClustering(readouts = NULL, clustering = NULL, sim = NULL, method = "foobar"))
})
