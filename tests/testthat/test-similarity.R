library(testthat)

context("Function: getSimilarityMatrix")

test_that("getSimilarityMatrix functionallity tests", {
  dat <- matrix(data = c(1, 2, 3, 4,
                         2.1, 3.9, 6.2, 7.5,
                         0, 1, 3, 2),
                nrow = 3, ncol = 4,
                byrow = TRUE,
                dimnames = list(c(1:3), c(1:4)))

  expect_warning(sim <- MultIS::getSimilarityMatrix(readouts = dat,
                                                    self = 1,
                                                    upper = TRUE,
                                                    method = "rsquared",
                                                    parallel = FALSE))

  expect_true(any(class(sim) == "matrix"))

  expect_equal(nrow(sim), nrow(dat))
  expect_equal(ncol(sim), nrow(dat))

  expect_true(all(diag(sim) == 1))
  expect_equal(colnames(sim), rownames(dat))
  expect_equal(rownames(sim), rownames(dat))

  expect_equal(length(sim), nrow(dat) ** 2)
})

test_that("getSimilarityMatrix functionallity tests (euclidean)", {
  dat <- matrix(data = c(1, 2, 3, 4,
                         2.1, 3.9, 6.2, 7.5,
                         0, 1, 3, 2),
                nrow = 3, ncol = 4,
                byrow = TRUE,
                dimnames = list(c(1:3), c(1:4)))

  sim <- MultIS::getSimilarityMatrix(readouts = dat,
                                     self = 1,
                                     upper = TRUE,
                                     method = "euclidean",
                                     parallel = FALSE)

  expect_true(any(class(sim) == "matrix"))

  expect_equal(nrow(sim), nrow(dat))
  expect_equal(ncol(sim), nrow(dat))

  expect_true(all(diag(sim) == 1))
  expect_equal(colnames(sim), rownames(dat))
  expect_equal(rownames(sim), rownames(dat))

  expect_equal(length(sim), nrow(dat) ** 2)
})

test_that("getSimilarityMatrix functionallity tests (maximum)", {
  dat <- matrix(data = c(1, 2, 3, 4,
                         2.1, 3.9, 6.2, 7.5,
                         0, 1, 3, 2),
                nrow = 3, ncol = 4,
                byrow = TRUE,
                dimnames = list(c(1:3), c(1:4)))

  sim <- MultIS::getSimilarityMatrix(readouts = dat,
                                     self = 1,
                                     upper = TRUE,
                                     method = "maximum",
                                     parallel = FALSE)

  expect_true(any(class(sim) == "matrix"))

  expect_equal(nrow(sim), nrow(dat))
  expect_equal(ncol(sim), nrow(dat))

  expect_true(all(diag(sim) == 1))
  expect_equal(colnames(sim), rownames(dat))
  expect_equal(rownames(sim), rownames(dat))

  expect_equal(length(sim), nrow(dat) ** 2)
})

test_that("getSimilarityMatrix functionallity tests (manhattan)", {
  dat <- matrix(data = c(1, 2, 3, 4,
                         2.1, 3.9, 6.2, 7.5,
                         0, 1, 3, 2),
                nrow = 3, ncol = 4,
                byrow = TRUE,
                dimnames = list(c(1:3), c(1:4)))

  sim <- MultIS::getSimilarityMatrix(readouts = dat,
                                     self = 1,
                                     upper = TRUE,
                                     method = "manhattan",
                                     parallel = FALSE)

  expect_true(any(class(sim) == "matrix"))

  expect_equal(nrow(sim), nrow(dat))
  expect_equal(ncol(sim), nrow(dat))

  expect_true(all(diag(sim) == 1))
  expect_equal(colnames(sim), rownames(dat))
  expect_equal(rownames(sim), rownames(dat))

  expect_equal(length(sim), nrow(dat) ** 2)
})

test_that("getSimilarityMatrix functionallity tests (canberra)", {
  dat <- matrix(data = c(1, 2, 3, 4,
                         2.1, 3.9, 6.2, 7.5,
                         0, 1, 3, 2),
                nrow = 3, ncol = 4,
                byrow = TRUE,
                dimnames = list(c(1:3), c(1:4)))

  sim <- MultIS::getSimilarityMatrix(readouts = dat,
                                     self = 1,
                                     upper = TRUE,
                                     method = "canberra",
                                     parallel = FALSE)

  expect_true(any(class(sim) == "matrix"))

  expect_equal(nrow(sim), nrow(dat))
  expect_equal(ncol(sim), nrow(dat))

  expect_true(all(diag(sim) == 1))
  expect_equal(colnames(sim), rownames(dat))
  expect_equal(rownames(sim), rownames(dat))

  expect_equal(length(sim), nrow(dat) ** 2)
})

test_that("getSimilarityMatrix functionallity tests (binary)", {
  dat <- matrix(data = c(1, 2, 3, 4,
                         2.1, 3.9, 6.2, 7.5,
                         0, 1, 3, 2),
                nrow = 3, ncol = 4,
                byrow = TRUE,
                dimnames = list(c(1:3), c(1:4)))

  sim <- MultIS::getSimilarityMatrix(readouts = dat,
                                     self = 1,
                                     upper = TRUE,
                                     method = "binary",
                                     parallel = FALSE)

  expect_true(any(class(sim) == "matrix"))

  expect_equal(nrow(sim), nrow(dat))
  expect_equal(ncol(sim), nrow(dat))

  expect_true(all(diag(sim) == 1))
  expect_equal(colnames(sim), rownames(dat))
  expect_equal(rownames(sim), rownames(dat))

  expect_equal(length(sim), nrow(dat) ** 2)
})

test_that("getSimilarityMatrix functionallity tests (minkowski)", {
  dat <- matrix(data = c(1, 2, 3, 4,
                         2.1, 3.9, 6.2, 7.5,
                         0, 1, 3, 2),
                nrow = 3, ncol = 4,
                byrow = TRUE,
                dimnames = list(c(1:3), c(1:4)))

  sim <- MultIS::getSimilarityMatrix(readouts = dat,
                                     self = 1,
                                     upper = TRUE,
                                     method = "minkowski",
                                     parallel = FALSE)

  expect_true(any(class(sim) == "matrix"))

  expect_equal(nrow(sim), nrow(dat))
  expect_equal(ncol(sim), nrow(dat))

  expect_true(all(diag(sim) == 1))
  expect_equal(colnames(sim), rownames(dat))
  expect_equal(rownames(sim), rownames(dat))

  expect_equal(length(sim), nrow(dat) ** 2)
})

test_that("getSimilarityMatrix functionallity tests (rsquared)", {
  dat <- matrix(data = c(1, 2, 3, 4,
                         2.1, 3.9, 6.2, 7.5,
                         0, 1, 3, 2),
                nrow = 3, ncol = 4,
                byrow = TRUE,
                dimnames = list(c(1:3), c(1:4)))

  expect_warning(sim <- MultIS::getSimilarityMatrix(readouts = dat,
                                                    self = 1,
                                                    upper = TRUE,
                                                    method = "rsquared",
                                                    parallel = FALSE))

  expect_true(any(class(sim) == "matrix"))

  expect_equal(nrow(sim), nrow(dat))
  expect_equal(ncol(sim), nrow(dat))

  expect_true(all(diag(sim) == 1))
  expect_equal(colnames(sim), rownames(dat))
  expect_equal(rownames(sim), rownames(dat))

  expect_equal(length(sim), nrow(dat) ** 2)
})

test_that("getSimilarityMatrix functionallity tests (no ISs)", {
  dat <- matrix(data = 0,
                nrow = 0, ncol = 3,
                byrow = TRUE,
                dimnames = list(NULL, c(1:3)))

  expect_warning(sim <- MultIS::getSimilarityMatrix(readouts = dat,
                                                    self = 1,
                                                    upper = TRUE,
                                                    method = "rsquared",
                                                    parallel = FALSE))

  expect_true(inherits(x = sim, "matrix"))

  expect_equal(nrow(sim), nrow(dat))
  expect_equal(ncol(sim), nrow(dat))

  expect_true(all(diag(sim) == 1))
  expect_equal(colnames(sim), rownames(dat))
  expect_equal(rownames(sim), rownames(dat))

  expect_equal(length(sim), nrow(dat) ** 2)
})

test_that("getSimilarityMatrix functionallity tests (readouts need to be a matrix and have names)", {
  dat <- data.frame("A" = c(1:3),
                    "B" = c(4:6),
                    "C" = c(7:9))
  expect_error(MultIS::getSimilarityMatrix(readouts = dat,
                                           self = 1,
                                           upper = TRUE,
                                           method = "rsquared",
                                           parallel = FALSE))

  dat <- matrix(data = c(1:9),
                nrow = 3, ncol = 3,
                byrow = TRUE,
                dimnames = list(NULL, c(1:3)))
  expect_error(MultIS::getSimilarityMatrix(readouts = dat,
                                           self = 1,
                                           upper = TRUE,
                                           method = "rsquared",
                                           parallel = FALSE))

  dat <- matrix(data = c(1:9),
                nrow = 3, ncol = 3,
                byrow = TRUE,
                dimnames = list(c(1:3), NULL))
  expect_error(MultIS::getSimilarityMatrix(readouts = dat,
                                           self = 1,
                                           upper = TRUE,
                                           method = "rsquared",
                                           parallel = FALSE))
})

test_that("getSimilarityMatrix functionallity tests (method must be given as a string or function)", {
  dat <- matrix(data = c(1:9),
                nrow = 3, ncol = 3,
                byrow = TRUE,
                dimnames = list(c(1:3), c(1:3)))
  expect_error(MultIS::getSimilarityMatrix(readouts = dat,
                                           self = 1,
                                           upper = TRUE,
                                           method = data.frame(),
                                           parallel = FALSE))
})
