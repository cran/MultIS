library(testthat)

context("File: filterFunctions.R")

test_that("convert_columnwise_percent functionallity tests", {
  dat <- matrix(data = c(0, 0, 15, 5,  3,
                         2, 0, 1,  5,  3,
                         3, 0, 1,  3,  0,
                         0, 0, 1,  10, 2),
                nrow = 4, ncol = 5, byrow = TRUE,
                dimnames = list(c("1", "2", "3", "4"),
                                c("10", "20", "30", "40", "45")))
  class(dat) <- c("matrix", "timeseries")

  n <- convert_columnwise_percent(dat)

  expect_identical(class(dat), class(n))
  expect_equal(nrow(n), nrow(dat))
  expect_equal(ncol(n), ncol(dat))
  expect_identical(colnames(n), colnames(dat))
  expect_identical(rownames(n), rownames(dat))

  expect_true(all(colSums(n) == 1 | colSums(n) == 0))
})

test_that("filter_atTP_min functionallity tests", {
  dat <- matrix(data = c(0, 0, 15, 5,  3,
                         2, 0, 1,  5,  3,
                         3, 0, 1,  3,  0,
                         0, 0, 1,  10, 2),
                nrow = 4, ncol = 5, byrow = TRUE,
                dimnames = list(c("1", "2", "3", "4"),
                                c("10", "20", "30", "40", "45")))
  class(dat) <- c("matrix", "timeseries")

  AT <- "45"
  MIN <- 2

  n <- filter_atTP_min(dat = dat, at = AT, min = MIN)

  expect_identical(class(dat), class(n))
  expect_equal(nrow(n), 3)
  expect_equal(ncol(n), ncol(dat))
  expect_identical(colnames(n), colnames(dat))
  expect_true(all(rownames(n) %in% rownames(dat)))

  expect_true(all(n[,AT] >= MIN))

  # Substring match
  AT <- "4"
  MIN <- 10

  n <- filter_atTP_min(dat = dat, at = AT, min = MIN)

  expect_identical(class(dat), class(n))
  expect_equal(nrow(n), 1)
  expect_equal(ncol(n), ncol(dat))
  expect_identical(colnames(n), colnames(dat))
  expect_true(all(rownames(n) %in% rownames(dat)))

  expect_true(all(rowSums(n[,grep(pattern = AT, x = colnames(n)), drop = FALSE], na.rm = TRUE) >= MIN))

  # No ISs
  AT <- "45"
  MIN <- 5

  n <- filter_atTP_min(dat = dat, at = AT, min = MIN)

  expect_identical(class(dat), class(n))
  expect_equal(nrow(n), 0)
  expect_equal(ncol(n), ncol(dat))
  expect_identical(colnames(n), colnames(dat))
  expect_true(all(rownames(n) %in% rownames(dat)))

  expect_true(all(n[,AT] >= MIN))
})

test_that("filter_atTP_biggestN functionallity tests", {
  dat <- matrix(data = c(0, 0, 15, 5,  3,
                         2, 0, 1,  5,  3,
                         3, 0, 1,  3,  0,
                         0, 0, 1,  10, 2),
                nrow = 4, ncol = 5, byrow = TRUE,
                dimnames = list(c("1", "2", "3", "4"),
                                c("10", "20", "30", "40", "45")))
  class(dat) <- c("matrix", "timeseries")

  AT <- "45"
  N <- 2L
  n <- filter_atTP_biggestN(dat = dat, at = AT, n = N)

  expect_identical(class(dat), class(n))
  expect_equal(nrow(n), N)
  expect_equal(ncol(n), ncol(dat))
  expect_identical(colnames(n), colnames(dat))
  expect_true(all(rownames(n) %in% rownames(dat)))

  # Substring match
  AT <- "4"
  N <- 2L
  n <- filter_atTP_biggestN(dat = dat, at = AT, n = N)

  expect_identical(class(dat), class(n))
  expect_equal(nrow(n), N + 1) # A tie between "1" and "2"
  expect_equal(ncol(n), ncol(dat))
  expect_identical(colnames(n), colnames(dat))
  expect_true(all(rownames(n) %in% rownames(dat)))

  # No ISs
  AT <- "45"
  N <- 0L
  n <- filter_atTP_biggestN(dat = dat, at = AT, n = N)

  expect_identical(class(dat), class(n))
  expect_equal(nrow(n), N)
  expect_equal(ncol(n), ncol(dat))
  expect_identical(colnames(n), colnames(dat))
  expect_true(all(rownames(n) %in% rownames(dat)))
})

test_that("filter_nrTP_min functionallity tests", {
  dat <- matrix(data = c(0, 0, 15, 0,  1, # 2
                         2, 0, 0,  0,  0, # 1
                         3, 2, 1,  NA,  NA, # 3
                         1, NA, 1,  10, 2), # 4
                nrow = 4, ncol = 5, byrow = TRUE,
                dimnames = list(c("1", "2", "3", "4"),
                                c("10", "20", "30", "40", "45")))
  class(dat) <- c("matrix", "timeseries")

  MIN <- 2L
  n <- filter_nrTP_min(dat = dat, min = MIN)

  expect_identical(class(dat), class(n))
  expect_equal(nrow(n), 3)
  expect_equal(ncol(n), ncol(dat))
  expect_identical(colnames(n), colnames(dat))
  expect_true(all(rownames(n) %in% rownames(dat)))

  # Filter out NA
  MIN <- 4L
  n <- filter_nrTP_min(dat = dat, min = MIN)

  expect_identical(class(dat), class(n))
  expect_equal(nrow(n), 1)
  expect_equal(ncol(n), ncol(dat))
  expect_identical(colnames(n), colnames(dat))
  expect_true(all(rownames(n) %in% rownames(dat)))

  # No ISs
  MIN <- 5L
  n <- filter_nrTP_min(dat = dat, min = MIN)

  expect_identical(class(dat), class(n))
  expect_equal(nrow(n), 0)
  expect_equal(ncol(n), ncol(dat))
  expect_identical(colnames(n), colnames(dat))
  expect_true(all(rownames(n) %in% rownames(dat)))
})

test_that("filter_names and filter_ISnames functionallity tests", {
  dat <- matrix(data = NA,
                nrow = 8, ncol = 5, byrow = TRUE,
                dimnames = list(c(as.character(1:8)),
                                c("10", "20", "30", "40", "45")))
  class(dat) <- c("matrix", "timeseries")

  n <- filter_ISnames(dat = dat, by = ".")

  expect_identical(class(dat), class(n))
  expect_equal(nrow(n), nrow(dat))
  expect_equal(ncol(n), ncol(dat))
  expect_identical(colnames(n), colnames(dat))
  expect_identical(rownames(n), rownames(dat))

  # Actually split, shortest unique prefix: is 8 characters long
  rownames(dat) <- c("ABBC_12dw(-)_213",
                     "ABBC_12sw(-)_213",
                     "ABBC_22sw(-)_21f(+)",
                     "AGBC_12sw(-)_313",
                     "AFBC_12sw(-)_213",
                     "AHBC_12sw(-)_213",
                     "AB3C_12sw(-)_213",
                     "ABB3_12sw(-)_213")
  n <- filter_ISnames(dat = dat, by = ".")

  expect_identical(class(dat), class(n))
  expect_equal(nrow(n), nrow(dat))
  expect_equal(ncol(n), ncol(dat))
  expect_identical(colnames(n), colnames(dat))
  expect_true(all(startsWith(rownames(dat), rownames(n))))

  # Actually split, shortest unique prefix: is up to the +/-
  rownames(dat) <- c("ABBC_12dw(-)_213",
                     "ABBC_12sw(-)_213",
                     "ABBC_22sw(-)_21f(+)",
                     "AGBC_12sw(-)_313",
                     "AGBC_12sw(+)_213",
                     "AHBC_12sw(-)_213",
                     "AB3C_12sw(-)_213",
                     "ABB3_12sw(-)_213")
  n <- filter_ISnames(dat = dat)

  expect_identical(class(dat), class(n))
  expect_equal(nrow(n), nrow(dat))
  expect_equal(ncol(n), ncol(dat))
  expect_identical(colnames(n), colnames(dat))
  expect_true(all(startsWith(rownames(dat), rownames(n))))

  # Actually split, a name is double, so there exits no shortest unique prefix
  rownames(dat) <- c("ABBC_12dw(-)_213",
                     "ABBC_12sw(-)_213",
                     "ABBC_22sw(-)_21f(+)",
                     "AGBC_12sw(-)_313",
                     "AGBC_12sw(-)_313",
                     "AHBC_12sw(-)_213",
                     "AB3C_12sw(-)_213",
                     "ABB3_12sw(-)_213")
  n <- filter_ISnames(dat = dat)

  expect_identical(class(dat), class(n))
  expect_equal(nrow(n), nrow(dat))
  expect_equal(ncol(n), ncol(dat))
  expect_identical(colnames(n), colnames(dat))
  expect_true(all(startsWith(rownames(dat), rownames(n))))
})

test_that("filter_match functionallity tests", {
  dat <- matrix(data = c(2, 4, 0, 3, 4, 5,
                         1, 1, 0, 3, 1, 2,
                         9, 3, 0, 2, 1, 0,
                         3, 9, 0, 1, 1, 9),
                nrow = 4, ncol = 6, byrow = TRUE,
                dimnames = list(c("1", "2", "3", "4"),
                                c("A_1", "A_2", "B_1", "B_2", "B_3", "C_1")))
  class(dat) <- c("matrix", "timeseries")

  n <- filter_match(dat = dat, match = "1")

  expect_identical(class(dat), class(n))
  expect_equal(nrow(n), nrow(dat))
  expect_equal(ncol(n), length(grep("1", colnames(dat))))
  expect_true(all(colnames(n) %in% colnames(dat)))
  expect_equal(rownames(n), rownames(dat))


  # Filter all rows
  dat <- matrix(data = c(2, 4, 0, 3, 4, 5,
                         1, 1, 0, 3, 1, 2,
                         9, 3, 0, 2, 1, 0,
                         3, 9, 0, 1, 1, 9),
                nrow = 4, ncol = 6, byrow = TRUE,
                dimnames = list(c("1", "2", "3", "4"),
                                c("A_1", "A_2", "B_1", "B_2", "B_3", "C_1")))
  class(dat) <- c("matrix", "timeseries")

  n <- filter_match(dat = dat, match = "foobar")

  expect_identical(class(dat), class(n))
  expect_equal(nrow(n), nrow(dat))
  expect_equal(ncol(n), 0)
  expect_true(all(colnames(n) %in% colnames(dat)))
  expect_equal(rownames(n), rownames(dat))
  expect_null(colnames(n))

  # Filter all rows
  dat <- matrix(data = c(2, 4, 0, 3, 4, 5,
                         1, 1, 0, 3, 1, 2,
                         9, 3, 0, 2, 1, 0,
                         3, 9, 0, 1, 1, 9),
                nrow = 4, ncol = 6, byrow = TRUE,
                dimnames = list(c("1", "2", "3", "4"),
                                c("A_1", "A_2", "B_1", "B_2", "B_3", "C_1")))
  class(dat) <- c("matrix", "timeseries")

  n <- filter_match(dat = dat, match = "")

  expect_identical(class(dat), class(n))
  expect_equal(nrow(n), nrow(dat))
  expect_equal(ncol(n), ncol(dat))
  expect_equal(colnames(n), colnames(dat))
  expect_equal(rownames(n), rownames(dat))
})

test_that("filter_measurement_names functionallity tests", {
  dat <- matrix(data = NA,
                nrow = 4, ncol = 5, byrow = TRUE,
                dimnames = list(c("1", "2", "3", "4"),
                                c("A_1_Z_999", "B_2_Y_998", "C_3_X_997", "D_4_W_996", "E_5_V_995")))
  class(dat) <- c("matrix", "timeseries")

  n <- filter_measurement_names(dat = dat, elems = c(1), by = "_")

  expect_identical(class(dat), class(n))
  expect_equal(nrow(n), nrow(dat))
  expect_equal(ncol(n), ncol(dat))
  for (i in 1:length(colnames(n)))
    expect_true(grepl(colnames(n)[i], colnames(dat)[i]))

  n <- filter_measurement_names(dat = dat, elems = c(1, 2), by = "_")

  expect_identical(class(dat), class(n))
  expect_equal(nrow(n), nrow(dat))
  expect_equal(ncol(n), ncol(dat))
  for (i in 1:length(colnames(n)))
    expect_true(grepl(colnames(n)[i], colnames(dat)[i]))

  n <- filter_measurement_names(dat = dat, elems = c(1, 3), by = "_")

  expect_identical(class(dat), class(n))
  expect_equal(nrow(n), nrow(dat))
  expect_equal(ncol(n), ncol(dat))
  for (i in 1:length(colnames(n)))
    expect_false(grepl(colnames(n)[i], colnames(dat)[i]))

  n <- filter_measurement_names(dat = dat, elems = c(4, 2, 1), by = "_")

  expect_identical(class(dat), class(n))
  expect_equal(nrow(n), nrow(dat))
  expect_equal(ncol(n), ncol(dat))
  for (i in 1:length(colnames(n)))
    expect_false(grepl(colnames(n)[i], colnames(dat)[i]))
})

test_that("filter_combine_measurement functionallity tests", {
  dat <- matrix(data = c(2, 4, 1, 3, 4, 5,
                         1, 1, 1, 1, 1, 1,
                         9, 3, 8, 2, 1, 0,
                         3, 9, 8, 1, 1, 9),
                nrow = 4, ncol = 6, byrow = TRUE,
                dimnames = list(c("1", "2", "3", "4"),
                                c("A_1", "A_2", "B_1", "B_2", "B_3", "C_1")))
  class(dat) <- c("matrix", "timeseries")

  n <- filter_measurement_names(dat = dat, elems = c(1), by = "_")
  n <- filter_combine_measurements(dat = n, pre.norm = FALSE, post.norm = FALSE)

  expect_identical(class(dat), class(n))
  expect_equal(nrow(n), nrow(dat))
  expect_equal(ncol(n), 3)

  n <- filter_measurement_names(dat = dat, elems = c(1), by = "_")
  n <- filter_combine_measurements(dat = n, pre.norm = TRUE, post.norm = FALSE)

  expect_identical(class(dat), class(n))
  expect_equal(nrow(n), nrow(dat))
  expect_equal(ncol(n), 3)
  expect_equal(colSums(n), c("A" = 2, "B" = 3, "C" = 1))

  n <- filter_measurement_names(dat = dat, elems = c(1), by = "_")
  n <- filter_combine_measurements(dat = n, pre.norm = FALSE, post.norm = TRUE)

  expect_identical(class(dat), class(n))
  expect_equal(nrow(n), nrow(dat))
  expect_equal(ncol(n), 3)
  expect_equal(colSums(n), c("A" = 1, "B" = 1, "C" = 1))
})

test_that("filter_zero_rows functionallity tests", {
  dat <- matrix(data = c(2, 4, 1, 3, 4, 5,
                         0, 0, 0, 0, 0, 0,
                         9, 3, 8, 2, 1, 0,
                         3, 9, 8, 1, 1, 9),
                nrow = 4, ncol = 6, byrow = TRUE,
                dimnames = list(c("1", "2", "3", "4"),
                                c("A_1", "A_2", "B_1", "B_2", "B_3", "C_1")))
  class(dat) <- c("matrix", "timeseries")

  n <- filter_zero_rows(dat = dat)

  expect_identical(class(dat), class(n))
  expect_equal(nrow(n), nrow(dat) - 1)
  expect_equal(ncol(n), ncol(dat))
  expect_equal(colnames(n), colnames(dat))
  expect_true(all(rownames(n) %in% rownames(dat)))

  # Filter all rows
  dat <- matrix(data = c(0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 0),
                nrow = 4, ncol = 6, byrow = TRUE,
                dimnames = list(c("1", "2", "3", "4"),
                                c("A_1", "A_2", "B_1", "B_2", "B_3", "C_1")))
  class(dat) <- c("matrix", "timeseries")

  n <- filter_zero_rows(dat = dat)

  expect_identical(class(dat), class(n))
  expect_equal(nrow(n), 0)
  expect_equal(ncol(n), ncol(dat))
  expect_equal(colnames(n), colnames(dat))
  expect_true(all(rownames(n) %in% rownames(dat)))
  expect_null(rownames(n))

  # Filter all rows
  dat <- matrix(data = c(2, 4, 1, 3, 4, 5,
                         0, 0, 0, 0, 0, 0,
                         0, 0, NA, 0, NA, NA,
                         3, 9, 8, 1, 1, 9),
                nrow = 4, ncol = 6, byrow = TRUE,
                dimnames = list(c("1", "2", "3", "4"),
                                c("A_1", "A_2", "B_1", "B_2", "B_3", "C_1")))
  class(dat) <- c("matrix", "timeseries")

  n <- filter_zero_rows(dat = dat)

  expect_identical(class(dat), class(n))
  expect_equal(nrow(n), 2)
  expect_equal(ncol(n), ncol(dat))
  expect_equal(colnames(n), colnames(dat))
  expect_true(all(rownames(n) %in% rownames(dat)))
})

test_that("filter_zero_columns functionallity tests", {
  dat <- matrix(data = c(2, 4, 0, 3, 4, 5,
                         1, 1, 0, 3, 1, 2,
                         9, 3, 0, 2, 1, 0,
                         3, 9, 0, 1, 1, 9),
                nrow = 4, ncol = 6, byrow = TRUE,
                dimnames = list(c("1", "2", "3", "4"),
                                c("A_1", "A_2", "B_1", "B_2", "B_3", "C_1")))
  class(dat) <- c("matrix", "timeseries")

  n <- filter_zero_columns(dat = dat)

  expect_identical(class(dat), class(n))
  expect_equal(nrow(n), nrow(dat))
  expect_equal(ncol(n), ncol(dat) - 1)
  expect_true(all(colnames(n) %in% colnames(dat)))
  expect_equal(rownames(n), rownames(dat))


  # Filter all rows
  dat <- matrix(data = c(0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 0),
                nrow = 4, ncol = 6, byrow = TRUE,
                dimnames = list(c("1", "2", "3", "4"),
                                c("A_1", "A_2", "B_1", "B_2", "B_3", "C_1")))
  class(dat) <- c("matrix", "timeseries")

  n <- filter_zero_columns(dat = dat)

  expect_identical(class(dat), class(n))
  expect_equal(nrow(n), nrow(dat))
  expect_equal(ncol(n), 0)
  expect_true(all(colnames(n) %in% colnames(dat)))
  expect_equal(rownames(n), rownames(dat))
  expect_null(colnames(n))

  # Filter all rows
  dat <- matrix(data = c(2, 4, 0, NA, 4, 5,
                         2, 4, 0, 0, 1, 0,
                         1, 2, 0, 0, 2, 3,
                         3, 9, 0, NA, 1, 9),
                nrow = 4, ncol = 6, byrow = TRUE,
                dimnames = list(c("1", "2", "3", "4"),
                                c("A_1", "A_2", "B_1", "B_2", "B_3", "C_1")))
  class(dat) <- c("matrix", "timeseries")

  n <- filter_zero_columns(dat = dat)

  expect_identical(class(dat), class(n))
  expect_equal(nrow(n), nrow(dat))
  expect_equal(ncol(n), 4)
  expect_true(all(colnames(n) %in% colnames(dat)))
  expect_equal(rownames(n), rownames(dat))
})
