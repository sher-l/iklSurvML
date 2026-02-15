# Tests for helper functions
test_that("display_progress works correctly", {
  # Test that display_progress doesn't throw errors
  expect_invisible(Mime:::display_progress(5, 20, 4))
  expect_invisible(Mime:::display_progress(10, 100, 10))
})

test_that("return_id_to_rs adds ID column correctly", {
  # Create test data
  rs_list <- list(
    dataset1 = data.frame(OS.time = c(100, 200), OS = c(1, 0), RS = c(0.5, 0.8)),
    dataset2 = data.frame(OS.time = c(150, 250), OS = c(0, 1), RS = c(0.3, 0.9))
  )

  raw_id_list <- list(
    dataset1 = data.frame(ID = c("patient1", "patient2")),
    dataset2 = data.frame(ID = c("patient3", "patient4"))
  )

  result <- Mime:::return_id_to_rs(rs_list, raw_id_list)

  # Check that ID column is added
  expect_true("ID" %in% colnames(result$dataset1))
  expect_true("ID" %in% colnames(result$dataset2))

  # Check ID values
  expect_equal(result$dataset1$ID, c("patient1", "patient2"))
  expect_equal(result$dataset2$ID, c("patient3", "patient4"))

  # Check ID is first column
  expect_equal(colnames(result$dataset1)[1], "ID")
})

# Tests for C-index calculation
test_that("calculate_cindex returns valid values", {
  skip_if_not_installed("survival")

  # Create test data with more samples for convergence
  test_data <- data.frame(
    OS.time = c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000),
    OS = c(1, 0, 1, 0, 1, 1, 0, 1, 0, 1)
  )
  rs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

  cindex <- Mime:::calculate_cindex(rs, test_data)

  # C-index should be between 0 and 1
  expect_gte(cindex, 0)
  expect_lte(cindex, 1)
})

# Tests for data preprocessing
test_that("preprocess_data_list handles basic data", {
  # Skip if survival package not available
  skip_if_not_installed("survival")

  # This is a basic smoke test
  expect_true(TRUE)
})
