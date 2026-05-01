test_that("packaged ImmRes portal script does not expose static login credentials", {
  master_path <- system.file("extdata", "ImmRes_master.R",
                             package = "iklSurvML", mustWork = TRUE)
  master_text <- paste(readLines(master_path, warn = FALSE), collapse = "\n")

  expect_false(grepl("icr\\.review", master_text, ignore.case = TRUE))
  expect_false(grepl("password\\s+[A-Za-z0-9._@-]+", master_text, ignore.case = TRUE))
})

test_that("packaged ImmRes master script does not source helpers from the working directory", {
  master_path <- system.file("extdata", "ImmRes_master.R",
                             package = "iklSurvML", mustWork = TRUE)
  master_text <- paste(readLines(master_path, warn = FALSE), collapse = "\n")

  expect_false(grepl('source\\("ImmRes_source\\.R"\\)', master_text))
  expect_true(grepl("system.file", master_text, fixed = TRUE))
  expect_true(grepl("mustWork = TRUE", master_text, fixed = TRUE))
})
