find_package_source_root <- function() {
  candidates <- c(
    getwd(),
    file.path(getwd(), ".."),
    file.path(getwd(), "..", ".."),
    file.path(getwd(), "00_pkg_src", "iklSurvML"),
    file.path(getwd(), "..", "00_pkg_src", "iklSurvML"),
    file.path(getwd(), "..", "..", "00_pkg_src", "iklSurvML"),
    testthat::test_path("..", "..")
  )

  candidates <- normalizePath(candidates, mustWork = FALSE)
  candidates[file.exists(file.path(candidates, "DESCRIPTION")) &
               dir.exists(file.path(candidates, "inst", "extdata"))][1]
}

test_that("dead duplicate IMPRES source tree is removed", {
  pkg_root <- find_package_source_root()
  skip_if(is.na(pkg_root), "package source root is unavailable")

  expect_false(dir.exists(file.path(pkg_root, "R", "IMPRES")))
})

test_that("runtime IMPRES extdata helpers do not hard-code the legacy RNG seed", {
  pkg_root <- find_package_source_root()
  skip_if(is.na(pkg_root), "package source root is unavailable")

  impres_files <- list.files(
    file.path(pkg_root, "inst", "extdata"),
    pattern = "^ImmRes.*\\.R$",
    full.names = TRUE
  )
  skip_if_not(length(impres_files) > 0, "IMPRES extdata source files are unavailable")

  impres_source <- unlist(lapply(impres_files, readLines, warn = FALSE), use.names = FALSE)

  expect_false(any(grepl("set\\.seed\\(1234\\)", impres_source)))
})

test_that("IMPRES source loader avoids search-path side effects", {
  pkg_root <- find_package_source_root()
  skip_if(is.na(pkg_root), "package source root is unavailable")

  source_text <- readLines(file.path(pkg_root, "inst", "extdata", "ImmRes_source.R"), warn = FALSE)
  expect_false(any(grepl("^\\s*library\\(", source_text)))
  expect_false(any(grepl("\\battach\\(", source_text)))

  before <- search()
  impres_env <- iklSurvML:::load_impres_extdata_helpers()
  after <- search()

  expect_equal(after, before)
  expect_true(exists("get.OE.bulk", envir = impres_env, inherits = FALSE))
})

test_that("deterministic IMPRES extdata helpers do not reset the caller RNG stream", {
  pkg_root <- find_package_source_root()
  skip_if(is.na(pkg_root), "package source root is unavailable")

  impres_dir <- file.path(pkg_root, "inst", "extdata")
  impres_env <- new.env(parent = globalenv())
  sys.source(file.path(impres_dir, "ImmRes_generic.R"), envir = impres_env)
  sys.source(file.path(impres_dir, "ImmRes_OE.R"), envir = impres_env)

  r <- list(
    tpm = matrix(seq_len(12), nrow = 4, dimnames = list(paste0("G", 1:4), paste0("S", 1:3))),
    genes = paste0("G", 1:4),
    rand.scores = matrix(0, nrow = 3, ncol = 1, dimnames = list(paste0("S", 1:3), "Sig"))
  )
  gene_sign <- list(Sig = c("G1", "G2"))

  had_rng <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  saved_rng <- if (had_rng) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    if (had_rng) {
      assign(".Random.seed", saved_rng, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  if (!had_rng) {
    invisible(stats::runif(1))
  }
  baseline_rng <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  expected_next <- stats::runif(3)

  assign(".Random.seed", baseline_rng, envir = .GlobalEnv)
  invisible(impres_env$get.OE.bulk.specific.tpm(r, gene_sign))
  observed_next <- stats::runif(3)

  expect_equal(observed_next, expected_next)
})
