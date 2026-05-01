# Load IMPRES helper scripts into the caller environment without attaching
# package namespaces or example data frames to the global search path.
.impres_env <- environment()
.impres_packages <- c(
  "matrixStats", "plyr", "ppcor", "survival", "ROCR", "Hmisc",
  "rms", "mixtools", "lme4", "lmerTest", "plotrix"
)
.impres_missing <- .impres_packages[
  !vapply(.impres_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(.impres_missing) > 0) {
  stop(
    paste0("Missing IMPRES helper packages: ",
           paste(.impres_missing, collapse = ", ")),
    call. = FALSE
  )
}

for (.impres_pkg in .impres_packages) {
  for (.impres_nm in getNamespaceExports(.impres_pkg)) {
    assign(.impres_nm, getExportedValue(.impres_pkg, .impres_nm), envir = .impres_env)
  }
}

.impres_files <- c(
  "ImmRes_output.R",
  "ImmRes_generic.R",
  "ImmRes_OE.R",
  "ImmRes_master.R",
  "ImmRes1_denovoCellTypeSig.R",
  "ImmRes2_immuneResistanceProgram.R",
  "ImmRes3_longitudinal.R",
  "ImmRes4_predictICBresponses.R",
  "ImmRes5_valCohort2.R",
  "ImmRes6_pancanDrug.R"
)

for (.impres_file in .impres_files) {
  sys.source(
    system.file("extdata", .impres_file, package = "iklSurvML", mustWork = TRUE),
    envir = .impres_env
  )
}

rm(
  .impres_env, .impres_packages, .impres_missing, .impres_pkg,
  .impres_nm, .impres_files, .impres_file
)
