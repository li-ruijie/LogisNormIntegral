## Table1.R: Reproduce the manuscript's certified-constants table (K = 1..11).
##
## Run from the package root or from inst/reproduction/:
##   Rscript inst/reproduction/Table1.R
## Or after installing the package:
##   Rscript -e 'source(system.file("reproduction/Table1.R", package = "LogisNormIntegral"))'
##
## Output: prints the manuscript Table 1 to stdout. With OUTPUT_CSV = TRUE it
## also writes `Table1.csv` to a temporary directory (tempdir()), never to the
## working directory, in line with CRAN policy.

OUTPUT_CSV <- TRUE

## Keep the printed table in one block regardless of the console default.
options(width = 120L)

if (!requireNamespace("LogisNormIntegral", quietly = TRUE)) {
  ## Allow running directly from a clone (without installing the package).
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("Install LogisNormIntegral first or run from a clone with devtools.")
  }
  pkg_root <- normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE)
  if (!dir.exists(pkg_root) || !file.exists(file.path(pkg_root, "DESCRIPTION"))) {
    pkg_root <- normalizePath(getwd(), mustWork = TRUE)
  }
  devtools::load_all(pkg_root)
}

K_range <- 1L:11L

## The K = 1..8 label names the case treated by Monahan and Stefanski, not the
## source of the shipped values. Those values are this package's own 100-digit
## Newton equioscillation solves (`mp.dps = 115`, see `inst/certs/deltas.csv`),
## which reproduce the published Table 18.4.1 constants to about eight
## significant figures. `inst/certs/README.md` makes the same distinction and
## records the full provenance.
pipelines <- c(
  rep("M&S case, 100-digit Newton re-solve", 8L),
  "CMA-ES + Remez seed, 100-digit Newton re-solve",
  "Perturbed-reference Remez seed, 100-digit Newton re-solve",
  "Sturm-seeded CMA-ES + chain LM + seq 200k"
)

format_delta <- function(K, delta) {
  ## Sixteen significant figures (fifteen decimal digits) for K = 9, 10, 11, and
  ## five significant figures (four decimal digits) otherwise, matching the
  ## manuscript table display. Sixteen is the ceiling worth printing here,
  ## because `minimax_delta` holds doubles, and beyond the 16th figure their
  ## digits drift from the certified constants, which the manuscript quotes to
  ## 19. The audited full-precision bounds are in `inst/certs/deltas.csv`.
  if (K >= 9) {
    formatC(delta, format = "e", digits = 15)
  } else {
    formatC(delta, format = "e", digits = 4)
  }
}

df <- data.frame(
  K        = K_range,
  Delta    = vapply(K_range,
                    function(k) format_delta(k, LogisNormIntegral::minimax_error(k)),
                    character(1)),
  Extrema  = 2L * K_range,
  Pipeline = pipelines,
  stringsAsFactors = FALSE
)

names(df) <- c("K", "Certified Delta_K^star (100-digit audit)", "Extrema", "Pipeline")

cat("Table 1: Certified minimax constants Delta_K^star for K = 1, ..., 11\n")
cat(strrep("=", 90), "\n", sep = "")
print(df, row.names = FALSE, right = FALSE)
cat(strrep("=", 90), "\n", sep = "")
cat("For K = 1..8 the pipeline names the case treated by Monahan and Stefanski\n",
    "(Balakrishnan 1992, Table 18.4.1). The values above are this package's own\n",
    "100-digit Newton equioscillation solves, which reproduce their published\n",
    "constants to about eight significant figures.\n", sep = "")

cat("\nMonotonic ratio check (Delta_{K+1} / Delta_K, K = 1..10):\n")
deltas <- vapply(K_range, function(k) LogisNormIntegral::minimax_error(k), numeric(1))
ratios <- deltas[-1L] / deltas[-length(deltas)]
for (i in seq_along(ratios)) {
  cat(sprintf("  K = %2d -> %2d : ratio = %.3f\n",
              K_range[i], K_range[i + 1L], ratios[i]))
}

if (isTRUE(OUTPUT_CSV)) {
  out_path <- file.path(tempdir(), "Table1.csv")
  tryCatch({
    write.csv(df, out_path, row.names = FALSE)
    cat(sprintf("\nWrote %s\n", out_path))
  }, error = function(e) {
    cat(sprintf("\n(could not write CSV: %s)\n", conditionMessage(e)))
  })
}
