## Prepare certified minimax weights and scales (K = 1..11 strict, K = 14
## to 17 relaxed BA, K = 22 signed) from the shipped inst/certs/ CSVs.
##
## Run: Rscript data-raw/minimax_constants.R
## Output: data/minimax_constants.rda

library(tools)

certs_dir <- file.path("inst", "certs")
if (!dir.exists(certs_dir)) {
  certs_dir <- file.path("..", "inst", "certs")
}
if (!dir.exists(certs_dir)) {
  stop("Cannot locate inst/certs/. Run from the package root.")
}

K_max <- 11L
minimax_constants <- vector("list", K_max)

## For K = 1..10 the shipped certificates carry the 100-digit Newton-refined
## constants (mp.dps = 115, residuals <= 5e-96, the largest being 4.3e-96 at
## K = 4 and the other nine below 1e-106). K = 11 carries the 50-digit
## constants of the chain-LM certification.
for (K in seq_len(K_max)) {
  csv_path <- file.path(certs_dir, sprintf("minimax_k%02d.csv", K))
  if (!file.exists(csv_path)) {
    stop(sprintf("Missing certified file: %s", csv_path))
  }
  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  stopifnot(all(c("i", "weight", "scale") %in% names(df)))
  df <- df[order(df$i), ]
  stopifnot(nrow(df) == K)
  stopifnot(all(df$weight > 0))
  stopifnot(abs(sum(df$weight) - 1) < 1e-12)
  stopifnot(all(df$scale > 0))
  stopifnot(length(unique(df$scale)) == K)

  minimax_constants[[K]] <- list(
    K      = K,
    weights = as.numeric(df$weight),
    scales = as.numeric(df$scale)
  )
}

## Delta_K* (minimax error), stored at full double precision (17 significant
## digits) from the 100-digit certificates, so minimax_error(K) returns the
## certified bound rather than a display rounding. Sources:
##   K=1..10: 100-digit Newton equioscillation solves (mp.dps = 115, residuals
##            <= 5e-96, the largest being 4.3e-96 at K = 4 and the other nine
##            below 1e-106). The shipped bounds are in inst/certs/deltas.csv.
##   K=11: Sturm-seeded CMA-ES + chain LM + seq200k LM at mp.dps = 100 (Python
##         mpmath), Newton-refined to alternation residuals 10^{-54}, matching
##         the certified-constants deposit in inst/certs/ to 19 digits
## History: K=1..8 were stored at 5 significant digits (nearest-rounded).
## For K=1, 3, 4 that rounding fell below the certified supremum, so the
## returned "bound" was exceeded by the shipped mixture.
minimax_delta <- c(
  9.4573091763710803e-3,
  5.1012691704537068e-4,
  4.3559344261310514e-5,
  4.7167335784457692e-6,
  5.9528866201192525e-7,
  8.3871636684423247e-8,
  1.2856920975179301e-8,
  2.1085613655930762e-9,
  3.6564505831665499e-10,
  6.6469230227209620e-11,
  1.2584376967160870e-11
)
stopifnot(length(minimax_delta) == K_max)

names(minimax_constants) <- paste0("K", seq_len(K_max))
names(minimax_delta) <- paste0("K", seq_len(K_max))

## Relaxed Bernstein-Akhiezer (BA) certs at K = 14..17 (positive class, all
## weights > 0, sum = 1). K=14, 15, and 16 are BA-(K+1) and K=17 is BA-K, which
## relax the strict 2K-equioscillation to >= K+1 dominant critical points
## (>= K at K = 17), at level >= 0.9 Delta and audited at mp.dps = 100 (Python
## mpmath). These extend the strict K<=11 hierarchy (see the cert-form table in
## inst/certs/README.md).
ba_specs <- list(
  list(K = 14L, file = "minimax_k14_BA_cert.csv", cert_form = "BA-(K+1)",
       delta = 7.232903421189757e-12),
  list(K = 15L, file = "minimax_k15_BA_cert.csv", cert_form = "BA-(K+1)",
       delta = 7.2100669764048355e-12),
  list(K = 16L, file = "minimax_k16_BA_cert.csv", cert_form = "BA-(K+1)",
       delta = 4.06514854806776278e-12),
  list(K = 17L, file = "minimax_k17_BA_cert.csv", cert_form = "BA-K",
       delta = 2.54350013158065048e-12)
)
for (spec in ba_specs) {
  csv_path <- file.path(certs_dir, spec$file)
  if (!file.exists(csv_path)) stop(sprintf("Missing BA cert file: %s", csv_path))
  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  stopifnot(all(c("i", "weight", "scale") %in% names(df)))
  df <- df[order(df$i), ]
  stopifnot(nrow(df) == spec$K)
  stopifnot(all(df$weight > 0))
  stopifnot(abs(sum(df$weight) - 1) < 1e-9)
  stopifnot(length(unique(df$scale)) == spec$K)
  key <- paste0("K", spec$K)
  minimax_constants[[key]] <- list(
    K = spec$K, weights = as.numeric(df$weight), scales = as.numeric(df$scale),
    signed = FALSE, cert_form = spec$cert_form
  )
  minimax_delta[[key]] <- spec$delta
}

## Signed K = 22 frontier (alternation-aware CMA-ES).
## Best certified bound on the integral: Delta = 3.604e-15, in the SIGNED class
## (Sum w = 1 with one negative atom, and Sum|w| ~ 1.024 is the evaluation-conditioning
## number). Cert form: signed BA-(K+1), with 30 strictly-alternating dominant
## extrema, audited at mp.dps = 100 (Python mpmath).
k22_path <- file.path(certs_dir, "minimax_k22_signed_BA_cert.csv")
if (!file.exists(k22_path)) {
  stop(sprintf("Missing signed K=22 frontier file: %s", k22_path))
}
df22 <- read.csv(k22_path, stringsAsFactors = FALSE)
stopifnot(all(c("i", "weight", "scale") %in% names(df22)))
df22 <- df22[order(df22$i), ]
stopifnot(nrow(df22) == 22)
stopifnot(abs(sum(df22$weight) - 1) < 1e-9)   # Sum w = 1 (signed, so not all weights are positive)
stopifnot(length(unique(df22$scale)) == 22)
minimax_constants[["K22"]] <- list(
  K               = 22L,
  weights         = as.numeric(df22$weight),
  scales          = as.numeric(df22$scale),
  signed          = TRUE,
  sum_abs_weights = sum(abs(as.numeric(df22$weight))),
  cert_form       = "signed BA-(K+1)"
)
minimax_delta[["K22"]] <- 3.60444176345590574e-15

## The audited bounds ship twice, as these literals and as inst/certs/deltas.csv,
## which records them at their certified precision (99 to 100 digits for
## K = 1..10, 60 for the rest). inst/certs/README.md states that as.numeric() of every CSV row
## equals minimax_error(K) exactly, so enforce it here rather than trusting two
## hand-transcribed copies to stay in step across a future re-audit.
deltas_csv <- file.path(certs_dir, "deltas.csv")
if (!file.exists(deltas_csv)) {
  stop(sprintf("Missing audited-bound file: %s", deltas_csv))
}
audited <- read.csv(deltas_csv, colClasses = "character", stringsAsFactors = FALSE)
stopifnot(all(c("K", "Delta") %in% names(audited)))
stopifnot(setequal(paste0("K", audited[["K"]]), names(minimax_delta)))
invisible(Map(
  function(K, txt) {
    key <- paste0("K", K)
    if (!identical(as.numeric(txt), minimax_delta[[key]])) {
      stop(sprintf(
        "minimax_delta[[%s]] = %.17e disagrees with %s (%.17e). Re-sync the two.",
        key, minimax_delta[[key]], deltas_csv, as.numeric(txt)))
    }
  },
  audited[["K"]], audited[["Delta"]]
))
message(sprintf("Checked %d bounds against %s", nrow(audited), deltas_csv))

save(minimax_constants, minimax_delta,
     file = file.path("data", "minimax_constants.rda"),
     version = 3, compress = "xz")

message(sprintf(
  "Wrote %d certified constants (K=1..%d strict, BA 14..17, signed 22) to data/minimax_constants.rda",
  length(minimax_constants), K_max))
