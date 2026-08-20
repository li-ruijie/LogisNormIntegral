## PirjolTiming.R: Per-point accuracy and cost of the corrected logisticNormal
##                 Poisson evaluator (Pirjol 2013) against the certified K = 11
##                 minimax mixture, on the eight-point grid of the upstream bug
##                 report (dan-pirjol/logisticNormal pull request #1).
##
## The upstream package is not shipped with this repository. Clone it first and
## pass the path to its R sources (the R/ directory of the clone, or a single
## file concatenating gPoisson, gRec, and phiLogisticNormal):
##
##   git clone https://github.com/dan-pirjol/logisticNormal
##   Rscript PirjolTiming.R logisticNormal/R [PirjolTiming.output.txt]
##
## The script applies the one-line loop-bound fix of pull request #1 in memory
## before any evaluation, so the upstream files are used unmodified on disk.
## Captured output: see PirjolTiming.output.txt next to this script.
##
## Conventions. Pirjol's phi(z, t) is the reflected integral 1 - I(z, t). The
## reference is adaptive Gauss-Kronrod quadrature of the defining integral at
## rel.tol = 1e-12. The mixture route evaluates the shipped K = 11 certificate
## (../certs/minimax_k11.csv) as 1 - sum(w * pnorm(s * z / sqrt(1 + s^2 * t))).
## Per-call costs are block means over REPS calls through an identical mapply
## harness for both routes, so the harness overhead is common to the two.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
    stop("usage: Rscript PirjolTiming.R <path to upstream R sources> [out.txt]")
}
src_path <- args[[1L]]
out_file <- if (length(args) >= 2L) args[[2L]] else "PirjolTiming.output.txt"

cmd_file <- commandArgs() |>
    Filter(function(a) startsWith(a, "--file="), x = _) |>
    (function(x) sub("^--file=", "", x[[1L]]))()
script_dir <- dirname(normalizePath(cmd_file))
cert_file <- file.path(script_dir, "..", "certs", "minimax_k11.csv")
stopifnot(file.exists(cert_file))

src_files <- if (dir.exists(src_path)) {
    file.path(src_path, c("gPoisson.R", "gRec.R", "phiLogisticNormal.R"))
} else {
    src_path
}
stopifnot(all(file.exists(src_files)))

BUG_LINE <- "for (j in 0:2*Nmax-1){"
FIX_LINE <- "for (j in 0:(2*Nmax-1)){"

patch_and_source <- function(files) {
    env <- new.env(parent = globalenv())
    hits <- Map(function(f) {
        txt <- readLines(f)
        tmp <- tempfile(fileext = ".R")
        writeLines(sub(BUG_LINE, FIX_LINE, txt, fixed = TRUE), tmp)
        sys.source(tmp, env)
        sum(grepl(BUG_LINE, txt, fixed = TRUE))
    }, files)
    n_sites <- Reduce(`+`, hits, accumulate = FALSE)
    if (n_sites != 1L) {
        stop(sprintf("expected exactly one loop-bound site upstream, found %d", n_sites))
    }
    env
}

upstream <- patch_and_source(src_files)
phi_upstream <- get("phiLogisticNormal", envir = upstream)
phi_pirjol <- function(z, t, n_terms) phi_upstream(z, t, n_terms)

## The unscaled integrand is integrated accurately on this script's grid, where |z|/sqrt(t) stays at or
## below sqrt(2), attained at (z, t) = (2, 2). A wider grid must rescale the integrand to a standard
## normal, as Coverage.R and Example1.R do, since integrate() can mis-handle a peak far from the origin
## on an infinite range.
phi_ref <- function(z, t) {
    integrate(function(x) 1 / (1 + exp(x)) * dnorm(x, z, sqrt(t)),
              -Inf, Inf, rel.tol = 1e-12, subdivisions = 4000L)[["value"]]
}

cert <- read.csv(cert_file)
mix_w <- cert[["weight"]]
mix_s <- cert[["scale"]]
phi_mix <- function(z, t) 1 - sum(mix_w * pnorm(mix_s * z / sqrt(1 + mix_s^2 * t)))

grid <- list(c(0.5, 0.5), c(1, 1), c(2, 2), c(1, 6.5),
             c(2, 10), c(1, 15), c(3, 25), c(5, 40))

err_row <- function(pt) {
    z <- pt[[1L]]
    t <- pt[[2L]]
    ref <- phi_ref(z, t)
    data.frame(z = z, t = t,
               err_n10 = abs(phi_pirjol(z, t, 10L) - ref),
               err_n30 = abs(phi_pirjol(z, t, 30L) - ref),
               err_mix = abs(phi_mix(z, t) - ref))
}
err_tab <- Map(err_row, grid) |> Reduce(rbind, x = _)

REPS <- 10000L
time_row <- function(pt) {
    z <- pt[[1L]]
    t <- pt[[2L]]
    per_call_us <- function(f) {
        t0 <- proc.time()[["elapsed"]]
        mapply(function(zz, tt) f(zz, tt), rep(z, REPS), rep(t, REPS))
        (proc.time()[["elapsed"]] - t0) / REPS * 1e6
    }
    data.frame(z = z, t = t,
               us_n10 = per_call_us(function(zz, tt) phi_pirjol(zz, tt, 10L)),
               us_n30 = per_call_us(function(zz, tt) phi_pirjol(zz, tt, 30L)),
               us_mix = per_call_us(phi_mix))
}
tim_tab <- Map(time_row, grid) |> Reduce(rbind, x = _)

stopifnot(max(err_tab[["err_n30"]]) < 1e-13)
## Guard at the certified K = 11 bound plus the absolute fp64 pad of
## tests/testthat/test-minimax.R, so a run that breached the bound could not
## pass. The printed line below displays the bound to four significant figures.
delta11 <- 1.2584376967160870e-11  # as.numeric of the K = 11 row of ../certs/deltas.csv
stopifnot(max(err_tab[["err_mix"]]) < delta11 + 1e-14)

sink(out_file)
cat("PirjolTiming.R output\n")
cat(sprintf("Date: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")))
## Record the platform rather than the machine name. The platform is what a
## reader needs to interpret the timings. The hostname identifies a private
## machine and has no reproduction value. The caller's path to the upstream
## clone has none either, so only the file names are recorded below.
cat(sprintf("Platform: %s\n", R.version$platform))
cat(sprintf("R: %s\n", R.version.string))
cat(sprintf("Upstream sources: %s (pull request #1 loop-bound fix applied in memory)\n",
            paste(basename(src_files), collapse = ", ")))
cat(sprintf("Timing: block mean over %d calls per point, identical mapply harness\n\n", REPS))
cat("Absolute error against the Gauss-Kronrod reference (rel.tol = 1e-12), phi = 1 - I:\n")
print(format(err_tab, digits = 3, scientific = TRUE), row.names = FALSE)
cat("\nPer-call cost in microseconds:\n")
print(format(tim_tab, digits = 3), row.names = FALSE)
cat(sprintf("\nMax |err|, Poisson N = 10:  %.3e\n", max(err_tab[["err_n10"]])))
cat(sprintf("Max |err|, Poisson N = 30:  %.3e\n", max(err_tab[["err_n30"]])))
cat(sprintf("Max |err|, K = 11 mixture:  %.3e (certified bound 1.258e-11)\n",
            max(err_tab[["err_mix"]])))
cat(sprintf(paste0("\nMedian per-call cost: Poisson N = 10: %.1f us, ",
                   "Poisson N = 30: %.1f us, K = 11 mixture: %.1f us\n"),
            median(tim_tab[["us_n10"]]), median(tim_tab[["us_n30"]]),
            median(tim_tab[["us_mix"]])))
sink()

cat(sprintf("done, output written to %s\n", out_file))
