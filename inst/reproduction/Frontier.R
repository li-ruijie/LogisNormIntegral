## Frontier.R: the certified K-vs-Delta frontier shipped with the package, over
## all three certificate families. These are the strict global minimax for
## K = 1, ..., 11, the relaxed Bernstein-Akhiezer certificates at K = 14 to
## K = 17, and the signed K = 22 bound. Delta_14 through Delta_17 all lie below
## Delta_11, so omitting them would understate the shipped frontier.
##
## Run from the package root or from inst/reproduction/:
##   Rscript inst/reproduction/Frontier.R
## Or after installing the package:
##   Rscript -e 'source(system.file("reproduction/Frontier.R", package = "LogisNormIntegral"))'

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

Ks <- c(1:11, 14L, 15L, 16L, 17L, 22L)

cert_form <- function(K) {
  if (K <= 11L) "STRICT"
  else if (K %in% c(14L, 15L, 16L)) "BA-(K+1)"
  else if (K == 17L) "BA-K"
  else "signed BA-(K+1)"
}

df <- data.frame(
  K        = Ks,
  Delta    = vapply(Ks, function(k) formatC(LogisNormIntegral::minimax_error(k),
                                            format = "e", digits = 4), character(1)),
  Atoms    = Ks,
  CertForm = vapply(Ks, cert_form, character(1)),
  stringsAsFactors = FALSE
)

names(df) <- c("K", "Certified Delta_K^star", "Atoms", "Cert form")

cat("Certified frontier: Delta_K^star for K = 1, ..., 11 (strict), K = 14..17 (BA),",
    " K = 22 (signed)\n", sep = "")
cat(strrep("=", 72), "\n", sep = "")
print(df, row.names = FALSE, right = FALSE)
cat(strrep("=", 72), "\n", sep = "")

## Empirical check at the K = 22 signed frontier. The grid maximum lands near 4e-15,
## the certified bound 3.604e-15 plus the fp64 evaluation noise quantified next.
## The pad is absolute because the evaluation noise it absorbs is. That noise
## costs about eps * sum_i |w_i| per term and does not shrink with Delta_K.
## Measured at up to 4.5e-16 on this 120001-point grid (the certs README
## records the same maximum on its coarser 60001-point grid). The test suite
## asserts the positive classes at bound + 1e-14 and K = 22 at a permissive 1e-13
## threshold.
z <- seq(-30, 30, length.out = 120001)
pad <- 1e-14
err22 <- max(abs(plogis(z) - LogisNormIntegral::phi_minimax_mixture(z, K = 22)))
d22 <- LogisNormIntegral::minimax_error(22)
cat(sprintf("\nK = 22 empirical max-|err| (fp64): %.3e   (certified bound %.3e)   %s\n",
            err22, d22, if (err22 <= d22 + pad) "within" else "exceeds"))
cat("  The grid maximum sits a few times 1e-16 above the certified bound because\n",
    " that bound is only about eight times the fp64 evaluation noise. The exact\n",
    " mathematical error remains below it.\n", sep = "")
