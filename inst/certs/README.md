# Certified minimax mixture data (100-digit audit, mp.dps = 100)

Each CSV records the certified K-atom mixture weights and scales for
the logistic-CDF approximation

$$
\hat{\Lambda}_K(z) \;=\; \sum_{i = 1}^K w_i^\star \, \Phi(s_i^\star z),
$$

with absolute error bounded uniformly by $\Delta_K^\star$ on $\mathbb{R}$.

## Files

| File                              | $K$ | Cert form           | $\Delta_K^\star$        |
|-----------------------------------|-----|---------------------|-------------------------|
| `minimax_k01.csv` ... `k08.csv`   | 1–8 | STRICT              | M&S Table 18.4.1        |
| `minimax_k09.csv`                 | 9   | STRICT              | $3.66 \times 10^{-10}$   |
| `minimax_k10.csv`                 | 10  | STRICT              | $6.65 \times 10^{-11}$   |
| `minimax_k11.csv`                 | 11  | STRICT              | $1.258 \times 10^{-11}$  |
| `minimax_k14_BA_cert.csv`         | 14  | BA-$(K{+}1)$        | $7.233 \times 10^{-12}$  |
| `minimax_k15_BA_cert.csv`         | 15  | BA-$(K{+}1)$        | $7.211 \times 10^{-12}$  |
| `minimax_k16_BA_cert.csv`         | 16  | BA-$(K{+}1)$        | $4.066 \times 10^{-12}$  |
| `minimax_k17_BA_cert.csv`         | 17  | BA-$K$              | $2.544 \times 10^{-12}$  |
| `minimax_k22_signed_BA_cert.csv`  | 22  | signed BA-$(K{+}1)$ | $3.605 \times 10^{-15}$  |
| `deltas.csv`                      | all | n/a                 | the audited bounds themselves |

## Column format

All mixture CSVs share the schema:

```
i,weight,scale
1,<w_1>,<s_1>
...
K,<w_K>,<s_K>
```

For $K = 1, \dots, 11$ and for the relaxed $K = 14, \dots, 17$ certificates the
weights sum to 1 and are all positive, and the scales are positive and
distinct. For the signed $K = 22$ frontier the weights still sum to 1 but
include one negative atom ($\sum_i |w_i| \approx 1.024$).

`deltas.csv` carries the audited error bounds at their certified precision, one
row per $K$, with the schema

```
K,Delta,digits,cert_form,source
```

in which `digits` is the number of significant digits recorded. Each `Delta` is
the supremum of the shipped mixture's absolute error for the exact decimal
constants in the paired CSV, verified at `mp.dps = 140` (Python mpmath) by
dense critical-point search with bisection
refinement, and rounded up at the last stored digit, so every stored value is a
true upper bound for its shipped mixture. It is the
full-precision counterpart of `minimax_error(K)`, which returns the same value
rounded to a double. Every row satisfies `as.numeric(Delta) == minimax_error(K)`
exactly, so the shipped double is the correctly rounded double of the audited
constant.

The mixture constants are stored at the audited precision: 100 digits for
$K = 1, \dots, 10$ (Newton equioscillation solves at `mp.dps = 115` with
residuals $\le 5 \times 10^{-96}$, the largest being $4.3 \times 10^{-96}$
at $K = 4$ and the other nine below $10^{-106}$), 50 digits for $K = 11$ (the chain-LM
certification), and about 40 digits for the signed $K = 22$ certificate
(`mp.dps = 100` audit in Python mpmath). The $K = 14, \dots, 17$ certificates
store the scales at 15 significant digits (the double-precision search
outputs, which are the exact defining values of the certified mixtures)
and the weights at about 60 digits from the `mp.dps = 100` audit.

## The signed K = 22 frontier

`minimax_k22_signed_BA_cert.csv` is the best certified bound,
$\Delta = 3.604 \times 10^{-15}$ (audited at `mp.dps = 100`), in the signed
class. It is a self-contained 22-atom mixture, the intermediate orders
$K = 12, \dots, 21$ not being required to reach it. The bound sits within a
factor of about thirty-three of the IEEE double-precision unit roundoff $2^{-53}$,
and about eight times above the fp64 evaluation noise measured below for
$\sum_i |w_i| \approx 1.024$.

## Provenance

- K = 1..8: Monahan & Stefanski 1992, *Normal Scale Mixture
  Approximations to F*(x) and Computation of the Logistic-Normal Integral*,
  Table 18.4.1. Reproduced for this package by 100-digit Newton
  equioscillation solves (`mp.dps = 115`).
- K = 9: seeded by CMA-ES with global escape, a Remez warm-start, and a
  256-bit MPFR Newton polish. Shipped constants from the 100-digit Newton
  equioscillation re-solve (`mp.dps = 115`), which gives
  $\Delta = 3.6564505831665499 \times 10^{-10}$.
- K = 10: seeded by perturbed-reference iterated Remez and a 256-bit MPFR
  Newton polish. Shipped constants from the same 100-digit Newton
  equioscillation re-solve (`mp.dps = 115`).
- K = 11: Sturm-seeded CMA-ES + chain Levenberg-Marquardt + sequential
  200k LM at `mp.dps = 100` (Python mpmath), Newton-refined to
  alternation residuals $10^{-54}$.
- K = 14: the 14-scale subselection of the K = 15 witness, weights
  reoptimised by constrained Remez and audited at `mp.dps = 100`
  (21 dominant extrema, strictly alternating).
- K = 15: 141-seed multistart re-run of the final solver with
  constrained-Remez weight reoptimisation, audited at `mp.dps = 100`
  (20 dominant extrema).
- K = 16, 17: scale configurations from the double-precision
  linear-programming grid pipeline (the `LP_SCALE` label in the `source`
  column of `deltas.csv`), weights audited at `mp.dps = 100`.
- K = 22: alternation-aware signed-class CMA-ES, audited
  at `mp.dps = 100` (30 strictly-alternating dominant extrema).

## Verification

The shipped data is built by `data-raw/minimax_constants.R` in the source
repository. To verify the cert mixture matches the certified $\Delta_K^\star$
to within double-precision noise:

```r
library(LogisNormIntegral)
z <- seq(-30, 30, length.out = 60001)
pad <- 1e-14                      # absolute, see the note below
for (K in c(1:11, 14, 15, 16, 17, 22)) {
  err <- max(abs(plogis(z) - phi_minimax_mixture(z, K = K)))
  cat(sprintf("K = %2d: empirical max-|err| = %.3e,  certified Delta = %.3e   %s\n",
              K, err, minimax_error(K),
              if (err <= minimax_error(K) + pad) "within" else "exceeds"))
}
```

The pad is absolute rather than relative because the effect it absorbs is.
Evaluating the mixture in double precision costs approximately
$\varepsilon \sum_i |w_i|$ per term, which does not shrink with
$\Delta_K^\star$, so a relative pad would fall below the noise at the
small-$\Delta$ end. On the 60001-point grid of the snippet above that noise has
been measured at up to $4.5 \times 10^{-16}$. The package test suite asserts
the same comparison with the same pad for the positive-class certificates
($K \le 11$ and $14$ to $17$, on its coarser test grid), and asserts the
$K = 22$ evaluation against a permissive $10^{-13}$ threshold rather than against
the mathematical bound.

Only $K = 22$ shows the noise at the level of its own bound. There the certified
bound $3.604 \times 10^{-15}$ is only about eight times that noise, so the grid
maximum lands a few times $10^{-16}$ above it (about $4.0 \times 10^{-15}$) while
the exact mathematical error remains below the bound. Every $K \le 11$ row agrees
with its certified constant to the printed precision.
