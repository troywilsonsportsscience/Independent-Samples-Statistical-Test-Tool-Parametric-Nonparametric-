# Independent-Samples Statistical Test Tool (R)
An all-in-one R script for automated independent-samples hypothesis testing with built-in assumption checks, safe transforms, effect sizes, and APA-style reporting.

## Features
- Assumption checks
    - Per-group normality (Shapiro–Wilk or Anderson–Darling)
    - Variance homogeneity via Levene’s test (median-centered)
- Automatic decision logic
    - Student’s t (equal variances)
    - Welch’s t (unequal variances or if Levene unavailable)
    - Mann–Whitney U (Wilcoxon rank-sum) if normality assumptions cannot be satisfied
- Safe transformations
    - Log, square root, inverse, Yeo–Johnson
    - Applied consistently across both groups using a shared shift/λ
    - Selects the transform that maximizes the minimum group normality p ≥ α
- Effect sizes
    - Parametric: Cohen’s d (pooled SD), Hedges’ g (bias-corrected d), Glass’s Δ (both variants), optional 95% CI for g (via MBESS)
    - Nonparametric: rank-biserial r (derived from Mann–Whitney Z)
- Outputs
    - APA-style text (with assumptions & interpretation)
    - Group descriptives (M, SD, Median, n)
    - Tidy tibble with test results & effect sizes
- Visualization
    - Histograms + densities by group
    - Q–Q plots by group
    - Transform comparison panels (optional)

---
## 1) Install
```r
# R >= 4.1 recommended
required_packages <- c("tidyverse", "car", "broom", "janitor", "rstatix")
install.packages(setdiff(required_packages, rownames(installed.packages())))

# Optional (for extra tests/intervals)
install.packages(c("nortest", "MBESS", "effsize"))

# The script includes a helper to auto-install missing packages
```
---

## 2) Example Data 
### A) Classic t-test (normal & equal variances)
```r
set.seed(101)
test_data <- tibble(
  group = factor(rep(c("A","B"), each = 30)),
  score = c(rnorm(30, mean = 72, sd = 8),
            rnorm(30, mean = 75, sd = 8))
)
```
### B) Welch's t-test (normal & unequal variance)
```r
set.seed(202)
test_data <- tibble(
  group = factor(rep(c("A","B"), each = 30)),
  score = c(rnorm(30, mean = 72, sd = 5),   # smaller variance
            rnorm(30, mean = 75, sd = 15))  # larger variance
)
```

### C) Mann–Whitney U (clearly non-normal; usually unfixable)
```r
set.seed(303)
test_data <- tibble(
  group = factor(rep(c("A","B"), each = 25)),
  score = c(rexp(25, rate = 0.1),
            rexp(25, rate = 0.1) + 10)
)
```

### D) Transform “rescue” (log transform usually fixes normality)
```r
set.seed(789)
test_data <- tibble(
  group = factor(rep(c("A","B"), each = 35)),
  score = c(
    rgamma(35, shape = 2, scale = 2),  # positively skewed
    rgamma(35, shape = 2, scale = 5)   # more skewed, higher variance
  )
)
```
---

## 3) Settings (Section 0)
```r
alpha_level     <- 0.05                   # significance
tail_type       <- "two.sided"            # "two.sided", "greater", "less"
metric_label    <- "exam score"           # label in outputs
group_variable  <- "group"                # factor with 2 levels
outcome_variable<- "score"                # numeric outcome

try_transform   <- TRUE                   # attempt log/sqrt/inv if non-normal
allowed_transforms <- c("log","sqrt","inv")

show_intermediate_plots <- TRUE           # raw/final histograms + boxplots
verbose         <- TRUE                   # print messages/tables
method_normality <- "shapiro"             # "shapiro", "ad", or "none"

report_cis_effect_sizes <- TRUE           # needs MBESS/effsize if available
use_glass_delta <- TRUE                   # recommend when variances unequal
glass_reference_group <- "A"              # SD reference for Glass’s Δ
```
---

## 5) Outputs
```r
Welch’s t-test showed that A (M = 72.10, SD = 5.02) differ from B (M = 75.30, SD = 15.20), t(43.7) = −1.98, p = .0531, 95% CI [−6.43, 0.04]. Cohen’s d (pooled) = −0.39; Hedges’ g = −0.38 (95% CI [−0.78, 0.01]); Welch’s d = −0.23; Glass’s Δ (ref A) = −0.64.
```
---
## 6) Methods & Rationale
- Normality: per-group Shapiro–Wilk by default; method_normality = "ad" (Anderson–Darling) optional; or "none" to skip tests and rely on robustness/visuals.
- Variance equality: Levene (median-centered; Brown–Forsythe) determines classic vs. Welch.
- Transforms: log, sqrt, inv applied safely (shifts & epsilons recorded) to salvage parametric assumptions where appropriate.
- Effect sizes:
    - Unequal variances → prefer Welch’s d or Glass’s Δ over pooled‐SD d.
    - Nonparametric → rank-biserial r and Cliff’s δ.
- Direction/sign: Effect signs are locked to the first factor level (reference group) for consistency.
---

## License
MIT License. Free for academic, research, and commercial use with attribution.

---
## Contact
Created and maintained by Troy Wilson with the assistance of OpenAI's ChatGPT to help structure code, explain statistical logic, and improve documentation clarity. For issues or feature requests, please use the GitHub Issues tab. 
