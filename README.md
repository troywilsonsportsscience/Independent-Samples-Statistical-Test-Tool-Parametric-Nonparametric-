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
required_packages <- c("tidyverse", "broom", "janitor", "ggplot2")
install.packages(setdiff(required_packages, rownames(installed.packages())))

# Optional (extra functionality: AD test, Levene, Yeo–Johnson, effect-size CI)
install.packages(c("nortest", "car", "bestNormalize", "MBESS"))

# The script includes a helper to auto-install missing packages
```
---
## 2) Data Input
- Source:
    - By default: file.choose() to import a CSV.
    - Or supply your own tibble.
- Required columns (first three columns are renamed automatically):
    - subject_id (ID)
    - group (factor with exactly two levels)
    - test_metric (numeric outcome)
- Cleaning:
    - Columns renamed with janitor::clean_names().
    - Non-finite and missing values dropped.
- Minimum N: at least 3 observations per group; Anderson–Darling requires n ≥ 8.

## 2.1) Example Data 
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
## 4) Methodology & Calculations
- Assumption Checks
    - Normality
        - By default uses Shapiro–Wilk per group (n ≥ 3).
        - Optionally Anderson–Darling (nortest, n ≥ 8 required).
        - If method_normality = "none", normality tests are skipped and data are assumed acceptable.
- Variance Homogeneity
    - Levene’s test (median-centered, via car).
    - If Levene’s p ≥ α → assume equal variances (Student’s t).
    - If Levene’s p < α or car not available → Welch’s t.
- Transformations
    - Candidate transforms: log, √, inverse, Yeo–Johnson (if bestNormalize installed).
    - Each transform is applied using a shared shift/epsilon/λ across both groups.
    - For each transform, group-wise normality p values are computed.The transform with the highest minimum p-value that passes α for both groups is chosen.
    - If no transform passes → fall back to Mann–Whitney U.
- Test Selection
    - Parametric (if assumptions satisfied or rescued):
        - Student’s t if equal variances.
        - Welch’s t otherwise.
    - Nonparametric:
        - Mann–Whitney U (Wilcoxon rank-sum).
        - Uses exact test if n1*n2 ≤ 50 and no ties; else asymptotic with continuity correction.
- Effect Sizes
    - All computed on the final analysis scale (raw or transformed).
    - Cohen’s d: mean difference ÷ pooled SD (computed even if Welch is chosen; note that it assumes homogeneity).
    - Hedges’ g: bias-corrected version of d. Optional 95% CI via MBESS::ci.smd().
    - Glass’s Δ: reported as two variants because no reference group is fixed:
        - Δ₁ = (M₁ − M₂) ÷ SD₁
        - Δ₂ = (M₁ − M₂) ÷ SD₂
        - Classic usage defines Δ relative to a “control” group; you can adapt the code to force a single reference.
    - Mann–Whitney U effect size: r = Z / √(n₁+n₂), where Z is derived from the U statistic.

---
## 5) Outputs
```r
Welch’s t-test showed that A (M = 72.10, SD = 5.02) differ from B (M = 75.30, SD = 15.20), t(43.7) = −1.98, p = .0531, 95% CI [−6.43, 0.04]. Cohen’s d (pooled) = −0.39; Hedges’ g = −0.38 (95% CI [−0.78, 0.01]); Welch’s d = −0.23; Glass’s Δ (ref A) = −0.64.
```
---
## License
MIT License. Free for academic, research, and commercial use with attribution.

---
## Contact
Created and maintained by Troy Wilson with the assistance of OpenAI's ChatGPT to help structure code, explain statistical logic, and improve documentation clarity. For issues or feature requests, please use the GitHub Issues tab. 
