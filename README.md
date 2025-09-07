# Independent-Samples-Statistical-Test-Tool-Parametric-Nonparametric-
This R script provides a fully automated framework for conducting **independent-samples hypothesis tests** to compare two groups on a continuous numeric outcome. It is designed for flexible reuse across datasets and projects in **data science, human performance, and research**.

---

## Features

- **Automatic assumption testing**:
  - Shapiro-Wilk test for normality (per group)
  - Levene’s test for homogeneity of variance
- **Dynamic test selection**:
  - Independent-samples t-test (equal variances)
  - Welch’s t-test (unequal variances)
  - Mann-Whitney U test (nonparametric fallback)
- **Optional transformations** (`log`, `sqrt`, `inv`) to attempt restoring normality
- **Visualization**:
  - Histograms with density overlays per group
  - Boxplots for group comparison
  - Annotated p-values for quick diagnostics
- **APA-style reporting**:
  - Clear test summaries with test statistics, p-values, and confidence intervals
  - Cohen’s *d* effect size (parametric only)
- **Tidy output tables** for downstream use

---

## Key Parameters (User-Editable)

- `alpha_level`: significance threshold (default = 0.05)  
- `tail_type`: `"two.sided"`, `"greater"`, or `"less"`  
- `metric_label`: friendly name for the variable (used in reporting)  
- `group_variable`: name of the grouping variable  
- `outcome_variable`: name of the numeric outcome variable  
- `try_transform`: attempt transformations if normality fails (`TRUE`/`FALSE`)  
- `allowed_transforms`: vector of transforms allowed (`c("log", "sqrt", "inv")`)  
- `show_intermediate_plots`: show assumption-check plots (`TRUE`/`FALSE`)  

---

## Dependencies

This tool requires the following R packages:

- [tidyverse](https://www.tidyverse.org/)  
- [car](https://cran.r-project.org/package=car)  
- [broom](https://cran.r-project.org/package=broom)  
- [rstatix](https://cran.r-project.org/package=rstatix)  
- [janitor](https://cran.r-project.org/package=janitor)  

The script includes a helper to **auto-install missing packages**.

---

## Usage

### 1. Load Your Data
```r
df <- read_csv("your_file.csv") %>%
  janitor::clean_names()

test_data <- df %>%
  select(group, score) %>%
  drop_na()
```
Or assign manully for testing/learning
```r
test_data <- tibble(
  group = c("A","A","B","B"),
  score = c(72.3, 75.1, 80.4, 83.2)
)
```

## Exampel Outputs
### Parametric (Independent-samples and Welch’s t-test)
An independent-samples t-test showed that Group A (M = 72.60, SD = 8.10) 
and Group B (M = 78.30, SD = 9.50) differed significantly in exam score, 
t(58) = -2.35, p = 0.022, 95% CI [-10.60, -0.85]. Cohen’s d = -0.61

### Nonparametric - Mann-Whitney U
A Mann-Whitney U test indicated a significant difference in exam score 
between groups, W = 234.0, p = 0.014. Group medians: A = 70.0, B = 78.0.

---

## Example Simulated Datasets
The script includes built-in examples to test different scenarios:
- Normal, equal variance
- Skewed distributions (non-normal)
- Unequal variances
Just uncomment the dataset you want to run.

---

## Notes
- Cohen’s d is computed with a pooled SD. For Welch’s t-test (unequal variances), interpretation may be less appropriate.
- Nonparametric tests (Mann-Whitney U) do not compute Cohen’s d; effect size must be estimated separately.

---

## License
MIT License. Free for academic, research, and commercial use with attribution.

---
## Contact
Created and maintained by Troy Wilson. For issues or feature requests, please use the GitHub Issues tab. This project was developed with the assistance of OpenAI's ChatGPT to help structure code, explain statistical logic, and improve documentation clarity.
