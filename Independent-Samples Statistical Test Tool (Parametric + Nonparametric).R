# ===============================================================
# Independent-Samples Statistical Test Tool (Parametric + Nonparametric)
# ===============================================================
# 
# Performs an automated independent-samples hypothesis test for two groups
# on a continuous outcome. The script:
#   - Checks normality (Shapiro-Wilk) and variance equality (Levene’s test)
#   - Selects the appropriate test (independent t-test, Welch’s t-test, or Mann–Whitney U)
#   - Outputs APA-style results and summary statistics
#   - Optionally applies data transformations and generates plots
#
# Full documentation, usage examples, and reporting guidelines:
# See the project README on GitHub.

# -------------------------------
# Helper: Install and Load Packages
# -------------------------------
# These packages are required for data manipulation, visualization, assumption checks, and APA-style reporting.

required_packages <- c("tidyverse", "car", "broom", "rstatix", "janitor")

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

invisible(lapply(required_packages, install_if_missing))

# ---- 0. Set-up ----
# Define test parameters (edit as needed for your analysis)
alpha_level <- 0.05 # Significance threshold for hypothesis tests

tail_type <- "two.sided" # Type of test: "two.sided", "greater", or "less"

metric_label <- "exam score" # Friendly label for the variable being tested (for APA reporting)

group_variable <- "group" # Name of the grouping variable

outcome_variable <- "score" # Name of the outcome/numeric variable

# Specify which transformations to allow
allowed_transforms <- c("log", "sqrt")  # Options: "log", "sqrt", "inv"


#--- Optional Settings

# Try transformations if normality assumption fails (TRUE = attempt to attempt to fix non-normality)
try_transform <- FALSE # FALSE defaults to Mann-Whitney U if normality fails.

# Show diagnostic plots during assumption testing (TRUE = show plots in 2A/2A.1)
show_intermediate_plots <- TRUE # FALSE suppresses plots during normality assumption checks.


# ---- 1. Real Data Input ----
# Load or assign your real dataset here
# Example structure: your_data <- read_csv("your_file.csv") OR assign directly
# Your data frame must include a grouping variable and a numeric outcome variable.


# Uncomment and edit below:
# path <- file.choose() # Select your CSV file
# df <- read_csv(path, show_col_types = FALSE) %>% clean_names() # Clean column names
#
# test_data <- df %>%
# select(all_of(c(group_variable, outcome_variable))) %>%
# drop_na() # Ensure no missing values in selected columns


# If you're pasting in a data frame manually, it should look like:
# test_data <- tibble(
# group = c("A", "A", "B", "B", ...),
# score = c(72.3, 75.1, 80.4, 83.2, ...)
# )


# ---- 1B. Simulated Example Data (Optional) ----
# Below are several test datasets to evaluate the tool under various data conditions.


# Simulated Normally Distributed, Equal Variance
# set.seed(101)
# test_data <- tibble(
#   group = factor(rep(c("A", "B"), each = 30)),
#   score = c(rnorm(30, mean = 70, sd = 8),
#             rnorm(30, mean = 78, sd = 8))
# )


# # Simulated Non-Normal (Skewed) Distribution
# set.seed(202)
# test_data <- tibble(
# group = factor(rep(c("A", "B"), each = 30)),
# score = c(rexp(30, rate = 1/70),
# rexp(30, rate = 1/80))
# )

# Simulated Non-Normal (Skewed) Distribution with Slight Shift - transform attempt TRUE
# set.seed(606)
# test_data <- tibble(
#   group = factor(rep(c("A", "B"), each = 30)),
#   score = c(rlnorm(30, meanlog = 4.2, sdlog = 0.25),   # Right-skewed, mild
#             rlnorm(30, meanlog = 4.3, sdlog = 0.25))   # Slightly shifted
# )

# Simulated Unequal Variance
set.seed(303)
test_data <- tibble(
  group = factor(rep(c("A", "B"), each = 30)),
  score = c(rnorm(30, mean = 72, sd = 5),
            rnorm(30, mean = 78, sd = 15))
)

# ---- 2. Assumption Checks ----
# 2A. Check Normality Within Each Group
cat("\n---- Assumption Check: Normality ----\n")
norm_results <- test_data %>%
  group_by(!!sym(group_variable)) %>%
  summarise(
    p_value = shapiro.test(!!sym(outcome_variable))$p.value
  ) %>%
  rename(group = 1)
print(norm_results)

# Output interpretation per group
for (i in 1:nrow(norm_results)) {
  g <- norm_results$group[i]
  p <- norm_results$p_value[i]
  result <- ifelse(p > alpha_level, "✅ Normal", "❌ Non-normal")
  cat(sprintf("Group %s: SW p = %.4f → %s\n", g, p, result))
}

normality_ok <- all(norm_results$p_value > alpha_level)
normality_ok_original <- normality_ok  # Save pre-transformation result - needed for APA reporting
cat("\nInterpretation: p >", alpha_level, "→ assumption of normality met; p <=", alpha_level, "→ possible non-normality.\n")

# Optional: show histograms for raw outcome variable by group
if (show_intermediate_plots) {
  cat("\n↪ Showing histograms of raw outcome variable by group...\n")
  
  annotated_data <- test_data %>%
    left_join(norm_results, by = group_variable) %>%
    mutate(facet_label = paste0(group, " (SW p = ", sprintf("%.3f", p_value), ")"))
  
  p_raw <- ggplot(annotated_data, aes_string(x = outcome_variable, fill = group_variable)) +
    geom_histogram(aes(y = ..density..), bins = 15, color = "black", alpha = 0.6, position = "identity") +
    geom_density(alpha = 0.3) +
    facet_wrap(~facet_label) +
    labs(title = "Histogram of Raw Outcome by Group with SW p-values", x = metric_label, y = "Density") +
    theme_minimal()
  print(p_raw)
}

# 2A.1. Optional Transformation Attempt (if enabled and normality is violated)
if (!normality_ok && try_transform) {
  cat("\n↪ Normality violation detected. Attempting transformations...\n")
  
  # Define and apply user-specified transforms dynamically
  available_transforms <- list(
    log = function(x) log(x),
    sqrt = function(x) sqrt(x),
    inv = function(x) 1 / x
  )
  
  # Filter to only allowed transforms
  selected_transforms <- available_transforms[names(available_transforms) %in% allowed_transforms]
  
  if (length(selected_transforms) == 0) {
    cat("\n⚠️ No valid transformations in `allowed_transforms`. Skipping transformation attempt.\n")
  } else {
    # Apply transforms and reshape
    transformed <- map_dfc(selected_transforms, ~ .x(test_data[[outcome_variable]])) %>%
      set_names(paste0(names(selected_transforms), "_score")) %>%
      bind_cols(test_data, .) %>%
      pivot_longer(
        cols = ends_with("_score"),
        names_to = "transform_type",
        values_to = "transformed_value"
      ) %>%
      mutate(transform_type = str_remove(transform_type, "_score"))
    
    # Run Shapiro-Wilk tests on each transformation
    transform_results <- transformed %>%
      group_by(transform_type, !!sym(group_variable)) %>%
      summarise(p_value = shapiro.test(transformed_value)$p.value, .groups = "drop")
    
    print(transform_results)
    
    # Optional diagnostic plots
    if (show_intermediate_plots) {
      cat("\n↪ Showing histograms for transformed values...\n")
      facet_data <- transformed %>%
        left_join(transform_results, by = c("transform_type", group_variable)) %>%
        mutate(facet_label = paste0(transform_type, " - ", .data[[group_variable]],
                                    " (SW p = ", sprintf("%.3f", p_value), ")"))
      
      p_trans <- ggplot(facet_data, aes(x = transformed_value)) +
        geom_histogram(aes(y = ..density..), bins = 15, fill = "skyblue", color = "black") +
        geom_density(col = "red", linetype = "dashed", size = 1) +
        facet_wrap(~facet_label, scales = "free", ncol = 2) +
        labs(title = "Distribution of Transformed Scores by Group and Method",
             x = "Transformed Value", y = "Density") +
        theme_minimal()
      print(p_trans)
    }
    
    # Select best transformation if it fixes all groups
    viable_transforms <- transform_results %>%
      group_by(transform_type) %>%
      summarise(min_p = min(p_value), .groups = "drop") %>%
      filter(min_p > alpha_level)
    
    if (nrow(viable_transforms) > 0) {
      best_transform <- viable_transforms %>% arrange(desc(min_p)) %>% slice(1) %>% pull(transform_type)
      test_data <- test_data %>%
        mutate(score = transformed %>%
                 filter(transform_type == best_transform) %>%
                 arrange(!!sym(group_variable)) %>%
                 pull(transformed_value))
      cat("\n✔️ Transformation successful using:", best_transform, "→ proceeding with parametric testing.\n")
      normality_ok <- TRUE
    } else {
      cat("\n❌ No transformation restored normality. Proceeding with nonparametric test.\n")
    }
  }
}

# ---- 2B. Check Homogeneity of Variance ----
if (normality_ok) {
  cat("\n---- Assumption Check: Equal Variance ----\n")
  formula_obj <- reformulate(group_variable, outcome_variable)
  levene_res <- leveneTest(formula_obj, data = test_data)
  print(levene_res)
  
  p_var <- levene_res$`Pr(>F)`[1]
  if (p_var > alpha_level) {
    cat(sprintf("\n✅ Equal variances confirmed (Levene's test p = %.4f > %.2f).\n", p_var, alpha_level))
  } else {
    cat(sprintf("\n❌ Variance assumption violated (Levene's test p = %.4f <= %.2f).\n", p_var, alpha_level))
    cat("↪ Welch's t-test will be used instead of classic t-test.\n")
  }
  
  cat("\nInterpretation: p >", alpha_level, "→ variances are equal; p <=", alpha_level, "→ variances differ.\n")
  
} else {
  cat("\n⏭ Skipping Levene's test due to violation of normality in one or both groups.\n")
}
variance_ok <- if (exists("levene_res")) levene_res$`Pr(>F)`[1] > alpha_level else FALSE


# ---- 2C. Visual Inspection ----
cat("\n---- Visualization: Distribution by Group ----\n")
cat("Use the plots below to visually inspect:")
cat("\n - Whether each group appears symmetric (supports normality)")
cat("\n - Whether group distributions have similar spread (supports homogeneity of variance)")
cat("\n - Whether there are extreme outliers or limited overlap between groups\n")

# Recalculate SW p-values for current (possibly transformed) score
sw_current <- test_data %>%
  group_by(!!sym(group_variable)) %>%
  summarise(
    p_value = shapiro.test(!!sym(outcome_variable))$p.value,
    .groups = "drop"
  )

# Annotate data with updated SW values
final_annotated <- test_data %>%
  left_join(sw_current, by = group_variable) %>%
  mutate(facet_label = paste0(.data[[group_variable]], " (SW p = ", sprintf("%.3f", p_value), ")"))

# Detect label of final test variable type
data_status <- if (exists("best_transform")) paste("Transformed (", best_transform, ")", sep = "") else "Raw"

# Histogram of final variable (raw or transformed)
p_final_hist <- ggplot(final_annotated, aes_string(x = outcome_variable, fill = group_variable)) +
  geom_histogram(aes(y = ..density..), bins = 15, alpha = 0.6, color = "black", position = "identity") +
  geom_density(alpha = 0.3) +
  facet_wrap(~facet_label) +
  labs(title = "Histogram of Final Test Variable by Group with SW p-values",
       subtitle = paste("Data used for test:", data_status),
       x = metric_label, y = "Density") +
  theme_minimal()
print(p_final_hist)

# Boxplot of final test variable
p_box <- ggplot(test_data, aes_string(x = group_variable, y = outcome_variable, fill = group_variable)) +
  geom_boxplot(alpha = 0.7, outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  labs(title = "Boxplot of Group Distributions (Final Test Variable)",
       subtitle = paste("Data used for test:", data_status),
       y = metric_label, x = "Group") +
  theme_minimal()
print(p_box)


# ---- Section 3: Selecting and Running the Appropriate Test ----
cat("\n---- Section 3: Selecting and Running the Appropriate Test ----\n")

# normality_ok is already set in 2A/2A.1
# Only calculate variance_ok if Levene's result exists
variance_ok <- if (exists("levene_res")) levene_res$`Pr(>F)`[1] > alpha_level else FALSE

# Select and run appropriate test based on assumptions
if (normality_ok && variance_ok) {
  test_method <- "classic_t"
  test_res <- t.test(
    formula = reformulate(group_variable, outcome_variable),
    data = test_data,
    var.equal = TRUE,
    alternative = tail_type
  )
  cat("✅ Normality and equal variances met → Running classic independent-samples t-test.\n")
  
} else if (normality_ok && !variance_ok) {
  test_method <- "welch_t"
  test_res <- t.test(
    formula = reformulate(group_variable, outcome_variable),
    data = test_data,
    var.equal = FALSE,
    alternative = tail_type
  )
  cat("⚠️ Normality met but unequal variances detected → Running Welch’s t-test.\n")
  
} else {
  test_method <- "nonparametric"
  test_res <- wilcox.test(
    formula = reformulate(group_variable, outcome_variable),
    data = test_data,
    alternative = tail_type
  )
  cat("❗ Normality assumption not met → Running Mann-Whitney U test (Wilcoxon rank-sum).\n")
}

cat("Test method selected:", test_method, "\n")

# ---- Section 4: Summary Statistics and APA Reporting ----
cat("\n---- Section 4: Summary Statistics and APA-style Output ----\n")

# Compute group-level descriptive statistics
summary_stats <- test_data %>%
  group_by(!!sym(group_variable)) %>%
  summarise(
    Mean = mean(!!sym(outcome_variable), na.rm = TRUE),
    SD = sd(!!sym(outcome_variable), na.rm = TRUE),
    Median = median(!!sym(outcome_variable), na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )
print(summary_stats)

# Output APA-style interpretation
if (test_method != "nonparametric") {
  tidy_res <- broom::tidy(test_res)
  
  df <- round(tidy_res$parameter, 1)
  tval <- round(tidy_res$statistic, 2)
  pval <- round(tidy_res$p.value, 4)
  ci_low <- round(tidy_res$conf.low, 2)
  ci_high <- round(tidy_res$conf.high, 2)
  
  # Compute Cohen's d using pooled SD
  mean_diff <- summary_stats$Mean[1] - summary_stats$Mean[2]
  pooled_sd <- sqrt(((summary_stats$n[1] - 1) * summary_stats$SD[1]^2 +
                       (summary_stats$n[2] - 1) * summary_stats$SD[2]^2) /
                      (summary_stats$n[1] + summary_stats$n[2] - 2))
  cohens_d <- round(mean_diff / pooled_sd, 2)
  
  # Determine test direction label for interpretation
  direction_text <- switch(
    tail_type,
    "two.sided" = "differ",
    "greater" = "is greater than",
    "less" = "is less than"
  )
  
  # Inform user if data was transformed before analysis
  if (exists("normality_ok_original") && !normality_ok_original && try_transform && exists("best_transform")) {
    cat(sprintf(
      "\nNote: The variable '%s' was transformed using '%s' to meet normality assumptions.\n",
      outcome_variable, best_transform
    ))
  }
  
  # Label test type for APA sentence
  test_type_label <- if (exists("best_transform") && !normality_ok_original && try_transform) {
    paste("A", best_transform, "transformed", ifelse(test_method == "welch_t", "Welch’s t-test", "independent-samples t-test"))
  } else {
    ifelse(test_method == "welch_t", "Welch’s t-test", "An independent-samples t-test")
  }
  
  cat(sprintf(
    "\n%s showed that Group %s (M = %.2f, SD = %.2f) and Group %s (M = %.2f, SD = %.2f) %s significantly in %s, t(%.1f) = %.2f, p = %.4f, 95%% CI [%.2f, %.2f]. Cohen’s d = %.2f\n",
    test_type_label,
    summary_stats[[group_variable]][1], summary_stats$Mean[1], summary_stats$SD[1],
    summary_stats[[group_variable]][2], summary_stats$Mean[2], summary_stats$SD[2],
    direction_text,
    metric_label, df, tval, pval, ci_low, ci_high, cohens_d
  ))
  
  if (test_method == "welch_t") {
    cat("Note: Cohen’s d is based on pooled SD and may be less appropriate with unequal variances (Welch’s test).\n")
  }
  
} else {
  stat <- round(test_res$statistic, 2)
  pval <- round(test_res$p.value, 4)
  
  cat(sprintf(
    "\nA Mann-Whitney U test indicated a %s difference in %s between groups, W = %.2f, p = %.4f. Group medians: %s = %.2f, %s = %.2f.\n",
    ifelse(pval < alpha_level, "significant", "non-significant"),
    metric_label,
    stat, pval,
    summary_stats[[group_variable]][1], summary_stats$Median[1],
    summary_stats[[group_variable]][2], summary_stats$Median[2]
  ))
  cat("ℹ️ Effect size (e.g., Cohen’s d) not computed for nonparametric tests.\n")
}


# ---- 5. Final Output Table (Optional Tidy Summary) ----
cat("\n---- Section 5: Tidy Output Table ----\n")

# Label transformation status
transform_status <- if (exists("best_transform") && !normality_ok_original && try_transform) {
  best_transform
} else {
  "None"
}

if (test_method != "nonparametric") {
  tidy_out <- tibble::tibble(
    Method = test_res$method,
    `Test Type` = ifelse(test_method == "classic_t", "Independent t-test", "Welch's t-test"),
    `Test Direction` = tail_type,
    `Group 1` = summary_stats[[group_variable]][1],
    `Group 2` = summary_stats[[group_variable]][2],
    `Mean 1` = round(summary_stats$Mean[1], 3),
    `SD 1` = round(summary_stats$SD[1], 3),
    `Mean 2` = round(summary_stats$Mean[2], 3),
    `SD 2` = round(summary_stats$SD[2], 3),
    `t(df)` = paste0("t(", round(test_res$parameter, 1), ") = ", round(test_res$statistic, 2)),
    `p-value` = round(test_res$p.value, 4),
    `95% CI` = paste0("[", round(test_res$conf.int[1], 2), ", ", round(test_res$conf.int[2], 2), "]"),
    `Transformed?` = transform_status
  )
  print(tidy_out)
  
} else {
  tidy_out <- tibble::tibble(
    Method = "Mann-Whitney U (Wilcoxon Rank-Sum Test)",
    `Test Direction` = tail_type,
    `Group 1` = summary_stats[[group_variable]][1],
    `Group 2` = summary_stats[[group_variable]][2],
    `Median 1` = round(summary_stats$Median[1], 3),
    `Median 2` = round(summary_stats$Median[2], 3),
    `W (U) statistic` = round(test_res$statistic, 3),
    `p-value` = round(test_res$p.value, 3),
    `Alpha Level` = alpha_level,
    `Significant?` = ifelse(test_res$p.value < alpha_level, "Yes", "No"),
    `Transformed?` = transform_status
  )
  print(tidy_out)
}

# End of Script
