# ===============================================================
# Independent-Samples Statistical Test Tool (Parametric + Nonparametric)
# ===============================================================
# - Checks normality (per group) & variance equality (Levene/Brown–Forsythe)
# - Chooses classic t, Welch t, or Mann–Whitney automatically
# - Optional transformations to rescue normality (safe log/sqrt/inv)
# - Outputs APA-style text + tidy summary
# - Effect sizes:
#     * Parametric: Cohen's d (pooled), Hedges' g, Welch's d; (optional) Glass's Δ when heteroscedastic
#     * Nonparametric: Rank-biserial r, Cliff's δ (optional CI if effsize is available)
# - Returns a structured result object invisibly for programmatic use
# ===============================================================

# -------------------------------
# Helper: Install and Load Packages
# -------------------------------
required_packages <- c("tidyverse", "car", "broom", "janitor", "rstatix")
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, dependencies = TRUE)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
invisible(lapply(required_packages, install_if_missing))
# Optional packages used if available (no hard failure if absent)
optional_packages <- c("nortest", "MBESS", "effsize")
invisible(lapply(optional_packages, function(p) if (!requireNamespace(p, quietly = TRUE)) try(suppressMessages(install.packages(p)), silent = TRUE)))

# ---- 0. Settings (edit for your analysis) ----
alpha_level     <- 0.05                # significance level
tail_type       <- "two.sided"         # "two.sided", "greater", "less"
metric_label    <- "exam score"        # friendly label for reporting
group_variable  <- "group"             # grouping column name (factor with 2 levels)
outcome_variable<- "score"             # numeric outcome column name

# Controls
try_transform   <- TRUE                # if TRUE and normality fails → attempt transforms
allowed_transforms <- c("log","sqrt","inv")  # any subset of "log","sqrt","inv"
show_intermediate_plots <- TRUE        # show hist/box plots & diagnostics?
verbose         <- TRUE                # print tables and messages
method_normality <- "shapiro"          # "shapiro", "ad", "none"
report_cis_effect_sizes <- TRUE        # requires MBESS (parametric) / effsize (nonparametric) if available

# (Optional) Effect-size option for unequal variances: choose a reference group for Glass's Δ
use_glass_delta <- TRUE
glass_reference_group <- "A"           # which group's SD to use if Glass's Δ is reported

# ---- 1. Data ----
path <- file.choose()  # Select CSV file
df <- read_csv(path, show_col_types = FALSE) %>% clean_names() # Convert to snake_case for ease of use
test_data <- df %>%
  select(subject_id, test_metric)  # or rename columns to match this format
                 
# ---- Sanity checks ----
stopifnot(all(c(group_variable, outcome_variable) %in% names(test_data)))
test_data <- test_data %>% select(all_of(c(group_variable, outcome_variable))) %>% drop_na()
test_data[[group_variable]] <- factor(test_data[[group_variable]])
if (nlevels(test_data[[group_variable]]) != 2) stop("This tool expects exactly 2 groups.")
# Lock a consistent reference for effect-size sign (first level)
test_data[[group_variable]] <- forcats::fct_relevel(test_data[[group_variable]], levels(test_data[[group_variable]])[1])
reference_group <- levels(test_data[[group_variable]])[1]

# ---- 2. Assumption Checks ----
if (verbose) cat("\n---- Assumption Check: Normality (per group) ----\n")

# Normality helpers
safe_shapiro <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  shapiro.test(x)$p.value
}
safe_ad <- function(x) {
  if (!requireNamespace("nortest", quietly = TRUE)) return(NA_real_)
  x <- x[is.finite(x)]
  if (length(x) < 8) return(NA_real_) # AD prefers larger n
  nortest::ad.test(x)$p.value
}

norm_fun <- switch(method_normality,
                   "shapiro" = safe_shapiro,
                   "ad"      = safe_ad,
                   "none"    = function(x) NA_real_,
                   safe_shapiro)

norm_results <- test_data %>%
  group_by(!!sym(group_variable)) %>%
  summarise(p_value = norm_fun(!!sym(outcome_variable)), .groups="drop") %>%
  rename(group = !!sym(group_variable))

if (verbose && method_normality != "none") print(norm_results)

if (method_normality == "none") {
  normality_ok <- TRUE
} else {
  if (verbose) {
    for (i in seq_len(nrow(norm_results))) {
      g <- as.character(norm_results$group[i])
      p <- norm_results$p_value[i]
      flag <- ifelse(is.na(p), "⚠️ insufficient n for test", ifelse(p > alpha_level, "✅ Normal", "❌ Non-normal"))
      cat(sprintf("Group %s: normality p = %s → %s\n", g, ifelse(is.na(p),"NA",sprintf("%.4f",p)), flag))
    }
    cat("\nInterpretation: p >", alpha_level, "→ normality ok; p <=", alpha_level, "→ possible non-normality.\n")
  }
  normality_ok <- all(is.na(norm_results$p_value) | norm_results$p_value > alpha_level)
}
normality_ok_original <- normality_ok

# Optional: quick histograms by group (raw)
if (show_intermediate_plots) {
  p_raw <-
    ggplot(test_data, aes_string(x = outcome_variable, fill = group_variable)) +
    geom_histogram(aes(y = ..density..), bins = 15, alpha = 0.6, color = "black", position = "identity") +
    geom_density(alpha = 0.3) +
    facet_wrap(as.formula(paste("~", group_variable))) +
    labs(title = "Histograms by Group (Raw Data)", x = metric_label, y = "Density") +
    theme_gray()
  print(p_raw)
}

# ---- 2A. Optional Transformation Attempt ----
# Safer transforms (record shift; prevent division by zero)
safe_log  <- function(x) { s <- -min(x, na.rm=TRUE) + 1e-6; attr(x,"shift") <- s; log(x + s) }
safe_sqrt <- function(x) sqrt(pmax(x, 0))
safe_inv  <- function(x) {
  s <- if (any(x <= 0, na.rm=TRUE)) -min(x, na.rm=TRUE) + 1e-6 else 0
  eps <- 1e-6
  attr(x,"shift") <- s
  1 / (x + s + eps)
}

if (!normality_ok && try_transform) {
  if (verbose) cat("\n↪ Normality violation detected. Attempting transformations...\n")
  available <- list(log = safe_log, sqrt = safe_sqrt, inv = safe_inv)
  selected  <- available[names(available) %in% allowed_transforms]
  if (length(selected) == 0) {
    if (verbose) cat("⚠️ No valid transformations listed in `allowed_transforms`. Skipping.\n")
  } else {
    tf_results <- lapply(names(selected), function(tf) {
      v <- selected[[tf]](test_data[[outcome_variable]])
      if (!all(is.finite(v))) return(list(name=tf, ok=FALSE, min_p=NA, vec=NULL))
      tmp <- test_data %>% mutate(.tf = v)
      ptab <- tmp %>%
        group_by(!!sym(group_variable)) %>%
        summarise(p_value = norm_fun(.tf), .groups="drop")
      ok <- if (method_normality == "none") TRUE else all(is.na(ptab$p_value) | ptab$p_value > alpha_level)
      list(name=tf, ok=ok, min_p=suppressWarnings(min(ptab$p_value, na.rm=TRUE)), vec=v)
    })
    viable <- Filter(function(x) isTRUE(x$ok), tf_results)
    if (length(viable) > 0) {
      best <- viable[[which.max(sapply(viable, function(x) x$min_p))]]
      test_data[[outcome_variable]] <- best$vec
      normality_ok <- TRUE
      assign("best_transform", best$name, inherits = TRUE)
      if (verbose) cat("✔️ Transformation successful using:", best$name, "→ proceeding with parametric testing.\n")
    } else {
      if (verbose) cat("❌ No transformation restored normality. Proceeding with nonparametric test.\n")
    }
  }
}

# ---- 2B. Homogeneity of Variance (Levene/Brown–Forsythe) ----
if (normality_ok) {
  if (verbose) cat("\n---- Assumption Check: Homogeneity of Variance (Levene; median-centered) ----\n")
  formula_obj <- reformulate(group_variable, outcome_variable)  # outcome ~ group
  levene_res  <- car::leveneTest(formula_obj, data = test_data, center = median)
  if (verbose) print(levene_res)
  p_var <- levene_res$`Pr(>F)`[1]
  variance_ok <- p_var > alpha_level
  if (verbose) cat(sprintf("\nLevene (median): F(%d,%d)=%.3f, p=%.4f → %s variances.\n",
                           levene_res$Df[1], levene_res$Df[2],
                           levene_res$`F value`[1], p_var,
                           ifelse(variance_ok,"equal","unequal")))
} else {
  if (verbose) cat("\n⏭ Skipping Levene's test (normality failed → using nonparametric test).\n")
  variance_ok <- FALSE
}

# ---- 2C. Visual Inspection (final variable state) ----
if (show_intermediate_plots) {
  state_lab <- if (exists("best_transform")) paste0("Transformed (", best_transform, ")") else "Raw"
  
  p_final_hist <-
    ggplot(test_data, aes_string(x = outcome_variable, fill = group_variable)) +
    geom_histogram(aes(y = ..density..), bins = 15, alpha = 0.6, color = "black", position = "identity") +
    geom_density(alpha = 0.3) +
    facet_wrap(as.formula(paste("~", group_variable))) +
    labs(title = "Final Variable by Group", subtitle = paste("Data used:", state_lab),
         x = metric_label, y = "Density") +
    theme_gray()
  print(p_final_hist)
  
  p_final_box <-
    ggplot(test_data, aes_string(x = group_variable, y = outcome_variable, fill = group_variable)) +
    geom_boxplot(alpha = 0.7, outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
    labs(title = "Boxplots (Final Variable)", subtitle = paste("Data used:", state_lab),
         y = metric_label, x = "Group") +
    theme_gray()
  print(p_final_box)
}

# ---- 3. Select & Run Test ----
if (verbose) cat("\n---- Section 3: Selecting and Running the Appropriate Test ----\n")
if (normality_ok && variance_ok) {
  test_method <- "classic_t"
  test_res <- t.test(reformulate(group_variable, outcome_variable), data = test_data,
                     var.equal = TRUE, alternative = tail_type)
  if (verbose) cat("✅ Normality & equal variances → Classic independent-samples t-test.\n")
} else if (normality_ok && !variance_ok) {
  test_method <- "welch_t"
  test_res <- t.test(reformulate(group_variable, outcome_variable), data = test_data,
                     var.equal = FALSE, alternative = tail_type)
  if (verbose) cat("⚠️ Normality ok but variances unequal → Welch’s t-test.\n")
} else {
  test_method <- "nonparametric"
  # Build vectors to decide exact vs approx and handle ties
  lvl <- levels(test_data[[group_variable]])
  x <- test_data %>% filter(!!sym(group_variable) == lvl[1]) %>% pull(!!sym(outcome_variable))
  y <- test_data %>% filter(!!sym(group_variable) == lvl[2]) %>% pull(!!sym(outcome_variable))
  any_ties <- any(duplicated(c(x, y)))
  exact_ok <- (length(x) * length(y) <= 5e4) && !any_ties
  test_res <- wilcox.test(x, y, alternative = tail_type, exact = exact_ok, correct = !exact_ok)
  if (verbose) cat(sprintf("❗ Normality failed → Mann–Whitney U (exact=%s).\n", ifelse(exact_ok,"TRUE","FALSE")))
}

# ---- 4. Descriptives & Effect Sizes ----
if (verbose) cat("\n---- Section 4: Summary Statistics & Effect Sizes ----\n")
summary_stats <- test_data %>%
  group_by(!!sym(group_variable)) %>%
  summarise(Mean = mean(!!sym(outcome_variable)),
            SD = sd(!!sym(outcome_variable)),
            Median = median(!!sym(outcome_variable)),
            n = n(), .groups="drop")
if (verbose) print(summary_stats)

# Convenience
g1 <- as.character(summary_stats[[group_variable]][1])
g2 <- as.character(summary_stats[[group_variable]][2])
m1 <- summary_stats$Mean[1]; s1 <- summary_stats$SD[1]; n1 <- summary_stats$n[1]
m2 <- summary_stats$Mean[2]; s2 <- summary_stats$SD[2]; n2 <- summary_stats$n[2]

# Effect size helpers
pooled_sd <- function(s1, s2, n1, n2) sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1 + n2 - 2))
welch_sd  <- function(s1, s2, n1, n2) sqrt((s1^2 / n1) + (s2^2 / n2))
hedges_g  <- function(d, n1, n2) { J <- 1 - 3/(4*(n1 + n2) - 9); J * d }

es_out <- list()
if (test_method %in% c("classic_t","welch_t")) {
  sp <- pooled_sd(s1, s2, n1, n2)
  sw <- welch_sd(s1, s2, n1, n2)
  d_pooled <- (m1 - m2) / sp
  d_welch  <- (m1 - m2) / sw
  g_pooled <- hedges_g(d_pooled, n1, n2)
  
  es_out$cohens_d_pooled <- round(d_pooled, 3)
  es_out$hedges_g        <- round(g_pooled, 3)
  es_out$cohens_d_welch  <- round(d_welch, 3)
  
  if (report_cis_effect_sizes && requireNamespace("MBESS", quietly = TRUE)) {
    ci_g <- try(MBESS::ci.smd(smd = d_pooled, n.1 = n1, n.2 = n2, conf.level = 0.95, Unbiased = TRUE), silent = TRUE)
    if (!inherits(ci_g, "try-error")) {
      es_out$hedges_g_CI <- paste0("[", round(ci_g$Lower.Conf.Limit.smd, 3), ", ", round(ci_g$Upper.Conf.Limit.smd, 3), "]")
    }
  }
  
  if (!variance_ok && use_glass_delta) {
    s_ref <- if (g1 == glass_reference_group) s1 else if (g2 == glass_reference_group) s2 else s1
    es_out$glass_delta <- round((m1 - m2)/s_ref, 3)
  }
} else {
  # Nonparametric: U, rank-biserial r, Cliff’s delta (lock sign to reference_group = first level)
  W <- as.numeric(test_res$statistic)
  U1 <- W - n1*(n1+1)/2
  r_rb <- 1 - (2*U1)/(n1*n2)
  cliff_delta <- (2*U1)/(n1*n2) - 1
  
  # Ensure effect > 0 means reference_group (g1) shows larger values
  # (g1 is already the first level / reference by construction above)
  es_out$rank_biserial_r <- round(r_rb, 3)
  es_out$cliffs_delta    <- round(cliff_delta, 3)
  
  if (report_cis_effect_sizes && requireNamespace("effsize", quietly = TRUE)) {
    cd <- try(effsize::cliff.delta(
      x = test_data %>% filter(!!sym(group_variable) == g1) %>% pull(!!sym(outcome_variable)),
      y = test_data %>% filter(!!sym(group_variable) == g2) %>% pull(!!sym(outcome_variable)),
      conf.level = 0.95
    ), silent = TRUE)
    if (!inherits(cd, "try-error") && !is.null(cd$conf.int)) {
      es_out$cliffs_delta_CI <- paste0("[", round(cd$conf.int[1], 3), ", ", round(cd$conf.int[2], 3), "]")
    }
  }
}
if (verbose) print(es_out)

# ---- 5. APA-style Output ----
if (verbose) cat("\n---- Section 5: APA-style Output ----\n")
if (test_method != "nonparametric") {
  tidy_res <- broom::tidy(test_res)
  df     <- round(tidy_res$parameter, 1)
  tval   <- round(tidy_res$statistic, 2)
  pval   <- round(tidy_res$p.value, 4)
  ci_low <- round(tidy_res$conf.low, 2)
  ci_high<- round(tidy_res$conf.high, 2)
  
  direction_text <- switch(tail_type,
                           "two.sided" = "differ",
                           "greater"   = "is greater than",
                           "less"      = "is less than")
  
  if (exists("normality_ok_original") && !normality_ok_original && try_transform && exists("best_transform")) {
    cat(sprintf("Note: '%s' was %s-transformed to meet normality.\n", outcome_variable, best_transform))
  }
  
  test_label <- if (test_method == "welch_t") "Welch’s t-test" else "independent-samples t-test"
  cat(sprintf(
    "%s showed that %s %s %s %s, t(%.1f) = %.2f, p = %.4f, 95%% CI [%.2f, %.2f].\n",
    tools::toTitleCase(test_label),
    g1, sprintf("(M = %.2f, SD = %.2f)", m1, s1),
    direction_text,
    paste(g2, sprintf("(M = %.2f, SD = %.2f)", m2, s2)),
    df, tval, pval, ci_low, ci_high
  ))
  
  # Report ES succinctly
  cat(sprintf("Cohen’s d (pooled) = %.3f; Hedges’ g = %.3f", es_out$cohens_d_pooled, es_out$hedges_g))
  if (!is.null(es_out$hedges_g_CI)) cat(sprintf(" (95%% CI %s)", es_out$hedges_g_CI))
  cat(sprintf("; Welch’s d = %.3f", es_out$cohens_d_welch))
  if (!variance_ok && use_glass_delta && !is.null(es_out$glass_delta)) {
    cat(sprintf("; Glass’s Δ (ref %s) = %.3f", glass_reference_group, es_out$glass_delta))
    cat("\nNote: With unequal variances, Glass’s Δ or Welch’s d is preferable to pooled-SD d.\n")
  } else cat("\n")
  
} else {
  stat <- as.numeric(test_res$statistic)
  pval <- round(test_res$p.value, 4)
  cat(sprintf(
    "A Mann–Whitney U test indicated a %s difference in %s between groups, W = %.2f, p = %.4f.\n",
    ifelse(pval < alpha_level, "significant", "non-significant"),
    metric_label, stat, pval
  ))
  cat(sprintf("Medians: %s = %.2f, %s = %.2f.\n", g1, summary_stats$Median[1], g2, summary_stats$Median[2]))
  cat(sprintf("Effect sizes: Rank-biserial r = %.3f; Cliff’s δ = %.3f",
              es_out$rank_biserial_r, es_out$cliffs_delta))
  if (!is.null(es_out$cliffs_delta_CI)) cat(sprintf(" (95%% CI %s)", es_out$cliffs_delta_CI))
  cat(".\n")
}

# ---- 6. Tidy Output Table ----
if (verbose) cat("\n---- Section 6: Tidy Output Table ----\n")
transform_status <- if (exists("best_transform") && !normality_ok_original && try_transform) best_transform else "None"

if (test_method != "nonparametric") {
  tidy_out <- tibble::tibble(
    Method = test_res$method,
    `Test Type` = ifelse(test_method == "classic_t", "Independent t-test", "Welch's t-test"),
    `Test Direction` = tail_type,
    `Group 1` = g1, `Group 2` = g2,
    `Mean 1` = round(m1, 3), `SD 1` = round(s1, 3),
    `Mean 2` = round(m2, 3), `SD 2` = round(s2, 3),
    `t(df)` = paste0("t(", round(test_res$parameter, 1), ") = ", round(test_res$statistic, 2)),
    `p-value` = round(test_res$p.value, 4),
    `95% CI` = paste0("[", round(test_res$conf.int[1], 2), ", ", round(test_res$conf.int[2], 2), "]"),
    `Transformed?` = transform_status,
    `Cohen_d_pooled` = es_out$cohens_d_pooled,
    `Hedges_g` = es_out$hedges_g,
    `Hedges_g_95CI` = if (!is.null(es_out$hedges_g_CI)) es_out$hedges_g_CI else NA_character_,
    `Cohen_d_welch` = es_out$cohens_d_welch,
    `Glass_Delta_ref` = if (!variance_ok && use_glass_delta) glass_reference_group else NA_character_,
    `Glass_Delta` = if (!variance_ok && use_glass_delta) es_out$glass_delta else NA_real_
  )
  if (verbose) print(tidy_out)
} else {
  tidy_out <- tibble::tibble(
    Method = "Mann–Whitney U (Wilcoxon Rank-Sum)",
    `Test Direction` = tail_type,
    `Group 1` = g1, `Group 2` = g2,
    `Median 1` = round(summary_stats$Median[1], 3),
    `Median 2` = round(summary_stats$Median[2], 3),
    `W (U-related)` = round(as.numeric(test_res$statistic), 3),
    `p-value` = round(test_res$p.value, 4),
    `Alpha` = alpha_level,
    `Significant?` = ifelse(test_res$p.value < alpha_level, "Yes", "No"),
    `Transformed?` = transform_status,
    `Rank-biserial_r` = es_out$rank_biserial_r,
    `Cliffs_delta` = es_out$cliffs_delta,
    `Cliffs_delta_95CI` = if (!is.null(es_out$cliffs_delta_CI)) es_out$cliffs_delta_CI else NA_character_
  )
  if (verbose) print(tidy_out)
}

# ---- 7. Return structured result (invisible) ----
result <- list(
  settings = list(
    alpha_level = alpha_level, tail_type = tail_type,
    transform = if (exists("best_transform")) best_transform else "none",
    method_normality = method_normality, report_cis_effect_sizes = report_cis_effect_sizes
  ),
  test_method = test_method,
  test_res = test_res,
  summary_stats = summary_stats,
  effect_sizes = es_out,
  tidy_out = tidy_out,
  plots = list(
    p_raw = if (exists("p_raw")) p_raw else NULL,
    p_final_hist = if (exists("p_final_hist")) p_final_hist else NULL,
    p_final_box  = if (exists("p_final_box"))  p_final_box  else NULL
  )
)
invisible(result) # Set to print for convenient data dump.

# ============================ End of Script ============================
