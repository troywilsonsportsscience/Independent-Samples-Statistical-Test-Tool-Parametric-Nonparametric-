# ===============================================================
# Independent-Samples Statistical Test Tool (Parametric + Nonparametric)
# ===============================================================
# Author: Troy Wilson
# Purpose: Automated independent-samples hypothesis testing with:
# - Group-wise normality checks (Shapiro or AD)
# - Safe transforms (log/sqrt/inv/Yeo–Johnson) applied consistently to both groups
# - Variance homogeneity (Levene’s) with Welch vs classic t selection
# - Mann–Whitney (Wilcoxon rank-sum) fallback
# - APA-style output + effect sizes (Hedges' g, Cohen's d, Glass's Δ; r for MWU)
# - Plots: histograms+densities by group, Q–Q by group, transform panels
# ===============================================================

# -------------------------------
# Helper: Install and Load Packages
# -------------------------------
required_packages <- c("tidyverse", "broom", "janitor", "ggplot2")
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, dependencies = TRUE)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
invisible(lapply(required_packages, install_if_missing))

# Optional (for AD normality; Levene; Yeo–Johnson; effect-size CI)
optional_pkgs <- c("nortest", "car", "bestNormalize", "MBESS")
invisible(lapply(optional_pkgs, function(p) if (!requireNamespace(p, quietly = TRUE)) {
  try(install.packages(p), silent = TRUE)
}))

# ===============================================================
# 0) Settings (edit as needed)
alpha_level <- 0.05
tail_type   <- "two.sided"  # "two.sided", "greater", "less"
metric_label <- "outcome"

# Controls
try_transform          <- TRUE
allowed_transforms     <- c("log","sqrt","inv","yeojohnson")
method_normality       <- "shapiro"  # "shapiro", "ad", "none"
use_levene             <- TRUE       # if FALSE (or car not available), default to Welch
show_plots             <- TRUE
show_transform_panel   <- TRUE
verbose                <- TRUE
report_cis_effect_size <- FALSE      # requires MBESS

# ===============================================================
# 1) Normality + Levene helpers
safe_sw <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  suppressWarnings(shapiro.test(x)$p.value)
}
safe_ad <- function(x) {
  if (!requireNamespace("nortest", quietly = TRUE)) return(NA_real_)
  x <- x[is.finite(x)]
  if (length(x) < 8) return(NA_real_)
  suppressWarnings(nortest::ad.test(x)$p.value)
}
norm_fun <- switch(method_normality,
                   "shapiro" = safe_sw,
                   "ad"      = safe_ad,
                   "none"    = function(x) NA_real_,
                   safe_sw)

safe_levene <- function(y, g) {
  if (!use_levene) return(NA_real_)
  if (!requireNamespace("car", quietly = TRUE)) return(NA_real_)
  df <- data.frame(y = y, g = as.factor(g))
  p <- tryCatch({
    car::leveneTest(y ~ g, data = df)$`Pr(>F)`[1]
  }, error = function(e) NA_real_)
  as.numeric(p)
}

# ===============================================================
# 2) Safe transforms (applied identically to both groups)
# We compute a global shift/epsilon using the full vector so both groups share the same mapping.
safe_log  <- function(x) {
  s <- -min(x, na.rm = TRUE) + 1e-6
  list(v = log(x + s), name = "log", s = s, e = 0, info = NULL)
}
safe_sqrt <- function(x) {
  s <- -min(x, na.rm = TRUE); if (s < 0) s <- 0
  list(v = sqrt(x + s), name = "sqrt", s = s, e = 0, info = NULL)
}
safe_inv  <- function(x) {
  s <- -min(x, na.rm = TRUE); if (s < 0) s <- 0
  e <- 1e-6
  list(v = 1/(x + s + e), name = "inv", s = s, e = e, info = NULL)
}
safe_yj <- function(x) {
  if (!requireNamespace("bestNormalize", quietly = TRUE)) return(NULL)
  obj <- bestNormalize::yeojohnson(x)
  list(v = obj$x.t, name = "yeojohnson", s = 0, e = 0, info = list(object = obj))
}

# ===============================================================
# 3) Example data (replace with your own)

# --- Load CSV and force expected column names (id, group, metric) ---
library(readr); library(dplyr); library(janitor)

path <- file.choose()
df_raw <- read_csv(path, show_col_types = FALSE) %>% clean_names()

if (ncol(df_raw) < 3) stop("CSV must have at least 3 columns: id, group, metric")

df <- df_raw %>%
  rename(
    subject_id  = 1,
    group       = 2,
    test_metric = 3
  ) %>%
  select(subject_id, group, test_metric) %>%
  drop_na(group, test_metric) %>%
  mutate(group = as.factor(group))

test_data <- df

# --- Or use built-in example ---
#set.seed(7)
# test_data <- tibble(
#   subject_id = 1:40,
#   group      = factor(rep(c("A","B"), each = 20)),
#   test_metric= c(rlnorm(20, meanlog = log(60), sdlog = 0.3),
#                  rlnorm(20, meanlog = log(72), sdlog = 0.3))
# )

# Required columns check
stopifnot(all(c("subject_id","group","test_metric") %in% names(test_data)))
test_data <- test_data %>% drop_na(group, test_metric) %>% clean_names()

# Ensure exactly two groups
if (n_distinct(test_data$group) != 2) stop("`group` must have exactly two levels.")
test_data$group <- droplevels(test_data$group)
g_levels <- levels(test_data$group)

# Split
x1 <- test_data %>% filter(group == g_levels[1]) %>% pull(test_metric)
x2 <- test_data %>% filter(group == g_levels[2]) %>% pull(test_metric)
n1 <- length(x1); n2 <- length(x2)
if (n1 < 3 || n2 < 3) stop("Each group needs at least 3 observations.")

# ===============================================================
# 4) Raw diagnostics: normality per group
p_norm_1 <- norm_fun(x1)
p_norm_2 <- norm_fun(x2)
norm_pass_raw <- if (method_normality == "none") TRUE else
  (is.finite(p_norm_1) && p_norm_1 >= alpha_level && is.finite(p_norm_2) && p_norm_2 >= alpha_level)

if (verbose && method_normality != "none") {
  cat(sprintf("Normality (raw): %s p=%.4f, %s p=%.4f → %s\n",
              g_levels[1], ifelse(is.na(p_norm_1), NA, p_norm_1),
              g_levels[2], ifelse(is.na(p_norm_2), NA, p_norm_2),
              ifelse(norm_pass_raw, "✅ pass both", "❌ at least one fails")))
}

# ===============================================================
# 5) Visualization (raw)
if (show_plots) {
  p_hist <-
    ggplot(test_data, aes(x = test_metric, fill = group)) +
    geom_histogram(aes(y = ..density..), bins = 15, alpha = 0.6, color = "black", position = "identity") +
    geom_density(alpha = 0.3) +
    facet_wrap(~ group, scales = "free") +
    labs(title = "Histograms by Group (Raw Data)", x = metric_label, y = "Density") +
    theme_minimal()
  print(p_hist)
  
  p_qq <-
    test_data %>%
    group_by(group) %>%
    summarise(val = list(test_metric), .groups = "drop") %>%
    tidyr::unnest_longer(val) %>%
    ggplot(aes(sample = val)) +
    stat_qq() + stat_qq_line(linetype = "dashed") +
    facet_wrap(~ group, scales = "free") +
    labs(title = "Q–Q Plots by Group (Raw)", x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_minimal()
  print(p_qq)
}

# ===============================================================
# 6) Transform audit & selection (if needed)
analysis_scale <- "raw"  # or name of transform
transform_meta <- list(name = "none", shift = 0, eps = 0, info = NULL,
                       norm_p1 = p_norm_1, norm_p2 = p_norm_2)

x1_use <- x1; x2_use <- x2

if (!norm_pass_raw && try_transform) {
  if (verbose) cat("↪ Trying transforms to satisfy group-wise normality...\n")
  
  # Candidate list (operate on combined vector to share the same mapping)
  X_all <- c(x1, x2)
  
  cand <- list()
  if ("log" %in% allowed_transforms)  cand$log  <- safe_log
  if ("sqrt" %in% allowed_transforms) cand$sqrt <- safe_sqrt
  if ("inv" %in% allowed_transforms)  cand$inv  <- safe_inv
  if ("yeojohnson" %in% allowed_transforms) cand$yeojohnson <- safe_yj
  
  tf_stats <- list()
  for (nm in names(cand)) {
    out <- cand[[nm]](X_all)
    if (is.null(out)) next
    v_all <- out$v
    v1 <- v_all[seq_len(n1)]
    v2 <- v_all[(n1+1):(n1+n2)]
    
    p1 <- norm_fun(v1); p2 <- norm_fun(v2)
    tf_stats[[nm]] <- list(name = nm, v1 = v1, v2 = v2, p1 = p1, p2 = p2,
                           s = out$s, e = out$e, info = out$info)
  }
  
  # Choose the transform that maximizes the *minimum* group p-value above alpha
  viable <- Filter(function(z) !is.na(z$p1) && !is.na(z$p2) && z$p1 >= alpha_level && z$p2 >= alpha_level, tf_stats)
  if (length(viable) > 0) {
    best_ix <- which.max(sapply(viable, function(z) min(z$p1, z$p2)))
    best <- viable[[best_ix]]
    x1_use <- best$v1; x2_use <- best$v2
    analysis_scale <- best$name
    transform_meta <- list(name = best$name, shift = best$s, eps = best$e,
                           info = best$info, norm_p1 = best$p1, norm_p2 = best$p2)
    if (verbose) cat(sprintf("✔️ Using %s transform (normality: %s p=%.4f; %s p=%.4f).\n",
                             best$name, g_levels[1], best$p1, g_levels[2], best$p2))
  } else {
    if (verbose) cat("❌ No transform yielded normality in both groups → will use nonparametric.\n")
  }
  
  # Optional: transform panels for comparison
  if (show_plots && show_transform_panel && length(tf_stats) > 0) {
    comp_df <- bind_rows(lapply(tf_stats, function(z) {
      tibble(
        group = rep(g_levels, times = c(length(z$v1), length(z$v2))),
        val   = c(z$v1, z$v2),
        scale = paste0(stringr::str_to_title(z$name),
                       "  (", g_levels[1], " p=", ifelse(is.finite(z$p1), sprintf("%.3f", z$p1), "NA"),
                       "; ", g_levels[2], " p=", ifelse(is.finite(z$p2), sprintf("%.3f", z$p2), "NA"), ")")
      )
    }))
    if (nrow(comp_df) > 0) {
      p_tf_hist <-
        ggplot(comp_df, aes(x = val, fill = group)) +
        geom_histogram(aes(y = ..density..), bins = 15, alpha = 0.6, color = "black", position = "identity") +
        geom_density(alpha = 0.3) +
        facet_wrap(~ scale, scales = "free") +
        labs(title = "Transform Candidates: Histograms by Group", x = "Value (scale-specific)", y = "Density") +
        theme_minimal()
      print(p_tf_hist)
      
      p_tf_qq <-
        ggplot(comp_df, aes(sample = val)) +
        stat_qq() + stat_qq_line(linetype = "dashed") +
        facet_grid(scale ~ group, scales = "free") +
        labs(title = "Transform Candidates: Q–Q by Group", x = "Theoretical Quantiles", y = "Sample Quantiles") +
        theme_minimal()
      print(p_tf_qq)
    }
  }
}

# ===============================================================
# 7) Final normality/Levene decision + test selection
use_parametric <- FALSE
use_mwu <- FALSE
equal_var <- NA

if (analysis_scale != "raw") {
  use_parametric <- TRUE
} else if (norm_pass_raw) {
  use_parametric <- TRUE
} else if (!try_transform) {
  use_parametric <- FALSE
} else {
  # No successful transform
  use_parametric <- FALSE
}

# Levene on the *final* scale if parametric; otherwise skip
levene_p <- NA_real_
if (use_parametric) {
  y_all <- c(x1_use, x2_use)
  g_all <- factor(rep(g_levels, times = c(n1, n2)), levels = g_levels)
  levene_p <- safe_levene(y_all, g_all)
  if (!is.na(levene_p)) {
    equal_var <- levene_p >= alpha_level
  } else {
    # If Levene unavailable, prefer Welch (robust) by default
    equal_var <- FALSE
  }
} else {
  use_mwu <- TRUE
}

# ===============================================================
# 8) Run tests
t_out <- NULL; w_out <- NULL
method_label <- NULL

if (use_parametric) {
  # Welch by default if unequal variances, else classic pooled t
  if (isTRUE(equal_var)) {
    t_out <- t.test(x1_use, x2_use, alternative = tail_type, var.equal = TRUE)
    method_label <- "Student t-test (equal variances)"
  } else {
    t_out <- t.test(x1_use, x2_use, alternative = tail_type, var.equal = FALSE)
    method_label <- "Welch t-test (unequal variances)"
  }
} else {
  # Mann–Whitney U (Wilcoxon rank-sum)
  # Use exact=TRUE only when small n and no ties
  any_ties <- any(duplicated(c(x1, x2)))
  exact_ok <- (n1 * n2 <= 50) && !any_ties
  w_out <- wilcox.test(x1, x2, alternative = tail_type, exact = exact_ok, correct = !exact_ok)
  method_label <- "Mann–Whitney U (Wilcoxon rank-sum)"
}

# ===============================================================
# 9) Descriptives & Effect sizes
desc <- tibble(
  group = g_levels,
  Mean  = c(mean(x1, na.rm = TRUE), mean(x2, na.rm = TRUE)),
  SD    = c(sd(x1, na.rm = TRUE),   sd(x2, na.rm = TRUE)),
  Median= c(median(x1, na.rm = TRUE), median(x2, na.rm = TRUE)),
  n     = c(n1, n2)
)

# Effect size helpers
J_correction <- function(n1, n2) 1 - 3/(4*(n1 + n2) - 9)  # unbiased g correction
hedges_g_ci <- function(d, n1, n2, conf.level = 0.95) {
  if (!requireNamespace("MBESS", quietly = TRUE)) return(NA_character_)
  out <- try(MBESS::ci.smd(smd = d, n.1 = n1, n.2 = n2, conf.level = conf.level, Unbiased = TRUE), silent = TRUE)
  if (inherits(out, "try-error")) return(NA_character_)
  sprintf("[%.3f, %.3f]", out$Lower.Conf.Limit.smd, out$Upper.Conf.Limit.smd)
}

effects <- list()

if (use_parametric) {
  m1 <- mean(x1_use); m2 <- mean(x2_use)
  s1 <- sd(x1_use);   s2 <- sd(x2_use)
  
  # Pooled SD for d/g (still commonly reported even with Welch)
  sp <- sqrt(((n1 - 1)*s1^2 + (n2 - 1)*s2^2) / (n1 + n2 - 2))
  d  <- (m1 - m2) / sp
  g  <- J_correction(n1, n2) * d
  
  effects$cohen_d  <- round(d, 3)
  effects$hedges_g <- round(g, 3)
  effects$glass_delta_1 <- round((m1 - m2) / s1, 3) # using grp1 SD
  effects$glass_delta_2 <- round((m1 - m2) / s2, 3) # using grp2 SD
  if (report_cis_effect_size) {
    effects$hedges_g_95CI <- hedges_g_ci(d, n1, n2)
  }
} else {
  # r = Z / sqrt(N). Compute Z from U (approximation).
  # Extract U from wilcox.test: statistic is W, which equals sum of ranks for first group.
  # Convert to U1 = W - n1*(n1+1)/2
  W  <- as.numeric(w_out$statistic)
  U1 <- W - n1*(n1 + 1)/2
  mu_U <- n1 * n2 / 2
  sigma_U <- sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
  # Continuity correction if not using exact
  cc <- if (!isTRUE(attr(w_out, "exact"))) 0.5 * sign(U1 - mu_U) else 0
  z  <- (U1 - mu_U - cc) / sigma_U
  r  <- as.numeric(z / sqrt(n1 + n2))
  effects$U <- round(U1, 2)
  effects$effect_r <- round(r, 3)
}

# ===============================================================
# 10) APA-style output
if (verbose) {
  cat("\n---- Summary Statistics ----\n")
  print(desc)
  
  # Diagnostics text
  norm_txt <- if (method_normality == "none") {
    "Normality not assessed by request."
  } else {
    sprintf("Normality: %s p=%.3f; %s p=%.3f.",
            g_levels[1], ifelse(is.finite(p_norm_1), p_norm_1, NA),
            g_levels[2], ifelse(is.finite(p_norm_2), p_norm_2, NA))
  }
  
  if (use_parametric) {
    cat("\n--- APA (t-test) ---\n")
    cat(sprintf("%s on %s (α = %.2f): t(%d) = %.2f, p = %.4f, 95%% CI [% .2f, % .2f].\n",
                method_label, if (analysis_scale=="raw") "raw data" else paste0(analysis_scale, " transform"),
                alpha_level, round(unname(t_out$parameter)), unname(t_out$statistic),
                t_out$p.value, t_out$conf.int[1], t_out$conf.int[2]))
    cat(sprintf("Effect sizes: Hedges' g = %.3f%s; Cohen's d = %.3f; Glass's Δ₁ = %.3f, Δ₂ = %.3f.\n",
                effects$hedges_g,
                if (!is.null(effects$hedges_g_95CI)) paste0(" (95% CI ", effects$hedges_g_95CI, ")") else "",
                effects$cohen_d, effects$glass_delta_1, effects$glass_delta_2))
    if (!is.na(levene_p)) {
      cat(sprintf("Levene’s test: p = %.3f → %s variances.\n",
                  levene_p, ifelse(levene_p >= alpha_level, "assumed equal", "Welch correction used")))
    } else {
      cat("Levene’s test unavailable → defaulting to Welch.\n")
    }
    if (analysis_scale == "raw") {
      cat(sprintf("Assumptions: %s\n", norm_txt))
    } else {
      cat(sprintf("Assumptions: Raw %s. Transform (%s) improved normality: %s p=%.3f; %s p=%.3f.\n",
                  "violated at least one group",
                  analysis_scale, g_levels[1], transform_meta$norm_p1, g_levels[2], transform_meta$norm_p2))
    }
  } else {
    cat("\n--- APA (Mann–Whitney U) ---\n")
    cat(sprintf("Mann–Whitney U test: U = %.2f, p = %.4f; effect size r = %.3f.\n",
                effects$U, w_out$p.value, effects$effect_r))
    cat(sprintf("Assumptions: %s Transforms %s restore normality in both groups → used nonparametric.\n",
                norm_txt, if (try_transform) "did not" else "were not attempted to"))
  }
  
  # Final quick interpretation
  sig_txt <- if ((use_parametric && t_out$p.value < alpha_level) ||
                 (!use_parametric && w_out$p.value < alpha_level)) "significant" else "not significant"
  cat(sprintf("\nInterpretation: Group difference was %s at α = %.2f.\n", sig_txt, alpha_level))
}

# ===============================================================
# 11) Final-scale plots (optional)
if (show_plots) {
  df_final <- tibble(
    group = rep(g_levels, times = c(n1, n2)),
    val   = c(x1_use, x2_use)
  )
  lbl <- if (analysis_scale=="raw") "raw" else paste0(analysis_scale, "-transformed")
  
  p_final_hist <-
    ggplot(df_final, aes(x = val, fill = group)) +
    geom_histogram(aes(y = ..density..), bins = 15, alpha = 0.6, color = "black", position = "identity") +
    geom_density(alpha = 0.3) +
    facet_wrap(~ group, scales = "free") +
    labs(title = paste("Final Scale Histograms by Group (", lbl, ")", sep=""),
         x = if (analysis_scale=="raw") metric_label else paste(metric_label, "(", lbl, ")", sep=" "),
         y = "Density") +
    theme_minimal()
  print(p_final_hist)
  
  p_final_qq <-
    df_final %>%
    group_by(group) %>% summarise(val = list(val), .groups = "drop") %>%
    tidyr::unnest_longer(val) %>%
    ggplot(aes(sample = val)) +
    stat_qq() + stat_qq_line(linetype = "dashed") +
    facet_wrap(~ group, scales = "free") +
    labs(title = paste0("Q–Q Plots by Group (", lbl, ")"),
         x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_minimal()
  print(p_final_qq)
}

# ===============================================================
# 12) Tidy result table
tidy_out <-
  if (use_parametric) {
    tibble::tibble(
      Method        = method_label,
      `Test Type`   = ifelse(analysis_scale == "raw", "t (raw)", paste0("t (", analysis_scale, ")")),
      Direction     = tail_type,
      Group1        = g_levels[1],
      Group2        = g_levels[2],
      n1            = n1,
      n2            = n2,
      Mean1_raw     = round(mean(x1, na.rm = TRUE), 3),
      SD1_raw       = round(sd(x1,   na.rm = TRUE), 3),
      Mean2_raw     = round(mean(x2, na.rm = TRUE), 3),
      SD2_raw       = round(sd(x2,   na.rm = TRUE), 3),
      `t(df)`       = paste0("t(", round(unname(t_out$parameter)), ") = ", round(unname(t_out$statistic), 3)),
      `p-value`     = round(t_out$p.value, 5),
      `95% CI`      = paste0("[", round(t_out$conf.int[1], 3), ", ", round(t_out$conf.int[2], 3), "]"),
      Levene_p      = ifelse(is.na(levene_p), NA, round(levene_p, 5)),
      Transform     = analysis_scale,
      Shift         = transform_meta$shift,
      Eps           = transform_meta$eps,
      Cohen_d       = effects$cohen_d,
      Hedges_g      = effects$hedges_g,
      Hedges_g_95CI = if (!is.null(effects$hedges_g_95CI)) effects$hedges_g_95CI else NA_character_,
      Glass_Delta1  = effects$glass_delta_1,
      Glass_Delta2  = effects$glass_delta_2
    )
  } else {
    tibble::tibble(
      Method        = method_label,
      Direction     = tail_type,
      Group1        = g_levels[1],
      Group2        = g_levels[2],
      n1            = n1,
      n2            = n2,
      U             = effects$U,
      `p-value`     = round(w_out$p.value, 5),
      Effect_r      = effects$effect_r
    )
  }

if (verbose) {
  cat("\n---- Tidy Output ----\n")
  print(tidy_out)
}
