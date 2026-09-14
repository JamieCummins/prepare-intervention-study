#!/usr/bin/env Rscript
#
# Simulated power analysis for the PREPARE intervention study
# ============================================================
#
# Design (as preregistered, September 2026):
#   * n = 400 students, individually randomised once into groups A and B.
#   * Week-wise alternating crossover over 10 topic-weeks: group A has early
#     PREPARE access on even topic-weeks, group B on odd topic-weeks. Every
#     week therefore has a concurrent randomised control on the same topic,
#     and total access is equal across groups.
#   * Weekly rhythm: lecture Monday -> PREPARE unlock Tuesday (early-access
#     group) -> Qualtrics assessment Friday -> release to the other group.
#   * Outcome: score on the weekly assessment for that week's lecture
#     (proportion of k items correct).
#   * Primary analysis (intention-to-treat), one lesson per week so
#     lesson == topic-week, which enters as a fixed effect with sum-to-zero
#     contrasts:
#
#         score ~ access + lesson + (1 + access | student)
#
#     fitted with lmerTest (Satterthwaite df), two-sided alpha = .05, with
#     the preregistered fallback ladder (drop the intercept-slope
#     correlation, then the slope) applied whenever a fit is singular or
#     fails to converge. See power_fit.R.
#
# Generative model (item-level binomial, so ceiling effects emerge naturally):
#
#   logit P(correct)_iwj = alpha + b_i + u_w + e_iw + (delta + s_i + d_w) * access_iw
#
#     b_i  ~ N(0, sigma_student)  student ability
#     u_w  ~ N(0, sigma_lesson)   lesson/topic difficulty
#     e_iw ~ N(0, sigma_sw)       student-by-week noise (overdispersion
#                                 relative to pure binomial)
#     s_i  ~ N(0, sigma_slope)    student-varying treatment effect
#     d_w  ~ N(0, sigma_delta)    topic-varying treatment effect (default 0)
#
#   alpha is calibrated numerically so that the *marginal* control-arm
#   accuracy equals base_acc, and delta so that the marginal treated-arm
#   accuracy equals base_acc + gain (accounting for the extra treated-arm
#   variance from s_i and d_w). Effects are therefore specified in accuracy
#   percentage points, the same scale as the fitted lmer estimate.
#
#   Missingness: each student-week assessment is completed with probability
#   `completion` (MCAR). An optional `dropout` parameter additionally makes
#   a fraction of students stop responding permanently from a uniformly
#   drawn week onward (default 0).
#
# Calibration: base_acc, sigma_student, sigma_lesson and sigma_sw were
# calibrated against the HS2025 weekly assessment data
# (analysis/historical_performance.R). sigma_slope has no data anchor
# (HS2025 had no intervention) and is swept in power_sensitivity.R; the
# default here is a moderate value so that the random slope in the primary
# model corresponds to real heterogeneity rather than a zero variance.
#
# Usage:
#   Rscript preregistration/power_simulation.R                 # full grid, 500 sims/cell
#   Rscript preregistration/power_simulation.R --sims 1000     # final prereg run
#   Rscript preregistration/power_simulation.R --quick         # ~1 min smoke test
#   Options: --sims N   --cores N   --seed N   --quick
#
# Outputs (written next to this script, in power_analysis/):
#   power_raw.csv       one row per simulated trial (written incrementally)
#   power_summary.csv   power, MCSE, bias, fallback/ceiling rates per cell
#   power_mde.csv       smallest swept effect reaching 80% power per design cell
#   power_curves.png    power curves
#   sessionInfo.txt     reproducibility record

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(parallel)
})

script_dir <- (function() {
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1]))) else getwd()
})()
source(file.path(script_dir, "power_fit.R"))

## ---------------------------------------------------------------------------
## Configuration
## ---------------------------------------------------------------------------

cfg <- list(
  n_students      = 400,
  n_weeks_grid    = c(8, 10, 12),         # 10 is the preregistered design; 8 and 12
                                          # bracket it (e.g., lost weeks)
  completion_grid = c(0.30, 0.60, 0.90),  # HS2025 voluntary in-week uptake was ~4-10%;
                                          # the study needs incentivised completion
  k_items_grid    = c(5, 8, 12),          # k=8 matches the planned assessment length and
                                          # the effective length of the HS2025 quizzes
  base_acc_grid   = 0.58,                 # HS2025 in-week mean; the August run showed
                                          # baseline accuracy (0.45-0.70) barely moves power
  gain_grid       = c(0, 0.02, 0.04, 0.06), # accuracy gain in the access arm; 0 = type-I error check

  # Variance components on the log-odds scale, calibrated against HS2025
  # in-week completions (analysis/historical_performance.R, 2026-08-05).
  sigma_student = 0.6,   # empirical (logit approx of prop-scale SD 0.145)
  sigma_lesson  = 0.5,   # empirical (logit approx of prop-scale SD 0.121)
  sigma_sw      = 0.2,   # chosen so total residual matches the empirical 0.181
                         # at k = 8 (binomial item noise supplies the rest)
  sigma_slope   = 0.2,   # student-varying treatment effect SD (~5 pp); no data
                         # anchor, swept in power_sensitivity.R
  sigma_delta   = 0.0,   # topic-varying treatment effect SD (sensitivity knob)
  dropout       = 0.0,   # fraction of students who stop responding permanently

  alpha_level = 0.05,
  n_sims      = 500,
  cores       = max(1L, parallel::detectCores() - 1L),
  master_seed = 20260905
)

args <- commandArgs(trailingOnly = TRUE)
val_after <- function(flag, cast) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) cast(args[i + 1]) else NULL
}
quick <- "--quick" %in% args
if (quick) {
  cfg$n_weeks_grid    <- 10
  cfg$completion_grid <- c(0.60, 1.00)
  cfg$k_items_grid    <- 8
  cfg$gain_grid       <- c(0, 0.04)
  cfg$n_sims          <- 25
}
opt <- val_after("--sims",  as.integer); if (!is.null(opt)) cfg$n_sims      <- opt
opt <- val_after("--cores", as.integer); if (!is.null(opt)) cfg$cores       <- opt
opt <- val_after("--seed",  as.integer); if (!is.null(opt)) cfg$master_seed <- opt

out_dir <- file.path(script_dir, "power_analysis", if (quick) "quick" else ".")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------------------
## Calibration: fixed effects on the log-odds scale from marginal accuracies
## ---------------------------------------------------------------------------

# E[ plogis(a + sigma*Z) ] for Z ~ N(0,1)
marginal_acc <- function(a, sigma) {
  integrate(function(z) plogis(a + sigma * z) * dnorm(z), -Inf, Inf)$value
}
calibrate_alpha <- function(target, sigma) {
  uniroot(function(a) marginal_acc(a, sigma) - target, c(-15, 15))$root
}

sigma_ctrl  <- sqrt(cfg$sigma_student^2 + cfg$sigma_lesson^2 + cfg$sigma_sw^2)
sigma_treat <- sqrt(sigma_ctrl^2 + cfg$sigma_slope^2 + cfg$sigma_delta^2)

grid <- expand.grid(
  n_weeks    = cfg$n_weeks_grid,
  completion = cfg$completion_grid,
  k_items    = cfg$k_items_grid,
  base_acc   = cfg$base_acc_grid,
  gain       = cfg$gain_grid,
  KEEP.OUT.ATTRS = FALSE
)
grid$scenario <- seq_len(nrow(grid))

grid$alpha <- vapply(grid$base_acc, calibrate_alpha, numeric(1), sigma = sigma_ctrl)
grid$delta <- vapply(grid$base_acc + grid$gain, calibrate_alpha, numeric(1),
                     sigma = sigma_treat) - grid$alpha

## ---------------------------------------------------------------------------
## Simulation and fitting
## ---------------------------------------------------------------------------

simulate_trial <- function(scen, cfg) {
  n <- cfg$n_students
  W <- scen$n_weeks
  group_A <- sample(rep(c(TRUE, FALSE), length.out = n))  # balanced randomisation

  b  <- rnorm(n, 0, cfg$sigma_student)
  s  <- rnorm(n, 0, cfg$sigma_slope)
  u  <- rnorm(W, 0, cfg$sigma_lesson)
  dw <- rnorm(W, 0, cfg$sigma_delta)

  student <- rep(seq_len(n), times = W)
  week    <- rep(seq_len(W), each = n)
  even    <- week %% 2L == 0L
  # A has access on even weeks, B on odd weeks
  access  <- as.integer(group_A[student] == even)

  eta   <- scen$alpha + b[student] + u[week] + rnorm(n * W, 0, cfg$sigma_sw) +
           (scen$delta + s[student] + dw[week]) * access
  score <- rbinom(n * W, scen$k_items, plogis(eta)) / scen$k_items

  keep <- runif(n * W) < scen$completion
  if (cfg$dropout > 0) {
    drops    <- runif(n) < cfg$dropout
    drop_wk  <- sample(seq_len(W), n, replace = TRUE)
    keep     <- keep & !(drops[student] & week > drop_wk[student])
  }
  data.frame(student = student[keep], lesson = week[keep],
             access = access[keep], score = score[keep])
}

run_one <- function(scen, cfg, seed) {
  set.seed(seed)
  dat  <- simulate_trial(scen, cfg)
  ctrl <- dat$score[dat$access == 0L]
  c(fit_primary(dat),
    n_obs        = nrow(dat),
    ceiling_ctrl = mean(ctrl == 1),   # share of perfect weekly scores, control arm
    sd_ctrl      = sd(ctrl))
}

## ---------------------------------------------------------------------------
## Run
## ---------------------------------------------------------------------------

n_scen <- nrow(grid)
message(sprintf("Power simulation: %d scenarios x %d sims = %s trials on %d cores%s",
                n_scen, cfg$n_sims, format(n_scen * cfg$n_sims, big.mark = ","),
                cfg$cores, if (quick) "  [QUICK MODE]" else ""))
message(sprintf("Output directory: %s", normalizePath(out_dir)))

set.seed(cfg$master_seed)
seed_mat <- matrix(sample.int(.Machine$integer.max, n_scen * cfg$n_sims),
                   nrow = n_scen)

raw_path <- file.path(out_dir, "power_raw.csv")
if (file.exists(raw_path)) file.remove(raw_path)

scen_cols <- c("scenario", "n_weeks", "completion", "k_items", "base_acc", "gain")
results <- vector("list", n_scen)
t0 <- Sys.time()

for (s in seq_len(n_scen)) {
  scen <- grid[s, ]
  sims <- mclapply(seq_len(cfg$n_sims),
                   function(i) run_one(scen, cfg, seed_mat[s, i]),
                   mc.cores = cfg$cores)
  m  <- as.data.frame(do.call(rbind, sims))
  df <- cbind(scen[rep(1L, nrow(m)), scen_cols, drop = FALSE],
              sim = seq_len(nrow(m)), m, row.names = NULL)
  results[[s]] <- df
  suppressWarnings(
    write.table(df, raw_path, sep = ",", row.names = FALSE,
                col.names = s == 1L, append = s > 1L)
  )
  el  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  eta <- el / s * (n_scen - s)
  message(sprintf("  scenario %3d/%d done  (elapsed %5.1f min, ~%5.1f min remaining)",
                  s, n_scen, el / 60, eta / 60))
}

raw <- bind_rows(results)

## ---------------------------------------------------------------------------
## Summaries
## ---------------------------------------------------------------------------

summ <- raw %>%
  group_by(scenario, n_weeks, completion, k_items, base_acc, gain) %>%
  summarise(
    n_sims        = n(),
    n_fit_failed  = sum(is.na(p)),
    power         = mean(p < cfg$alpha_level, na.rm = TRUE),
    power_mcse    = sqrt(power * (1 - power) / sum(!is.na(p))),
    mean_est_pp   = 100 * mean(est, na.rm = TRUE),   # bias check: should ~= 100*gain
    sd_est_pp     = 100 * sd(est, na.rm = TRUE),
    singular_full = mean(singular_full, na.rm = TRUE), # model 1 singular
    used_full     = mean(model_used == 1, na.rm = TRUE),
    used_uncorr   = mean(model_used == 2, na.rm = TRUE),
    used_intonly  = mean(model_used == 3, na.rm = TRUE),
    mean_sd_slope = mean(sd_slope, na.rm = TRUE),      # recovered student-slope SD
    ceiling_ctrl  = mean(ceiling_ctrl, na.rm = TRUE),
    mean_sd_ctrl  = mean(sd_ctrl, na.rm = TRUE),
    mean_n_obs    = mean(n_obs),
    .groups = "drop"
  ) %>%
  mutate(d_equiv = ifelse(gain > 0, gain / mean_sd_ctrl, NA_real_),
         extra_correct_per_quiz = gain * k_items)  # what the gain means in answers

write.csv(summ, file.path(out_dir, "power_summary.csv"), row.names = FALSE)

# Smallest swept effect reaching 80% power, per design cell
mde <- summ %>%
  filter(gain > 0) %>%
  group_by(n_weeks, completion, k_items, base_acc) %>%
  summarise(
    mde80_pp = if (any(power >= 0.80)) 100 * min(gain[power >= 0.80]) else NA_real_,
    .groups = "drop"
  )
write.csv(mde, file.path(out_dir, "power_mde.csv"), row.names = FALSE)

message("\n--- Type-I error check (gain = 0; should be ~= 0.05) ---")
print(as.data.frame(summ %>% filter(gain == 0) %>%
        select(n_weeks, completion, k_items, base_acc, power, power_mcse)),
      row.names = FALSE, digits = 3)

message("\n--- Smallest swept effect (accuracy pp) reaching 80% power ---")
print(as.data.frame(mde), row.names = FALSE, digits = 3)

message("\n--- Fallback ladder usage (share of trials analysed with each model) ---")
print(as.data.frame(summ %>%
        select(n_weeks, completion, k_items, gain, used_full, used_uncorr, used_intonly)),
      row.names = FALSE, digits = 3)

message("\n--- Ceiling diagnostics (share of perfect control-arm quiz scores) ---")
print(as.data.frame(summ %>% filter(gain == 0) %>%
        select(n_weeks, completion, k_items, base_acc, ceiling_ctrl)),
      row.names = FALSE, digits = 3)

## ---------------------------------------------------------------------------
## Power curves
## ---------------------------------------------------------------------------

plot_dat <- summ %>%
  mutate(
    completion_lab = factor(sprintf("%d%%", round(100 * completion)),
                            levels = sprintf("%d%%", sort(round(100 * unique(completion))))),
    items_lab    = factor(paste(k_items, "items"),
                          levels = paste(sort(unique(k_items)), "items")),
    weeks_lab    = factor(paste(n_weeks, "weeks"),
                          levels = paste(sort(unique(n_weeks)), "weeks")),
    base_acc_lab = factor(sprintf("Baseline accuracy %d%%", round(100 * base_acc)),
                          levels = sprintf("Baseline accuracy %d%%",
                                           sort(round(100 * unique(base_acc)))))
  )

# Sequential single-hue ramp (completion is ordered); CVD-checked.
ramp <- c("#6BAED6", "#2E7EBB", "#084594")
completion_cols <- setNames(ramp[seq_along(levels(plot_dat$completion_lab))],
                            levels(plot_dat$completion_lab))

p <- ggplot(plot_dat,
            aes(x = 100 * gain, y = power,
                colour = completion_lab, linetype = items_lab)) +
  geom_hline(yintercept = 0.80, colour = "grey70", linewidth = 0.4) +
  geom_hline(yintercept = cfg$alpha_level, colour = "grey85", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  facet_grid(base_acc_lab ~ weeks_lab) +
  scale_colour_manual(values = completion_cols, name = "Weekly completion") +
  scale_linetype_discrete(name = "Quiz length") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title    = "PREPARE intervention study: simulated power",
    subtitle = sprintf(paste0("n = %d students, week-wise alternating crossover, %d sims/cell\n",
                              "score ~ access + topic_week + (1 + access | student), ",
                              "with preregistered fallback ladder"),
                       cfg$n_students, cfg$n_sims),
    x = "True effect of access (accuracy percentage points)",
    y = sprintf("Power (two-sided alpha = %.2f)", cfg$alpha_level)
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey92"),
    legend.position  = "bottom",
    plot.title.position = "plot"
  )

ggsave(file.path(out_dir, "power_curves.png"), p,
       width = 9, height = 4.5, dpi = 300, bg = "white")

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

message(sprintf("\nDone in %.1f min. Results in %s",
                as.numeric(difftime(Sys.time(), t0, units = "mins")),
                normalizePath(out_dir)))
