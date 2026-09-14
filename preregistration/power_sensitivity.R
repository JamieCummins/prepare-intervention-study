#!/usr/bin/env Rscript
#
# Sensitivity checks for the power-simulation assumptions that the HS2025
# historical data cannot verify (companion to power_simulation.R).
#
#   A. Topic-varying treatment effects. The main simulation assumes the
#      access effect is constant across lessons (sigma_delta = 0). Here
#      sigma_delta is swept over {0, 0.1, 0.2} on the log-odds scale
#      (approx. 0, 2.4, 4.9 pp SD of per-topic effects) to price the
#      power cost of heterogeneity. gain = 0 rows show that a model
#      without an access-by-topic term (as preregistered) rejects more
#      often than alpha when the effect truly varies across topics: the
#      preregistered test concerns the average effect across the topics
#      studied, not a population of topics.
#
#   B. Non-random assessment completion. The main simulation assumes
#      MCAR completion, independent of the access arm. Here completion
#      follows P(complete) = plogis(c0 + 0.5*b_i + dcomp*access):
#      diligent/able students complete more (0.5*b_i, as the HS2025
#      self-selection suggests), and access itself may boost completion
#      (dcomp in {0, 0.3, 0.6} logits ~ 0, +7, +13 pp). Ability selection
#      alone cancels between arms; an access-dependent boost pulls
#      marginal (weaker) students into the access arm's observed sample,
#      biasing the ITT estimate downward. gain = 0 rows show whether
#      that distortion also corrupts the type-I error rate.
#
#   D. Topic-week random slope. Check A shows what happens when the
#      preregistered model omits a topic-week slope while the true effect
#      varies across topics. Here the same data are analysed with a model
#      that adds (1 + access | lesson) (fit_topic_slope in power_fit.R) to
#      quantify the type-I error protection and the power cost of the
#      extra variance component with only 10 topic-weeks.
#
#   E. Topic-week as a random intercept. The same grid as check A analysed
#      with topic-week entered as a random intercept rather than a fixed
#      effect (fit_random_topic in power_fit.R), to confirm that the two
#      specifications give the same estimate, power and
#      (conditional-on-topics) type-I behaviour.
#
#   C. Student-varying treatment effects. The primary model includes a
#      random slope for access by student, but sigma_slope has no data
#      anchor. Here it is swept over {0, 0.2, 0.4} (approx. 0, 5, 10 pp
#      SD of per-student effects); the run reports power, bias and how
#      often the preregistered fallback ladder is triggered.
#
# Fixed design cell: n = 400, 10 weeks, k = 8 items, base_acc = .58,
# marginal completion ~= .60, calibrated variance components, and the
# preregistered model with fallback ladder (power_fit.R).
#
# Usage: Rscript preregistration/power_sensitivity.R  [--sims N] [--cores N] [--checks A,B,C,D,E]

suppressPackageStartupMessages({
  library(lme4); library(lmerTest); library(parallel)
})

script_dir <- (function() {
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1]))) else getwd()
})()
source(file.path(script_dir, "power_fit.R"))

args <- commandArgs(trailingOnly = TRUE)
val_after <- function(flag, cast) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) cast(args[i + 1]) else NULL
}
N_SIMS <- val_after("--sims", as.integer);  if (is.null(N_SIMS)) N_SIMS <- 400
CORES  <- val_after("--cores", as.integer); if (is.null(CORES))  CORES  <- max(1L, detectCores() - 1L)
CHECKS <- val_after("--checks", function(x) strsplit(x, ",")[[1]]); if (is.null(CHECKS)) CHECKS <- c("A", "B", "C", "D", "E")

cfg <- list(n_students = 400, n_weeks = 10, k_items = 8, base_acc = 0.58,
            completion = 0.60, sigma_student = 0.6, sigma_lesson = 0.5,
            sigma_sw = 0.2, sigma_slope = 0.2, alpha_level = 0.05,
            master_seed = 20260905)

marginal_p <- function(a, sigma) {
  integrate(function(z) plogis(a + sigma * z) * dnorm(z), -Inf, Inf)$value
}
calibrate <- function(target, sigma) {
  uniroot(function(a) marginal_p(a, sigma) - target, c(-15, 15))$root
}

sigma_ctrl <- with(cfg, sqrt(sigma_student^2 + sigma_lesson^2 + sigma_sw^2))
alpha0     <- calibrate(cfg$base_acc, sigma_ctrl)
delta_for  <- function(gain, sigma_slope, sigma_delta) {
  calibrate(cfg$base_acc + gain,
            sqrt(sigma_ctrl^2 + sigma_slope^2 + sigma_delta^2)) - alpha0
}
# completion intercept: ability selection strength 0.5 -> latent SD 0.5*sigma_student
c0 <- calibrate(cfg$completion, 0.5 * cfg$sigma_student)

run_one <- function(gain, sigma_delta, dcomp, sigma_slope, seed, fitter = fit_primary) {
  set.seed(seed)
  n <- cfg$n_students; W <- cfg$n_weeks
  group_A <- sample(rep(c(TRUE, FALSE), length.out = n))
  b  <- rnorm(n, 0, cfg$sigma_student)
  s  <- rnorm(n, 0, sigma_slope)
  u  <- rnorm(W, 0, cfg$sigma_lesson)
  dw <- rnorm(W, 0, sigma_delta)
  student <- rep(seq_len(n), times = W); week <- rep(seq_len(W), each = n)
  access  <- as.integer(group_A[student] == (week %% 2L == 0L))
  eta   <- alpha0 + b[student] + u[week] + rnorm(n * W, 0, cfg$sigma_sw) +
           (delta_for(gain, sigma_slope, sigma_delta) + s[student] + dw[week]) * access
  score <- rbinom(n * W, cfg$k_items, plogis(eta)) / cfg$k_items
  keep  <- runif(n * W) < plogis(c0 + 0.5 * b[student] + dcomp * access)
  dat <- data.frame(student = student[keep], lesson = week[keep],
                    access = access[keep], score = score[keep])
  c(fitter(dat),
    comp_acc = mean(keep[access == 1L]), comp_ctl = mean(keep[access == 0L]))
}

grid <- rbind(
  expand.grid(check = "A_topic_heterogeneity", gain = c(0, 0.02, 0.04),
              sigma_delta = c(0, 0.1, 0.2), dcomp = 0, sigma_slope = cfg$sigma_slope,
              stringsAsFactors = FALSE),
  expand.grid(check = "B_differential_completion", gain = c(0, 0.04),
              sigma_delta = 0, dcomp = c(0, 0.3, 0.6), sigma_slope = cfg$sigma_slope,
              stringsAsFactors = FALSE),
  expand.grid(check = "C_student_heterogeneity", gain = c(0, 0.02),
              sigma_delta = 0, dcomp = 0, sigma_slope = c(0, 0.2, 0.4),
              stringsAsFactors = FALSE),
  expand.grid(check = "D_topic_slope_model", gain = c(0, 0.02, 0.04),
              sigma_delta = c(0, 0.1, 0.2), dcomp = 0, sigma_slope = cfg$sigma_slope,
              stringsAsFactors = FALSE),
  expand.grid(check = "E_random_topic_model", gain = c(0, 0.02, 0.04),
              sigma_delta = c(0, 0.1, 0.2), dcomp = 0, sigma_slope = cfg$sigma_slope,
              stringsAsFactors = FALSE)
)
grid <- grid[substr(grid$check, 1, 1) %in% CHECKS, ]
rownames(grid) <- NULL

set.seed(cfg$master_seed)
seeds <- matrix(sample.int(.Machine$integer.max, nrow(grid) * N_SIMS), nrow(grid))
t0 <- Sys.time()
res <- lapply(seq_len(nrow(grid)), function(g) {
  sims <- mclapply(seq_len(N_SIMS), function(i)
    run_one(grid$gain[g], grid$sigma_delta[g], grid$dcomp[g], grid$sigma_slope[g],
            seeds[g, i],
            fitter = switch(grid$check[g],
                            D_topic_slope_model = fit_topic_slope,
                            E_random_topic_model = fit_random_topic,
                            fit_primary)),
    mc.cores = CORES)
  m <- do.call(rbind, sims)
  ok <- !is.na(m[, "p"])
  message(sprintf("  cell %2d/%d done (%.1f min)", g, nrow(grid),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  data.frame(grid[g, ],
    n_ok        = sum(ok),
    power       = mean(m[ok, "p"] < cfg$alpha_level),
    mean_est_pp = 100 * mean(m[ok, "est"]),
    bias_pp     = 100 * (mean(m[ok, "est"]) - grid$gain[g]),
    used_full     = mean(m[ok, "model_used"] == 1),
    used_uncorr   = mean(m[ok, "model_used"] == 2),
    used_intonly  = mean(m[ok, "model_used"] == 3),
    mean_sd_slope = mean(m[ok, "sd_slope"], na.rm = TRUE),
    completion_access  = mean(m[ok, "comp_acc"]),
    completion_control = mean(m[ok, "comp_ctl"]))
})
res <- do.call(rbind, res)

out <- file.path(script_dir, "power_analysis",
                 if (setequal(CHECKS, c("A", "B", "C", "D", "E"))) "sensitivity_checks.csv"
                 else sprintf("sensitivity_checks_%s.csv", paste(sort(CHECKS), collapse = "")))
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
write.csv(res, out, row.names = FALSE)
print(res, row.names = FALSE, digits = 3)
