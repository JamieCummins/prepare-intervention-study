#!/usr/bin/env Rscript
#
# Historical performance on the weekly Statistik I exercises (HS2025)
# ===================================================================
#
# Purpose: per-lesson average performance from last year's Qualtrics exports
# (data/old-reference-data/), to check for ceiling effects and to calibrate
# the variance components of preregistration/power_simulation.R.
#
# Data notes (established 2026-08-05):
#   * One CSV per weekly exercise (Uebung 1-13), standard Qualtrics export
#     with 3 header rows. `Pseudo` = student pseudonym, `SC0` = built-in
#     Qualtrics total score, question columns match ^U<w>F<k>.
#   * Exercise 13 is a comprehensive term-review test (all topics), reported
#     separately, not as a lesson.
#   * 25-35% of pseudonyms have multiple finished attempts (retakes, often
#     January = exam revision), and 80-100% of first attempts happen after
#     Dec 1. Performance therefore uses only students who completed the
#     exercise IN THE WEEK OF THE LESSON: first finished attempt within
#     7 days of the survey release. Releases were consecutive Wednesdays
#     (matching the practicals; verified from response timing, Sep 24 -
#     Dec 10 2025, exactly 7 days apart), so the window is Wed-Tue.
#     All-first-attempts-any-date is kept as a comparison column.
#   * The exports contain no answer key. Each week's maximum score is taken
#     as the OBSERVED maximum over all finished attempts (retakes included,
#     which helps: January revision retakes give many shots at a perfect
#     score). Where >= 5 students attain it exactly, it is treated as the
#     true maximum; otherwise (weeks 3, 8, 13) it is a lower bound, making
#     the reported mean percentage an UPPER bound for those weeks - which
#     only strengthens the no-ceiling conclusion. A least-squares
#     decomposition of SC0 over per-question answer patterns cannot pin the
#     maximum down further (students with a unique pattern on every question
#     make per-pattern values non-identified), but its ~1e-14 residual is
#     kept as a check that SC0 is a deterministic additive function of the
#     exercise answers alone (i.e. the self-report items are unscored).
#
# Usage:  Rscript analysis/historical_performance.R
# Outputs in analysis/historical_performance/:
#   lesson_performance.csv   per-lesson summary (the main deliverable)
#   max_diagnostics.csv      reconstructed vs observed maxima + residuals
#   first_attempts_long.csv  pseudonym x lesson first-attempt scores (reuse)
#   lesson_performance.png   figure
#   calibration_notes.txt    variance components for the power simulation

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(MASS, exclude = "select")
  library(lme4)
})

script_dir <- (function() {
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1]))) else getwd()
})()
data_dir <- normalizePath(file.path(script_dir, "..", "data", "old-reference-data"))
out_dir  <- file.path(script_dir, "historical_performance")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- sort(list.files(data_dir, pattern = "^[0-9]{2} Statistik I.*\\.csv$",
                         full.names = TRUE))
stopifnot(length(files) == 13)

IN_WEEK_DAYS <- 7   # lesson week: survey release (Wed) through the next Tuesday
SNAP_TOL <- 0.02   # snap reconstructed maxima to the nearest 0.25 within this

snap_quarter <- function(x) {
  s <- round(x * 4) / 4
  ifelse(abs(x - s) < SNAP_TOL, s, x)
}

## ---------------------------------------------------------------------------
## Read one export; reconstruct the week maximum; first-attempt scores
## ---------------------------------------------------------------------------

read_week <- function(fn) {
  wk <- as.integer(substr(basename(fn), 1, 2))
  # single parse: rows 2-3 (question text / ImportId) contain embedded
  # newlines, so header rows must be sliced as CSV records, never via skip=
  raw <- read.csv(fn, header = FALSE, colClasses = "character",
                  check.names = FALSE, encoding = "UTF-8")
  hdr <- as.character(raw[1, ])
  txt <- as.character(raw[2, ])
  df  <- raw[-(1:3), , drop = FALSE]
  names(df) <- hdr

  # topic label from the "verstanden" self-report question text, if present
  vi <- which(hdr == "verstanden_1")
  topic <- if (length(vi) == 1) {
    sub("^.*\\? - ", "", txt[vi])
  } else NA_character_

  fin <- df[df$Finished == "True" & nzchar(trimws(df$Pseudo)), ]
  fin$score  <- as.numeric(fin$SC0)
  fin$pseudo <- tolower(trimws(fin$Pseudo))
  fin$start  <- as.POSIXct(fin$StartDate, tz = "UTC")
  stopifnot(!anyNA(fin$score), !anyNA(fin$start))

  ## --- week maximum: observed max over ALL finished attempts --------------
  obsmax   <- snap_quarter(max(fin$score))
  n_at_obs <- sum(abs(fin$score - max(fin$score)) < 1e-6)
  max_pts  <- obsmax
  max_exact <- n_at_obs >= 5   # else: lower bound -> mean pct is an upper bound

  ## --- additivity check: SC0 is a function of exercise answers only -------
  qcols <- hdr[grepl(sprintf("^U%dF[0-9]+", wk), hdr) &
               !grepl("time", hdr, ignore.case = TRUE) &
               !grepl("Click", hdr)]
  qbase <- sub("^(U[0-9]+F[0-9]+).*", "\\1", qcols)
  Xs <- lapply(split(qcols, qbase), function(cs) {
    pat <- do.call(paste, c(fin[cs], list(sep = "\r")))
    u <- unique(pat)
    m <- outer(pat, u, `==`) * 1
    m
  })
  X  <- do.call(cbind, Xs)
  co <- as.vector(MASS::ginv(X) %*% fin$score)
  max_resid <- max(abs(as.vector(X %*% co) - fin$score))

  ## --- first finished attempt per pseudonym -------------------------------
  first <- fin %>%
    arrange(start) %>%
    distinct(pseudo, .keep_all = TRUE) %>%
    mutate(
      lesson  = wk,
      prop    = pmin(score / max_pts, 1),
      in_week = as.numeric(difftime(start, min(fin$start), units = "days"))
                <= IN_WEEK_DAYS
    ) %>%
    dplyr::select(lesson, pseudo, start, score, prop, in_week)

  list(
    first = first,
    diag  = data.frame(
      lesson = wk,
      topic = gsub("\\s*[\r\n]+\\s*", " - ", topic),
      n_questions = length(unique(qbase)),
      release = format(min(fin$start), "%Y-%m-%d"),
      max_used = max_pts, n_at_max = n_at_obs, max_exact = max_exact,
      additivity_resid = signif(max_resid, 3),
      n_finished = nrow(fin), n_first = nrow(first),
      n_in_week = sum(first$in_week)
    )
  )
}

weeks <- lapply(files, read_week)
long  <- bind_rows(lapply(weeks, `[[`, "first"))
diagn <- bind_rows(lapply(weeks, `[[`, "diag"))

write.csv(diagn, file.path(out_dir, "max_diagnostics.csv"), row.names = FALSE)
write.csv(long,  file.path(out_dir, "first_attempts_long.csv"), row.names = FALSE)

message("--- Week maxima and additivity check ---")
print(diagn %>% dplyr::select(-topic), row.names = FALSE)

## ---------------------------------------------------------------------------
## Per-lesson performance (lessons 1-12; exercise 13 reported separately)
## ---------------------------------------------------------------------------

summarise_perf <- function(d) {
  d %>% summarise(
    n            = n(),
    mean_pct     = 100 * mean(prop),
    sd_pct       = 100 * sd(prop),
    median_pct   = 100 * median(prop),
    q25_pct      = 100 * quantile(prop, .25),
    q75_pct      = 100 * quantile(prop, .75),
    pct_perfect  = 100 * mean(prop >= 1 - 1e-6),
    pct_ge90     = 100 * mean(prop >= 0.9),
    .groups = "drop"
  )
}

perf <- long %>% filter(in_week) %>% group_by(lesson) %>% summarise_perf() %>%
  left_join(diagn %>% dplyr::select(lesson, topic, max_used, max_exact),
            by = "lesson")

perf_any <- long %>% group_by(lesson) %>% summarise_perf() %>%
  dplyr::select(lesson, n_anytime = n, mean_pct_anytime = mean_pct)

perf <- left_join(perf, perf_any, by = "lesson") %>%
  relocate(lesson, topic)

write.csv(perf, file.path(out_dir, "lesson_performance.csv"), row.names = FALSE)

message("\n--- Per-lesson performance, completed in the lesson week (% of max) ---")
print(perf %>% dplyr::select(lesson, n, mean_pct, sd_pct, median_pct,
                             pct_perfect, pct_ge90, n_anytime, mean_pct_anytime) %>%
        mutate(across(where(is.numeric), ~ round(.x, 1))),
      row.names = FALSE)

## ---------------------------------------------------------------------------
## Figure (lessons 1-12)
## ---------------------------------------------------------------------------

p12 <- perf %>% filter(lesson <= 12) %>%
  mutate(se = sd_pct / sqrt(n),
         lo = mean_pct - 1.96 * se, hi = mean_pct + 1.96 * se)

plot_dat <- bind_rows(
  p12 %>% transmute(lesson, value = mean_pct, lo, hi,
                    metric = "Mean score (% of max)"),
  p12 %>% transmute(lesson, value = pct_perfect, lo = NA, hi = NA,
                    metric = "Students with perfect score (%)")
)

cols <- c("Mean score (% of max)"            = "#084594",
          "Students with perfect score (%)"  = "#6BAED6")

fig <- ggplot(plot_dat, aes(lesson, value, colour = metric)) +
  # data= subset: ggplot2 4.0.3 drops the whole layer if any colour-group
  # is all-NA on ymin/ymax, even with na.rm = TRUE
  geom_ribbon(aes(ymin = lo, ymax = hi), data = subset(plot_dat, !is.na(lo)),
              fill = "#084594", alpha = .18, colour = NA) +
  geom_line(linewidth = .7) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_x_continuous(breaks = 1:12) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(
    title    = "Statistik I HS2025: performance on weekly exercises",
    subtitle = sprintf(paste0("Students who completed the exercise in the ",
                              "lesson week (release Wed + 7 days),\nscored ",
                              "against the observed week maximum. ",
                              "n = %d-%d per lesson; band = 95%% CI."),
                       min(p12$n), max(p12$n)),
    x = "Lesson / exercise", y = "Percent"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey92"),
        legend.position = "bottom", plot.title.position = "plot")

ggsave(file.path(out_dir, "lesson_performance.png"), fig,
       width = 8.5, height = 5, dpi = 300, bg = "white")

## ---------------------------------------------------------------------------
## Calibration for preregistration/power_simulation.R (lessons 1-12)
## ---------------------------------------------------------------------------

cal_dat <- long %>% filter(lesson <= 12, in_week) %>%
  mutate(pseudo = factor(pseudo), lesson_f = factor(lesson))
fit <- lmer(prop ~ 1 + (1 | pseudo) + (1 | lesson_f), data = cal_dat)
vc  <- as.data.frame(VarCorr(fit))
sd_student <- vc$sdcor[vc$grp == "pseudo"]
sd_lesson  <- vc$sdcor[vc$grp == "lesson_f"]
sd_resid   <- vc$sdcor[vc$grp == "Residual"]
pbar       <- mean(cal_dat$prop)
to_logit   <- 1 / (pbar * (1 - pbar))   # delta-method slope at the mean

cohort <- long %>% filter(lesson <= 12) %>% distinct(pseudo) %>% nrow()
attend <- long %>% filter(lesson <= 12, in_week) %>% count(lesson) %>%
  mutate(rate = n / cohort)

cal <- c(
  sprintf(paste0("Calibration from HS2025 in-week completions, lessons 1-12 ",
                 "(n obs = %d from %d students; cohort of ever-active pseudonyms = %d)"),
          nrow(cal_dat), nlevels(droplevels(cal_dat$pseudo)), cohort),
  "  CAUTION: small, self-selected in-week sample; treat as rough anchors.",
  sprintf("  Marginal mean score           : %.3f  (power-sim base_acc)", pbar),
  "  Gaussian LMM on proportion scale: prop ~ 1 + (1|student) + (1|lesson)",
  sprintf("    SD student   = %.3f   -> logit approx %.2f  (power-sim sigma_student)", sd_student, sd_student * to_logit),
  sprintf("    SD lesson    = %.3f   -> logit approx %.2f  (power-sim sigma_lesson)",  sd_lesson,  sd_lesson  * to_logit),
  sprintf("    SD residual  = %.3f   (includes item-sampling noise; power-sim sigma_sw", sd_resid),
  "                            should be set BELOW its logit equivalent because the",
  "                            binomial item draw already contributes this noise)",
  sprintf("  Perfect-score share overall   : %.1f%%", 100 * mean(cal_dat$prop >= 1 - 1e-6)),
  sprintf("  In-week completion, median    : %d students/lesson (%.0f%% of cohort)",
          round(median(attend$n)), 100 * median(attend$rate)),
  "  NOTE: voluntary in-week uptake was this low last year; the study's",
  "  completion assumptions must come from the incentivised assessment design,",
  "  not from these historical rates."
)
writeLines(cal, file.path(out_dir, "calibration_notes.txt"))
message("\n", paste(cal, collapse = "\n"))
