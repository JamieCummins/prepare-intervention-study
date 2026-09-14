# Primary analysis model for the PREPARE intervention study, with the
# pre-specified fallback ladder (shared by power_simulation.R and
# power_sensitivity.R; sourced, not run directly).
#
# Preregistered model (RQ1, intention-to-treat). Topic-week enters as a
# fixed effect with sum-to-zero contrasts, so the access coefficient is
# the average effect across the topic-weeks in the data:
#
#   1. score ~ 1 + access + topic_week + (1 + access |  student)
#
# If model 1 is singular or does not converge, the intercept-slope
# correlation is dropped:
#
#   2. score ~ 1 + access + topic_week + (1 + access || student)
#
# If model 2 is singular or does not converge, the random slope is dropped:
#
#   3. score ~ 1 + access + topic_week + (1 | student)
#
# The access coefficient, its Satterthwaite p-value and the level of the
# ladder that was used are returned. In the simulation scripts the
# topic-week variable is called `lesson` (`lesson_f` once made a factor).

primary_formulas <- list(
  score ~ access + lesson_f + (1 + access |  student),
  score ~ access + lesson_f + (1 + access || student),
  score ~ access + lesson_f + (1 | student)
)

fit_is_clean <- function(fit) {
  !lme4::isSingular(fit) &&
    length(fit@optinfo$conv$lme4$messages) == 0L &&
    isTRUE(fit@optinfo$conv$opt == 0L)
}

slope_sd <- function(fit) {
  vc <- as.data.frame(lme4::VarCorr(fit))
  i  <- which(vc$grp == "student" & vc$var1 == "access" & is.na(vc$var2))
  if (length(i)) vc$sdcor[i[1]] else NA_real_
}

fit_primary <- function(dat) {
  dat$lesson_f <- factor(dat$lesson)
  contrasts(dat$lesson_f) <- contr.sum(nlevels(dat$lesson_f))
  out <- c(est = NA_real_, se = NA_real_, p = NA_real_,
           model_used = NA_real_, singular_full = NA_real_, sd_slope = NA_real_)
  singular_full <- NA_real_
  for (m in seq_along(primary_formulas)) {
    fit <- try(suppressMessages(suppressWarnings(
      lmerTest::lmer(primary_formulas[[m]], data = dat)
    )), silent = TRUE)
    if (inherits(fit, "try-error")) next
    if (m == 1L) singular_full <- as.numeric(lme4::isSingular(fit))
    if (!fit_is_clean(fit) && m < length(primary_formulas)) next
    cc <- tryCatch(coef(summary(fit)), error = function(e) NULL)
    if (is.null(cc) || !"access" %in% rownames(cc)) next
    out[] <- c(unname(cc["access", "Estimate"]),
               unname(cc["access", "Std. Error"]),
               unname(cc["access", "Pr(>|t|)"]),
               m, singular_full, slope_sd(fit))
    return(out)
  }
  out["singular_full"] <- singular_full
  out
}

# Alternative ladder used only in sensitivity check D (power_sensitivity.R):
# a random slope for access by topic-week in addition to the student slope.
# Correlations are dropped if a fit is singular or does not converge; the
# topic-week slope itself is never dropped, so a zero slope variance is
# accepted as a (singular) fit rather than triggering a further fallback.
topic_slope_formulas <- list(
  score ~ access + (1 + access |  student) + (1 + access |  lesson),
  score ~ access + (1 + access || student) + (1 + access || lesson),
  score ~ access + (1 | student) + (1 + access || lesson)
)

fit_topic_slope <- function(dat) {
  out <- c(est = NA_real_, se = NA_real_, p = NA_real_,
           model_used = NA_real_, singular_full = NA_real_, sd_slope = NA_real_)
  singular_full <- NA_real_
  for (m in seq_along(topic_slope_formulas)) {
    fit <- try(suppressMessages(suppressWarnings(
      lmerTest::lmer(topic_slope_formulas[[m]], data = dat)
    )), silent = TRUE)
    if (inherits(fit, "try-error")) next
    if (m == 1L) singular_full <- as.numeric(lme4::isSingular(fit))
    if (!fit_is_clean(fit) && m < length(topic_slope_formulas)) next
    cc <- tryCatch(coef(summary(fit)), error = function(e) NULL)
    if (is.null(cc) || !"access" %in% rownames(cc)) next
    vc <- as.data.frame(lme4::VarCorr(fit))
    i  <- which(vc$grp == "lesson" & vc$var1 == "access" & is.na(vc$var2))
    out[] <- c(unname(cc["access", "Estimate"]),
               unname(cc["access", "Std. Error"]),
               unname(cc["access", "Pr(>|t|)"]),
               m, singular_full, if (length(i)) vc$sdcor[i[1]] else NA_real_)
    return(out)
  }
  out["singular_full"] <- singular_full
  out
}

# Alternative used in sensitivity check E (power_sensitivity.R): topic-week
# as a random intercept instead of a fixed effect (the model considered
# before the fixed-effect specification was adopted), with the same student
# ladder, to confirm that the two specifications behave the same.
random_topic_formulas <- list(
  score ~ access + (1 + access |  student) + (1 | lesson),
  score ~ access + (1 + access || student) + (1 | lesson),
  score ~ access + (1 | student) + (1 | lesson)
)

fit_random_topic <- function(dat) {
  out <- c(est = NA_real_, se = NA_real_, p = NA_real_,
           model_used = NA_real_, singular_full = NA_real_, sd_slope = NA_real_)
  singular_full <- NA_real_
  for (m in seq_along(random_topic_formulas)) {
    fit <- try(suppressMessages(suppressWarnings(
      lmerTest::lmer(random_topic_formulas[[m]], data = dat)
    )), silent = TRUE)
    if (inherits(fit, "try-error")) next
    if (m == 1L) singular_full <- as.numeric(lme4::isSingular(fit))
    if (!fit_is_clean(fit) && m < length(random_topic_formulas)) next
    cc <- tryCatch(coef(summary(fit)), error = function(e) NULL)
    if (is.null(cc) || !"access" %in% rownames(cc)) next
    out[] <- c(unname(cc["access", "Estimate"]),
               unname(cc["access", "Std. Error"]),
               unname(cc["access", "Pr(>|t|)"]),
               m, singular_full, slope_sd(fit))
    return(out)
  }
  out["singular_full"] <- singular_full
  out
}
