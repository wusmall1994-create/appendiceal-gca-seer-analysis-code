# R verification of Task 2 (race/ethnicity-stratified survival) previously
# computed in Python by revision_analysis/run_phase1.py.
#
# Study: appendiceal goblet cell adenocarcinoma (GCA) case-definition research
# in SEER-17, 2000-2023 (revision phase).
#
# Replicates the exact Python conventions in run_phase1.py:
#   - cohort: ICD-O-3 Hist/behav starts with "8243" (D1), diagnosed 2000-2018
#   - survival: t = min(Survival months, 60); event = (Vital status == "Dead") & (surv_m <= 60)
#   - no special handling of Survival months == 0 (kept as-is)
#   - KM point estimate; CI: Python used log-log CI (Greenwood), z = 1.96 exactly.
#     R survfit log-log CI uses z = qnorm(0.975) = 1.959964 -> reported both ways.
#   - log-rank: standard k-group (survdiff)
#   - Cox: univariate, Breslow ties (Python used Breslow; R default is Efron -> set ties="breslow")
#
# Required private input (not distributed; see README):
#   raw/gca_case_listing_seer17_2000_2023.csv
#
# Output: task2_r_km.csv + task2_r_tests.csv under outputs/revision_analysis/
# (excluded from version control), plus a console summary.
#
# Dependencies: R survival package.

library(survival)

SRC <- file.path("raw", "gca_case_listing_seer17_2000_2023.csv")
OUTDIR <- file.path("outputs", "revision_analysis")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("R version:", R.version.string, "\n")
cat("survival version:", as.character(packageVersion("survival")), "\n\n")

df <- read.csv(SRC, colClasses = "character", check.names = FALSE)
d1 <- df[startsWith(df[["ICD-O-3 Hist/behav, malignant"]], "8243"), ]
stopifnot(nrow(d1) == 2571)
cat("D1 cohort n =", nrow(d1), "(expect 2571)\n")

RACE <- "Race recode (with detailed Asian and Native Hawaiian other PI)"
API_LBL <- c("Chinese", "Filipino", "Hawaiian", "Asian Indian, Pakistani",
             "Other Asian American", "Korean", "Vietnamese", "Japanese",
             "Laotian", "Samoan", "Guamanian/Chamorro", "Other Pacific Islander")
race4 <- function(x) {
  ifelse(x == "White", "White",
  ifelse(x == "Black", "Black",
  ifelse(x %in% API_LBL, "Asian/Pacific Islander", "Other/Unknown")))
}
d1$race4 <- race4(d1[[RACE]])
d1$hisp  <- ifelse(startsWith(d1[["Origin recode NHIA (Hispanic, Non-Hisp)"]], "Spanish"),
                   "Hispanic", "Non-Hispanic")
d1$year  <- as.integer(d1[["Year of diagnosis"]])

d1$surv_m <- suppressWarnings(as.numeric(d1[["Survival months"]]))
d1$dead   <- as.integer(d1[["Vital status recode (study cutoff used)"]] == "Dead")

sv <- d1[d1$year <= 2018 & !is.na(d1$surv_m), ]
sv$t <- pmin(sv$surv_m, 60)
sv$e <- as.integer(sv$dead == 1 & sv$surv_m <= 60)
cat("2000-2018 survival subset n =", nrow(sv), "\n\n")

# ---------- KM at 60 months ----------
km_at60 <- function(sub) {
  fit <- survfit(Surv(t, e) ~ 1, data = sub, conf.type = "log-log")
  s <- summary(fit, times = 60, extend = TRUE)
  S <- s$surv; se <- s$std.err           # se = S * sqrt(greenwood var)
  # R native log-log CI (z = qnorm(0.975))
  lo_r <- s$lower; hi_r <- s$upper
  # Python-exact log-log CI (z = 1.96 literally)
  if (!is.na(S) && S > 0 && S < 1 && se > 0) {
    theta <- log(-log(S))
    se_th <- se / (S * abs(log(S)))
    lo_py <- exp(-exp(theta + 1.96 * se_th))
    hi_py <- exp(-exp(theta - 1.96 * se_th))
  } else { lo_py <- NA_real_; hi_py <- NA_real_ }
  # R default (conf.type="log") CI for reference
  fit_log <- survfit(Surv(t, e) ~ 1, data = sub, conf.type = "log")
  s_log <- summary(fit_log, times = 60, extend = TRUE)
  data.frame(n = nrow(sub), events = sum(sub$e),
             S60 = S,
             lo_loglog_R = lo_r, hi_loglog_R = hi_r,
             lo_loglog_py196 = lo_py, hi_loglog_py196 = hi_py,
             lo_log_R = s_log$lower, hi_log_R = s_log$upper)
}

rows <- list()
groups_race <- c("White", "Black", "Asian/Pacific Islander")
for (g in groups_race) {
  r <- km_at60(sv[sv$race4 == g, ])
  rows[[paste0("Race ", g)]] <- cbind(comparison = paste0("Race ", g, ": 5-year OS (KM)"), r)
}

# log-rank: White vs Black vs API
sv3 <- sv[sv$race4 %in% groups_race, ]
sv3$race4 <- factor(sv3$race4, levels = groups_race)
lr <- survdiff(Surv(t, e) ~ race4, data = sv3)
chi2_race <- lr$chisq; df_race <- length(lr$n) - 1
p_race <- pchisq(chi2_race, df_race, lower.tail = FALSE)

# ---------- Cox (Breslow ties, univariate) ----------
cox_pair <- function(dat, grp, ref, cmp) {
  dd <- dat[dat[[grp]] %in% c(ref, cmp), ]
  dd$g <- factor(dd[[grp]], levels = c(ref, cmp))
  fit <- coxph(Surv(t, e) ~ g, data = dd, ties = "breslow")
  co <- summary(fit)$coefficients
  ci <- summary(fit)$conf.int
  list(hr = unname(ci[1, "exp(coef)"]), lo = unname(ci[1, "lower .95"]),
       hi = unname(ci[1, "upper .95"]), p = unname(co[1, "Pr(>|z|)"]),
       n = nrow(dd), events = sum(dd$e))
}
cx_b  <- cox_pair(sv, "race4", "White", "Black")
cx_a  <- cox_pair(sv, "race4", "White", "Asian/Pacific Islander")

# ---------- Ethnicity ----------
for (g in c("Hispanic", "Non-Hispanic")) {
  r <- km_at60(sv[sv$hisp == g, ])
  rows[[g]] <- cbind(comparison = paste0(g, ": 5-year OS (KM)"), r)
}
lr2 <- survdiff(Surv(t, e) ~ hisp, data = sv)
chi2_h <- lr2$chisq; p_h <- pchisq(chi2_h, 1, lower.tail = FALSE)
cx_h <- cox_pair(sv, "hisp", "Non-Hispanic", "Hispanic")

# ---------- assemble ----------
km_df <- do.call(rbind, rows)
rownames(km_df) <- NULL

tests_df <- data.frame(
  comparison = c("Race: White vs Black vs API (log-rank)",
                 "Ethnicity: Hispanic vs Non-Hispanic (log-rank)",
                 "Cox HR Black vs White (unadjusted, Breslow)",
                 "Cox HR API vs White (unadjusted, Breslow)",
                 "Cox HR Hispanic vs Non-Hispanic (unadjusted, Breslow)"),
  stat = c(sprintf("chi2(%d)=%.3f", df_race, chi2_race),
           sprintf("chi2(1)=%.3f", chi2_h),
           sprintf("HR %.2f (95%% CI %.2f-%.2f)", cx_b$hr, cx_b$lo, cx_b$hi),
           sprintf("HR %.2f (95%% CI %.2f-%.2f)", cx_a$hr, cx_a$lo, cx_a$hi),
           sprintf("HR %.2f (95%% CI %.2f-%.2f)", cx_h$hr, cx_h$lo, cx_h$hi)),
  p = c(p_race, p_h, cx_b$p, cx_a$p, cx_h$p),
  n = c(NA, NA, cx_b$n, cx_a$n, cx_h$n),
  events = c(NA, NA, cx_b$events, cx_a$events, cx_h$events)
)

cat("===== KM 5-year OS (60 months) =====\n")
print(km_df, digits = 6)
cat("\n===== Tests =====\n")
print(tests_df, digits = 8)

# full-precision p values
cat("\nFull-precision p values:\n")
cat(sprintf("race log-rank p   = %.15g\n", p_race))
cat(sprintf("hisp log-rank p   = %.15g\n", p_h))
cat(sprintf("cox Black p       = %.15g\n", cx_b$p))
cat(sprintf("cox API p         = %.15g\n", cx_a$p))
cat(sprintf("cox Hispanic p    = %.15g\n", cx_h$p))
cat("\nFull-precision HRs:\n")
cat(sprintf("Black vs White     HR=%.10f CI %.10f-%.10f\n", cx_b$hr, cx_b$lo, cx_b$hi))
cat(sprintf("API vs White       HR=%.10f CI %.10f-%.10f\n", cx_a$hr, cx_a$lo, cx_a$hi))
cat(sprintf("Hisp vs Non-Hisp   HR=%.10f CI %.10f-%.10f\n", cx_h$hr, cx_h$lo, cx_h$hi))

write.csv(km_df, file.path(OUTDIR, "task2_r_km.csv"), row.names = FALSE)
write.csv(tests_df, file.path(OUTDIR, "task2_r_tests.csv"), row.names = FALSE)
cat("\nCSV written.\n")

# ---------- group sizes cross-check ----------
cat("\nGroup sizes / events (2000-2018 subset):\n")
print(aggregate(e ~ race4, sv, function(x) c(n = length(x), events = sum(x))))
print(aggregate(e ~ hisp, sv, function(x) c(n = length(x), events = sum(x))))
cat("Total D1 by race4:\n"); print(table(d1$race4))
