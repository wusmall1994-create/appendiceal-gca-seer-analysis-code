options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
  library(patchwork)
  library(svglite)
})

root <- getwd()
if (!dir.exists(file.path(root, "raw")) && !dir.exists(file.path(root, "raw data"))) {
  script_arg <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1])
  script_dir <- dirname(normalizePath(script_arg, winslash = "/", mustWork = FALSE))
  root <- if (basename(script_dir) == "R") dirname(script_dir) else script_dir
}
raw_dir <- file.path(root, "raw")
if (!dir.exists(raw_dir)) raw_dir <- file.path(root, "raw data")
out_dir <- file.path(root, "outputs", "formal_next_phase")
fig_dir <- Sys.getenv("GCA_FIG_DIR", unset = file.path(out_dir, "figures"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

inc_file <- file.path(raw_dir, "gca_incidence_age_hist_year_seer17_2000_2023.csv")
case_file <- file.path(raw_dir, "gca_case_listing_seer17_2000_2023.csv")
stopifnot(file.exists(inc_file), file.exists(case_file))

definitions <- list(
  D1_strict_GCA = 8243L,
  D2_legacy_GCC_adenocarcinoid = c(8243L, 8245L),
  D3_historical_expanded_spectrum = c(8243L, 8244L, 8245L)
)
periods <- list(`2000-2004` = 2000:2004, `2005-2009` = 2005:2009,
                `2010-2014` = 2010:2014, `2015-2019` = 2015:2019,
                `2020-2023` = 2020:2023, `2000-2023` = 2000:2023)
std_pop <- c(`00`=3794901, `01-04`=15191619, `05-09`=19919840,
             `10-14`=20056779, `15-19`=19819518, `20-24`=18257225,
             `25-29`=17722067, `30-34`=19511370, `35-39`=22179956,
             `40-44`=22479229, `45-49`=19805793, `50-54`=17224359,
             `55-59`=13307234, `60-64`=10654272, `65-69`=9409940,
             `70-74`=8725574, `75-79`=7414559, `80-84`=4900234,
             `85-89`=2678567, `90+`=1580606)
std_w <- std_pop / sum(std_pop)

write_csv <- function(x, name) write.csv(x, file.path(out_dir, name), row.names = FALSE, fileEncoding = "UTF-8")
extract_hist <- function(x) suppressWarnings(as.integer(substr(x, 1, 4)))
age_key <- function(x) {
  y <- trimws(tolower(gsub("years|year", "", x)))
  y <- gsub(" ", "", y)
  y[y %in% c("00", "0", "<1")] <- "00"
  y[grepl("^90", y)] <- "90+"
  y
}
age_mid <- function(x) {
  z <- age_key(x)
  ans <- rep(NA_real_, length(z))
  ans[z == "00"] <- 0.5
  ans[z == "90+"] <- 92
  ok <- grepl("^[0-9]{2}-[0-9]{2}$", z)
  ans[ok] <- rowMeans(do.call(rbind, strsplit(z[ok], "-", fixed = TRUE)) |> apply(2, as.numeric))
  ans
}
stage_group <- function(x) {
  z <- tolower(x)
  ans <- rep("Unknown", length(z))
  ans[grepl("localized", z)] <- "Localized"
  ans[grepl("regional", z)] <- "Regional"
  ans[grepl("distant", z)] <- "Distant"
  factor(ans, levels = c("Localized", "Regional", "Distant", "Unknown"))
}
period_group <- function(y) {
  cut(y, breaks = c(1999, 2004, 2009, 2014, 2019, 2023),
      labels = c("2000-2004", "2005-2009", "2010-2014", "2015-2019", "2020-2023"))
}

# Fay-Feuer gamma interval for a directly standardized rate.
direct_rate <- function(cases, pop, weights = std_w[names(cases)], alpha = .05) {
  ok <- is.finite(cases) & is.finite(pop) & pop > 0 & is.finite(weights)
  d <- cases[ok]; n <- pop[ok]; w <- weights[ok]
  r <- sum(w * d / n)
  v <- sum(w^2 * d / n^2)
  wm <- max(w / n)
  if (r <= 0 || v <= 0) return(c(rate=0, se=0, lcl=0, ucl=-log(alpha) * wm * 100000))
  lcl <- qgamma(alpha/2, shape = r^2/v, scale = v/r)
  ucl <- qgamma(1-alpha/2, shape = (r + wm)^2/(v + wm^2), scale = (v + wm^2)/(r + wm))
  c(rate=r*100000, se=sqrt(v)*100000, lcl=lcl*100000, ucl=ucl*100000)
}

inc <- read.csv(inc_file, check.names = FALSE, encoding = "UTF-8-BOM")
inc$histology <- extract_hist(inc[["ICD-O-3 Hist/behav, malignant"]])
inc$year <- suppressWarnings(as.integer(inc[["Year of diagnosis"]]))
inc$age <- age_key(inc[["Age recode with <1 year olds and 90+"]])
inc$cases <- suppressWarnings(as.numeric(inc[["Count"]]))
inc$population <- suppressWarnings(as.numeric(inc[["Population"]]))
inc <- inc[inc$histology %in% 8243:8245 & inc$year %in% 2000:2023 & inc$age %in% names(std_w), ]

make_cells <- function(years, codes) {
  z <- inc[inc$year %in% years & inc$histology %in% codes, ]
  counts <- aggregate(cases ~ year + age, z, sum)
  pops <- aggregate(population ~ year + age, z, max)
  a <- merge(counts, pops, by = c("year", "age"), all = TRUE)
  a$cases[is.na(a$cases)] <- 0
  counts2 <- aggregate(cases ~ age, a, sum)
  pops2 <- aggregate(population ~ age, a, sum)
  merge(counts2, pops2, by = "age", all = TRUE)
}

annual_rows <- list(); k <- 1
for (dn in names(definitions)) for (yr in 2000:2023) {
  cell <- make_cells(yr, definitions[[dn]])
  names_c <- cell$age; d <- setNames(cell$cases, names_c); n <- setNames(cell$population, names_c)
  est <- direct_rate(d, n, std_w[names_c])
  annual_rows[[k]] <- data.frame(definition=dn, year=yr, cases=sum(d),
    age_adjusted_rate_per_100k=est["rate"], se=est["se"], lcl_95=est["lcl"], ucl_95=est["ucl"])
  k <- k + 1
}
annual <- do.call(rbind, annual_rows)
write_csv(annual, "formal_incidence_annual_fay_feuer.csv")

period_rows <- list(); k <- 1; set.seed(20260818)
for (pn in names(periods)) {
  yrs <- periods[[pn]]
  code_cells <- lapply(c(8243L,8244L,8245L), function(h) make_cells(yrs, h))
  ages <- names(std_w)
  getv <- function(cell, nm) setNames(cell[[nm]], cell$age)[ages] |> unname()
  count_mat <- cbind(`8243`=getv(code_cells[[1]],"cases"), `8244`=getv(code_cells[[2]],"cases"), `8245`=getv(code_cells[[3]],"cases"))
  count_mat[is.na(count_mat)] <- 0
  pop <- getv(code_cells[[1]], "population"); w <- unname(std_w[ages]); denom <- pop
  B <- 10000L
  sim43 <- matrix(rpois(B*length(ages), rep(count_mat[,"8243"], each=B)), nrow=B)
  sim44 <- matrix(rpois(B*length(ages), rep(count_mat[,"8244"], each=B)), nrow=B)
  sim45 <- matrix(rpois(B*length(ages), rep(count_mat[,"8245"], each=B)), nrow=B)
  std <- function(m) as.vector(m %*% (w/denom) * 100000)
  sims <- list(D1_strict_GCA=std(sim43), D2_legacy_GCC_adenocarcinoid=std(sim43+sim45),
               D3_historical_expanded_spectrum=std(sim43+sim44+sim45))
  base_sim <- sims[[1]]
  for (dn in names(definitions)) {
    codes <- definitions[[dn]]; idx <- match(as.character(codes), colnames(count_mat))
    d <- rowSums(count_mat[,idx,drop=FALSE]); names(d) <- ages; names(pop) <- ages
    est <- direct_rate(d, pop, std_w)
    inflation <- 100*(est["rate"] / direct_rate(setNames(count_mat[,"8243"],ages), pop, std_w)["rate"] - 1)
    infl_sim <- 100*(sims[[dn]]/base_sim - 1); infl_sim <- infl_sim[is.finite(infl_sim)]
    diff_sim <- sims[[dn]] - base_sim
    period_rows[[k]] <- data.frame(period=pn, definition=dn, cases=sum(d),
      age_adjusted_rate_per_100k=est["rate"], lcl_95=est["lcl"], ucl_95=est["ucl"],
      rate_inflation_vs_D1_percent=inflation,
      inflation_lcl_95=if(dn=="D1_strict_GCA") 0 else quantile(infl_sim,.025,na.rm=TRUE),
      inflation_ucl_95=if(dn=="D1_strict_GCA") 0 else quantile(infl_sim,.975,na.rm=TRUE),
      rate_difference_vs_D1=est["rate"]-direct_rate(setNames(count_mat[,"8243"],ages),pop,std_w)["rate"],
      difference_lcl_95=if(dn=="D1_strict_GCA") 0 else quantile(diff_sim,.025,na.rm=TRUE),
      difference_ucl_95=if(dn=="D1_strict_GCA") 0 else quantile(diff_sim,.975,na.rm=TRUE))
    k <- k + 1
  }
}
period_rates <- do.call(rbind, period_rows)
write_csv(period_rates, "formal_incidence_period_definition_effect.csv")

# Weighted one-breakpoint log-linear trend for strict GCA.
trend_dat <- annual[annual$definition == "D1_strict_GCA" & annual$age_adjusted_rate_per_100k > 0,]
trend_dat$log_rate <- log(trend_dat$age_adjusted_rate_per_100k)
trend_dat$wt <- 1/(trend_dat$se/trend_dat$age_adjusted_rate_per_100k)^2
fits <- lapply(2005:2018, function(bp) {
  z <- trend_dat; z$hinge <- pmax(0,z$year-bp)
  fit <- lm(log_rate ~ year + hinge, data=z, weights=wt)
  list(bp=bp, fit=fit, aic=AIC(fit))
})
best <- fits[[which.min(sapply(fits, `[[`, "aic"))]]
b <- coef(best$fit); V <- vcov(best$fit)
s1 <- b["year"]; se1 <- sqrt(V["year","year"])
s2 <- b["year"]+b["hinge"]; se2 <- sqrt(V["year","year"]+V["hinge","hinge"]+2*V["year","hinge"])
apc <- function(s,se) c(APC=100*(exp(s)-1), LCL=100*(exp(s-1.96*se)-1), UCL=100*(exp(s+1.96*se)-1))
overall_fit <- lm(log_rate ~ year, data=trend_dat, weights=wt)
so <- coef(overall_fit)["year"]; seo <- sqrt(vcov(overall_fit)["year","year"])
ao <- apc(so,seo); a1 <- apc(s1,se1); a2 <- apc(s2,se2)
trend_out <- rbind(
  data.frame(model="Overall WLS", segment="2000-2023", breakpoint=NA, APC=unname(ao[1]), LCL_95=unname(ao[2]), UCL_95=unname(ao[3]), AIC=AIC(overall_fit)),
  data.frame(model="One-breakpoint WLS", segment=paste0("2000-",best$bp), breakpoint=best$bp, APC=unname(a1[1]), LCL_95=unname(a1[2]), UCL_95=unname(a1[3]), AIC=best$aic),
  data.frame(model="One-breakpoint WLS", segment=paste0(best$bp+1,"-2023"), breakpoint=best$bp, APC=unname(a2[1]), LCL_95=unname(a2[2]), UCL_95=unname(a2[3]), AIC=best$aic))
write_csv(trend_out, "strict_GCA_weighted_segmented_trend.csv")

cases <- read.csv(case_file, check.names = FALSE, encoding = "UTF-8-BOM")
cases$histology <- extract_hist(cases[["ICD-O-3 Hist/behav, malignant"]])
cases$year <- as.integer(cases[["Year of diagnosis"]])
cases$age_num <- age_mid(cases[["Age recode with <1 year olds and 90+"]])
cases$period <- period_group(cases$year)
cases$stage <- stage_group(cases[["Combined Summary Stage with Expanded Regional Codes (2004+)"]])
cases$time_months <- suppressWarnings(as.numeric(cases[["Survival Days"]])) / 30.4375
cases$event_os <- as.integer(cases[["Vital status recode (study cutoff used)"]] == "Dead")
cases$event_css <- as.integer(cases[["SEER cause-specific death classification"]] == "Dead (attributable to this cancer dx)")
css_class <- cases[["SEER cause-specific death classification"]]
other_class <- cases[["SEER other cause of death classification"]]
vital_class <- cases[["Vital status recode (study cutoff used)"]]
cases$cause_status <- ifelse(
  vital_class == "Alive", 0L,
  ifelse(css_class == "Dead (attributable to this cancer dx)", 1L,
         ifelse(other_class == "Dead (attributable to causes other than this cancer dx)", 2L, NA_integer_))
)
cases$hist_positive <- cases[["Diagnostic Confirmation"]] == "Positive histology"
cases$one_primary <- cases[["Sequence number"]] == "One primary only"
cases$age_band <- cut(cases$age_num, breaks=c(-Inf,49,64,74,Inf), labels=c("<50","50-64","65-74","75+"))
cases$sex2 <- factor(cases[["Sex"]], levels=c("Female","Male"))
race <- cases[["Race recode (with detailed Asian and Native Hawaiian other PI)"]]
# Sparse detailed race levels generated separation in bootstrap Cox samples;
# collapse them before modeling while retaining the original export unchanged.
cases$race_group <- factor(ifelse(race=="White","White","Other"), levels=c("White","Other"))
cases$hispanic <- factor(ifelse(grepl("Non-Spanish",cases[["Origin recode NHIA (Hispanic, Non-Hisp)"]]),"Non-Hispanic",
                                ifelse(grepl("Spanish|Hispanic",cases[["Origin recode NHIA (Hispanic, Non-Hisp)"]]),"Hispanic","Unknown")),
                         levels=c("Non-Hispanic","Hispanic","Unknown"))

gca <- cases[cases$histology %in% c(8243,8244,8245),]
composition <- as.data.frame(table(gca$period, gca$histology), stringsAsFactors=FALSE)
names(composition) <- c("period","histology","cases")
composition <- composition[composition$period != "" & composition$cases > 0,]
composition$period_total <- ave(composition$cases, composition$period, FUN=sum)
composition$share_percent <- 100*composition$cases/composition$period_total
write_csv(composition, "coding_composition_by_period_formal.csv")

stage_rows <- list(); k <- 1
for (dn in names(definitions)) {
  z <- cases[cases$year >= 2004 & cases$year <= 2023 & cases$histology %in% definitions[[dn]],]
  tb <- as.data.frame(table(z$stage), stringsAsFactors=FALSE); names(tb) <- c("stage","cases")
  tb$definition <- dn; tb$percent <- 100*tb$cases/sum(tb$cases); stage_rows[[k]] <- tb; k <- k+1
}
stage_out <- do.call(rbind, stage_rows)
write_csv(stage_out, "stage_distribution_2004_2023.csv")

grade_missing <- function(z) {
  old <- z$year <= 2017
  miss_old <- grepl("Unknown|Blank", z[["Grade Recode (thru 2017)"]], ignore.case=TRUE)
  miss_new <- grepl("Unknown|Blank", z[["Derived Summary Grade 2018 (2018+)"]], ignore.case=TRUE)
  mean(ifelse(old,miss_old,miss_new),na.rm=TRUE)*100
}
missing_rows <- list(); k <- 1
for (h in c(8243,8244,8245)) for (pn in levels(gca$period)) {
  z <- gca[gca$histology==h & gca$period==pn,]
  if (!nrow(z)) next
  missing_rows[[k]] <- data.frame(histology=h,period=pn,n=nrow(z),median_age=median(z$age_num,na.rm=TRUE),
    female_percent=100*mean(z[["Sex"]]=="Female"),positive_histology_percent=100*mean(z$hist_positive),
    stage_unknown_or_blank_percent=100*mean(z$stage=="Unknown"),grade_missing_percent=grade_missing(z))
  k <- k+1
}
write_csv(do.call(rbind,missing_rows), "coding_period_case_characteristics_and_missingness.csv")

km_estimates <- function(z, event_name, targets=c(12,36,60)) {
  if (!nrow(z)) return(NULL)
  sf <- survfit(Surv(time_months, z[[event_name]]) ~ 1, data=z, conf.type="log-log")
  ss <- summary(sf, times=targets, extend=TRUE)
  data.frame(month=targets, survival=ss$surv, lcl_95=ss$lower, ucl_95=ss$upper)
}
surv_rows <- list(); k <- 1
base_surv <- cases[cases$year <= 2018 & is.finite(cases$time_months),]
for (dn in names(definitions)) for (ep in c("event_os","event_css")) {
  z <- base_surv[base_surv$histology %in% definitions[[dn]],]
  q <- km_estimates(z,ep); q$definition <- dn; q$endpoint <- ifelse(ep=="event_os","OS","CSS"); q$n <- nrow(z)
  surv_rows[[k]] <- q; k <- k+1
}
km_def <- do.call(rbind,surv_rows); write_csv(km_def,"KM_survival_nested_definitions_2000_2018.csv")

sens_specs <- list(
  Base_2000_2018 = function(z) z$year<=2018,
  Positive_histology = function(z) z$year<=2018 & z$hist_positive,
  One_primary_only = function(z) z$year<=2018 & z$one_primary,
  Common_stage_era_2004_2018 = function(z) z$year>=2004 & z$year<=2018,
  Positive_histology_one_primary = function(z) z$year<=2018 & z$hist_positive & z$one_primary
)
sens_rows <- list(); k <- 1
for (sn in names(sens_specs)) for (dn in names(definitions)) for (ep in c("event_os","event_css")) {
  keep <- sens_specs[[sn]](cases) & cases$histology %in% definitions[[dn]] & is.finite(cases$time_months)
  z <- cases[keep,]; q <- km_estimates(z,ep,60)
  sens_rows[[k]] <- data.frame(sensitivity=sn,definition=dn,endpoint=ifelse(ep=="event_os","OS","CSS"),n=nrow(z),
    survival_5y=q$survival,lcl_95=q$lcl_95,ucl_95=q$ucl_95); k<-k+1
}
write_csv(do.call(rbind,sens_rows),"survival_definition_sensitivity_5y.csv")

# Exclusive-histology Cox models; 8244 and 8245 are comparisons, not reclassified GCA.
coxdat <- gca[gca$year>=2004 & gca$year<=2018 & gca$hist_positive & is.finite(gca$time_months) & is.finite(gca$age_num),]
coxdat$hist_group <- factor(as.character(coxdat$histology), levels=c("8243","8244","8245"))
coxdat$year_c <- (coxdat$year-2010)/10
cox_formula <- function(event) as.formula(paste0("Surv(time_months,",event,") ~ hist_group + sex2 + race_group + stage + year_c + strata(age_band)"))
fit_os <- coxph(cox_formula("event_os"), data=coxdat, x=TRUE, model=TRUE, ties="efron")
fit_css <- coxph(cox_formula("event_css"), data=coxdat, x=TRUE, model=TRUE, ties="efron")
tidy_cox <- function(fit, endpoint) {
  s <- summary(fit); cf <- as.data.frame(s$coefficients); ci <- as.data.frame(s$conf.int)
  data.frame(endpoint=endpoint,term=rownames(cf),HR=ci[["exp(coef)"]],LCL_95=ci[["lower .95"]],UCL_95=ci[["upper .95"]],p_value=cf[["Pr(>|z|)"]],row.names=NULL)
}
cox_out <- rbind(tidy_cox(fit_os,"OS"),tidy_cox(fit_css,"CSS"))
write_csv(cox_out,"exclusive_histology_adjusted_cox_models.csv")
ph_out <- rbind(data.frame(endpoint="OS",term=rownames(cox.zph(fit_os)$table),cox.zph(fit_os)$table,row.names=NULL),
                data.frame(endpoint="CSS",term=rownames(cox.zph(fit_css)$table),cox.zph(fit_css)$table,row.names=NULL))
write_csv(ph_out,"cox_proportional_hazards_tests.csv")

# Fine-Gray subdistribution hazards and prespecified effect-modification tests.
# Deaths with missing/unknown cause are excluded from the primary competing-risk analysis.
coxdat$age65 <- factor(ifelse(coxdat$age_num < 65, "<65", ">=65"), levels=c("<65", ">=65"))
coxdat$int_8244_age65 <- as.integer(coxdat$hist_group == "8244" & coxdat$age65 == ">=65")
coxdat$int_8245_age65 <- as.integer(coxdat$hist_group == "8245" & coxdat$age65 == ">=65")
coxdat$case_id <- seq_len(nrow(coxdat))
fgdat <- coxdat[!is.na(coxdat$cause_status), ]
fgdat$event_mstate <- factor(fgdat$cause_status, levels=0:2,
                             labels=c("censor", "cancer_death", "other_death"))

event_counts <- as.data.frame(table(fgdat$hist_group, fgdat$event_mstate), stringsAsFactors=FALSE)
names(event_counts) <- c("histology", "event", "n")
unknown_counts <- aggregate(case_id ~ hist_group, coxdat[is.na(coxdat$cause_status), ], length)
names(unknown_counts) <- c("histology", "n")
unknown_counts$event <- "unknown_cause_death_excluded"
event_counts <- rbind(event_counts, unknown_counts[,c("histology","event","n")])
write_csv(event_counts, "finegray_event_counts.csv")

fit_finegray <- function(rhs) {
  fg_formula <- as.formula(paste0("Surv(time_months, event_mstate) ~ ", rhs, " + case_id"))
  expanded <- finegray(fg_formula, data=fgdat, etype="cancer_death", id=case_id)
  fit_formula <- as.formula(paste0("Surv(fgstart, fgstop, fgstatus) ~ ", rhs, " + cluster(case_id)"))
  fit <- coxph(fit_formula, data=expanded, weights=fgwt, ties="efron", x=TRUE, model=TRUE)
  list(fit=fit, expanded=expanded)
}

fg_main <- fit_finegray("hist_group + age_band + sex2 + race_group + stage + year_c")
fg_out <- tidy_cox(fg_main$fit, "Cancer-specific death, Fine-Gray")
write_csv(fg_out, "finegray_adjusted_subdistribution_hazards.csv")

# Independent implementation check using cmprsk::crr on the same complete-case
# design matrix and event coding (1=cancer death, 2=other death, 0=censored).
if (requireNamespace("cmprsk", quietly=TRUE)) {
  mm <- model.matrix(~ hist_group + age_band + sex2 + race_group + stage + year_c, data=fgdat)[,-1,drop=FALSE]
  crr_fit <- cmprsk::crr(ftime=fgdat$time_months, fstatus=fgdat$cause_status,
                         cov1=mm, failcode=1, cencode=0)
  crr_se <- sqrt(diag(crr_fit$var))
  crr_check <- data.frame(term=names(crr_fit$coef), sHR=exp(crr_fit$coef),
                          LCL_95=exp(crr_fit$coef-1.96*crr_se),
                          UCL_95=exp(crr_fit$coef+1.96*crr_se),
                          p_value=2*pnorm(-abs(crr_fit$coef/crr_se)), row.names=NULL)
  write_csv(crr_check, "finegray_cmprsk_crosscheck.csv")
}

joint_wald <- function(fit, pattern, endpoint, modifier) {
  ix <- grep(pattern, names(coef(fit)), fixed=TRUE)
  if (!length(ix)) stop("No interaction coefficients matched: ", pattern)
  b <- coef(fit)[ix]
  V <- vcov(fit)[ix,ix,drop=FALSE]
  message("Interaction terms [",endpoint,"; ",modifier,"]: ",paste(names(b),collapse=", "),
          "; covariance rank=",qr(V)$rank,"/",length(ix))
  stat <- as.numeric(t(b) %*% solve(V, b))
  data.frame(endpoint=endpoint, modifier=modifier, test="Global interaction Wald test",
             df=length(ix), statistic=stat, p_value=pchisq(stat,length(ix),lower.tail=FALSE))
}

contrast_rows <- function(fit, endpoint, modifier, levels_modifier, interaction_suffix) {
  V <- vcov(fit); b <- coef(fit); out <- list(); k <- 1
  for (h in c("8244","8245")) {
    base_name <- paste0("hist_group",h)
    for (lev in levels_modifier) {
      L <- setNames(rep(0,length(b)),names(b)); L[base_name] <- 1
      if (lev != levels_modifier[1]) {
        candidates <- if (interaction_suffix == "explicit_age65") {
          paste0("int_", h, "_age65")
        } else {
          c(paste0(base_name,":",interaction_suffix), paste0(interaction_suffix,":",base_name))
        }
        int_name <- candidates[candidates %in% names(b)][1]
        if (is.na(int_name)) stop("Interaction coefficient not found for ",base_name," and ",interaction_suffix)
        L[int_name] <- 1
      }
      est <- sum(L*b); se <- sqrt(as.numeric(t(L)%*%V%*%L))
      out[[k]] <- data.frame(endpoint=endpoint, modifier=modifier, subgroup=lev,
                             comparison=paste0(h," vs 8243"), HR=exp(est),
                             LCL_95=exp(est-1.96*se), UCL_95=exp(est+1.96*se)); k <- k+1
    }
  }
  do.call(rbind,out)
}

# OS interaction models. The two explicit age-interaction indicators avoid a redundant
# reference-group parameter while retaining four-band age-stratified baseline hazards.
fit_os_age_int <- coxph(Surv(time_months,event_os) ~ hist_group + int_8244_age65 + int_8245_age65 + sex2 + race_group + stage + year_c + strata(age_band),
                        data=coxdat, ties="efron", x=TRUE, model=TRUE)
fit_os_sex_int <- coxph(Surv(time_months,event_os) ~ hist_group*sex2 + race_group + stage + year_c + strata(age_band),
                        data=coxdat, ties="efron", x=TRUE, model=TRUE)

# Fine-Gray interaction models use the same competing event definition and covariate set.
fg_age_int <- fit_finegray("hist_group + int_8244_age65 + int_8245_age65 + age_band + sex2 + race_group + stage + year_c")
fg_sex_int <- fit_finegray("hist_group*sex2 + age_band + race_group + stage + year_c")

interaction_tests <- rbind(
  joint_wald(fit_os_age_int,"int_824","Overall survival","Age (<65 vs >=65)"),
  joint_wald(fit_os_sex_int,":sex2Male","Overall survival","Sex"),
  joint_wald(fg_age_int$fit,"int_824","Cancer-specific death, Fine-Gray","Age (<65 vs >=65)"),
  joint_wald(fg_sex_int$fit,":sex2Male","Cancer-specific death, Fine-Gray","Sex")
)
write_csv(interaction_tests,"histology_age_sex_global_interaction_tests.csv")

interaction_estimates <- rbind(
  contrast_rows(fit_os_age_int,"Overall survival","Age",c("<65",">=65"),"explicit_age65"),
  contrast_rows(fit_os_sex_int,"Overall survival","Sex",c("Female","Male"),"sex2Male"),
  contrast_rows(fg_age_int$fit,"Cancer-specific death, Fine-Gray","Age",c("<65",">=65"),"explicit_age65"),
  contrast_rows(fg_sex_int$fit,"Cancer-specific death, Fine-Gray","Sex",c("Female","Male"),"sex2Male")
)
write_csv(interaction_estimates,"histology_age_sex_interaction_stratum_estimates.csv")

standardized_5y <- function(fit, data, endpoint, B=200) {
  calc <- function(f, refdata) {
    bh <- basehaz(f, centered=FALSE)
    if ("strata" %in% names(bh)) {
      h60 <- tapply(seq_len(nrow(bh)), bh$strata, function(ii) max(bh$hazard[ii][bh$time[ii] <= 60], na.rm=TRUE))
      names(h60) <- sub("^.*=", "", names(h60))
      H0 <- unname(h60[match(as.character(refdata$age_band), names(h60))])
    } else H0 <- rep(max(bh$hazard[bh$time<=60],na.rm=TRUE), nrow(refdata))
    sapply(levels(data$hist_group), function(h) {
      nd <- refdata; nd$hist_group <- factor(h,levels=levels(data$hist_group))
      lp <- predict(f,newdata=nd,type="lp",reference="zero")
      mean(exp(-H0*exp(lp)),na.rm=TRUE)
    })
  }
  point <- calc(fit,data); boots <- matrix(NA_real_,B,length(point)); set.seed(8243 + nchar(endpoint)); boot_errors <- character()
  for (i in seq_len(B)) {
    ii <- sample.int(nrow(data),replace=TRUE); dd <- data[ii,]
    f <- try(coxph(formula(fit),data=dd,ties="efron",x=TRUE,model=TRUE),silent=TRUE)
    if (!inherits(f,"try-error")) boots[i,] <- tryCatch(calc(f,data),error=function(e) {boot_errors <<- c(boot_errors,conditionMessage(e)); rep(NA_real_,length(point))})
  }
  if (length(boot_errors)) writeLines(unique(boot_errors), file.path(out_dir,paste0("standardization_bootstrap_errors_",endpoint,".txt")))
  data.frame(endpoint=endpoint,histology=names(point),standardized_survival_5y=as.numeric(point),
             lcl_95=apply(boots,2,quantile,.025,na.rm=TRUE),ucl_95=apply(boots,2,quantile,.975,na.rm=TRUE),
             bootstrap_replicates=B)
}
std_surv <- rbind(standardized_5y(fit_os,coxdat,"OS"),standardized_5y(fit_css,coxdat,"CSS"))
write_csv(std_surv,"exclusive_histology_marginal_standardized_5y_survival.csv")

flow <- data.frame(
  step=c("All exported appendiceal tumors","Historical spectrum 8243/8244/8245","Positive histologic confirmation",
         "Survival cohort diagnosed 2000-2018","Adjusted cohort diagnosed 2004-2018 with covariates"),
  n=c(nrow(cases),nrow(gca),sum(gca$hist_positive),sum(gca$year<=2018 & is.finite(gca$time_months)),nrow(coxdat)))
write_csv(flow,"cohort_flow.csv")

# ---------------- Publication figures (R backend only) ----------------
pal <- c(D1_strict_GCA="#24557A",D2_legacy_GCC_adenocarcinoid="#D8873A",D3_historical_expanded_spectrum="#8E6C9E")
hist_pal <- c(`8243`="#24557A",`8244`="#8E6C9E",`8245`="#D8873A")
theme_pub <- theme_classic(base_size=8,base_family="Arial") + theme(axis.line=element_line(linewidth=.35),axis.ticks=element_line(linewidth=.35),
  legend.title=element_blank(),legend.position="top",strip.background=element_blank(),strip.text=element_text(face="bold"),plot.title=element_text(face="bold",size=9))
p1a <- ggplot(annual,aes(year,age_adjusted_rate_per_100k,color=definition,fill=definition))+
  geom_ribbon(aes(ymin=lcl_95,ymax=ucl_95),alpha=.10,color=NA)+geom_line(linewidth=.7)+geom_point(size=.8)+
  scale_color_manual(values=pal,labels=c("Strict 8243","8243 + 8245","8243 + 8244 + 8245"))+
  scale_fill_manual(values=pal,labels=c("Strict 8243","8243 + 8245","8243 + 8244 + 8245"))+
  labs(x="Year of diagnosis",y="Age-adjusted rate per 100,000",title="a  Reported incidence depends on case definition")+theme_pub+
  theme(legend.position="bottom",plot.margin=margin(4,4,10,4))
comp_plot <- composition; comp_plot$histology <- factor(comp_plot$histology,levels=c("8245","8244","8243"))
p1b <- ggplot(comp_plot,aes(period,share_percent,fill=histology))+geom_col(width=.72)+
  scale_fill_manual(values=hist_pal,breaks=c("8243","8244","8245"),labels=c("8243 strict GCA","8244 mixed","8245 adenocarcinoid"))+
  labs(x="Diagnosis period",y="Share of expanded spectrum (%)",title="b  Coding composition")+theme_pub+theme(axis.text.x=element_text(angle=25,hjust=1),legend.position="bottom",plot.margin=margin(8,4,4,4))
infl <- period_rates[period_rates$period!="2000-2023" & period_rates$definition!="D1_strict_GCA",]
p1c <- ggplot(infl,aes(period,rate_inflation_vs_D1_percent,color=definition,group=definition))+
  geom_hline(yintercept=0,color="grey70")+geom_errorbar(aes(ymin=inflation_lcl_95,ymax=inflation_ucl_95),width=.12,linewidth=.35)+geom_point(size=1.8)+geom_line(linewidth=.5)+
  scale_color_manual(values=pal[-1],labels=c("8243 + 8245","Expanded spectrum"))+
  labs(x="Diagnosis period",y="Rate change vs strict 8243 (%)",title="c  Definition effect over time")+theme_pub+theme(axis.text.x=element_text(angle=25,hjust=1),legend.position="bottom",plot.margin=margin(8,4,4,4))
fig1 <- p1a / (p1b | p1c) + plot_layout(heights=c(1.15,1))

km_long <- list(); km_risk <- list(); kk<-1
risk_times <- seq(0,120,24)
for(h in levels(coxdat$hist_group)) {
  z<-coxdat[coxdat$hist_group==h,]; sf<-survfit(Surv(time_months,event_os)~1,data=z); ss<-summary(sf)
  km_long[[kk]]<-data.frame(time=ss$time,survival=ss$surv,histology=h)
  sr<-summary(sf,times=risk_times,extend=TRUE)
  km_risk[[kk]]<-data.frame(time=risk_times,n_risk=sr$n.risk,histology=h)
  kk<-kk+1
}
km_long<-do.call(rbind,km_long)
km_risk<-do.call(rbind,km_risk)
km_risk$histology<-factor(km_risk$histology,levels=c("8245","8244","8243"))
km_risk$hjust<-ifelse(km_risk$time==0,0,ifelse(km_risk$time==120,1,.5))
write_csv(km_risk,"exclusive_histology_numbers_at_risk.csv")
p2a<-ggplot(km_long,aes(time,survival,color=histology))+geom_step(linewidth=.75)+coord_cartesian(ylim=c(.2,1))+
  scale_x_continuous(limits=c(0,120),breaks=risk_times,expand=expansion(mult=c(0,.01)))+
  scale_color_manual(values=hist_pal,labels=c("8243 strict GCA","8244 mixed","8245 adenocarcinoid"))+
  labs(x=NULL,y="Overall survival",title="a  Exclusive-code survival differs")+theme_pub+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position="top",plot.margin=margin(4,4,0,4))
p2risk<-ggplot(km_risk,aes(time,histology,label=n_risk,color=histology,hjust=hjust))+geom_text(size=2.25,show.legend=FALSE)+
  scale_x_continuous(limits=c(0,120),breaks=risk_times,expand=expansion(mult=c(0,.01)))+
  scale_y_discrete(labels=c(`8245`="8245",`8244`="8244",`8243`="8243"))+
  scale_color_manual(values=hist_pal)+labs(x="Months since diagnosis",y="No. at risk")+theme_pub+
  theme(legend.position="none",axis.line.y=element_blank(),axis.ticks.y=element_blank(),
        panel.grid.major.x=element_line(color="grey90",linewidth=.25),plot.margin=margin(0,4,4,4))
forest<-cox_out[cox_out$endpoint=="OS" & grepl("hist_group",cox_out$term),]
forest$label<-c("8244 vs 8243","8245 vs 8243")[match(forest$term,c("hist_group8244","hist_group8245"))]
p2b<-ggplot(forest,aes(HR,label))+geom_vline(xintercept=1,color="grey60",linetype=2)+geom_errorbar(aes(xmin=LCL_95,xmax=UCL_95),width=.15,orientation="y")+geom_point(size=2,color="#24557A")+
  scale_x_log10()+labs(x="Adjusted hazard ratio (log scale)",y=NULL,title="b  Covariate-adjusted OS")+theme_pub+theme(legend.position="none")
stdp<-std_surv[std_surv$endpoint=="OS",]; stdp$histology<-factor(stdp$histology,levels=c("8243","8244","8245"))
p2c<-ggplot(stdp,aes(histology,100*standardized_survival_5y,fill=histology))+geom_col(width=.6)+geom_errorbar(aes(ymin=100*lcl_95,ymax=100*ucl_95),width=.15)+
  scale_fill_manual(values=hist_pal)+labs(x="Exclusive histology code",y="Standardized 5-year OS (%)",title="c  Common-covariate survival")+theme_pub+theme(legend.position="none")
p2left <- p2a / p2risk + plot_layout(heights=c(3.4,1.15))
fig2 <- p2left | (p2b / p2c) + plot_layout(widths=c(1.35,1))

# Supplementary Figure S1: competing-risk estimate and prespecified exploratory
# effect-modification analyses. Global interaction P values are shown in facets.
fg_forest <- fg_out[grepl("hist_group",fg_out$term),]
fg_forest$comparison <- c("8244 vs 8243","8245 vs 8243")[match(fg_forest$term,c("hist_group8244","hist_group8245"))]
pS1a <- ggplot(fg_forest,aes(HR,comparison,color=comparison))+
  geom_vline(xintercept=1,color="grey60",linetype=2)+
  geom_errorbar(aes(xmin=LCL_95,xmax=UCL_95),width=.16,orientation="y",linewidth=.45)+
  geom_point(size=2)+scale_x_log10()+scale_color_manual(values=c("8244 vs 8243"="#8E6C9E","8245 vs 8243"="#D8873A"))+
  labs(x="Adjusted subdistribution hazard ratio (log scale)",y=NULL,
       title="a  Cancer-specific death with other death competing")+theme_pub+theme(legend.position="none")

int_plot <- interaction_estimates
int_plot$endpoint_short <- ifelse(int_plot$endpoint=="Overall survival","OS","Fine-Gray")
global_lookup <- interaction_tests
global_lookup$endpoint_short <- ifelse(global_lookup$endpoint=="Overall survival","OS","Fine-Gray")
facet_lab <- setNames(sprintf("%s (P-interaction = %.3f)",global_lookup$endpoint_short,global_lookup$p_value),
                      paste(global_lookup$modifier,global_lookup$endpoint_short,sep="|"))
make_interaction_panel <- function(mod,title) {
  z <- int_plot[int_plot$modifier==mod,]
  z$subgroup <- ifelse(z$subgroup==">=65","65 or older",z$subgroup)
  if (mod == "Age") z$subgroup <- factor(z$subgroup, levels=c("<65","65 or older"))
  z$facet <- facet_lab[paste(ifelse(mod=="Age","Age (<65 vs >=65)","Sex"),z$endpoint_short,sep="|")]
  ggplot(z,aes(subgroup,HR,color=comparison,group=comparison))+
    geom_hline(yintercept=1,color="grey60",linetype=2)+
    geom_errorbar(aes(ymin=LCL_95,ymax=UCL_95),position=position_dodge(width=.42),width=.15,linewidth=.4)+
    geom_point(position=position_dodge(width=.42),size=1.8)+scale_y_log10()+
    scale_color_manual(values=c("8244 vs 8243"="#8E6C9E","8245 vs 8243"="#D8873A"))+
    facet_wrap(~facet,nrow=1)+labs(x=NULL,y="Adjusted HR or sHR (log scale)",title=title)+
    theme_pub+theme(legend.position="bottom",axis.text.x=element_text(angle=18,hjust=1))
}
pS1b <- make_interaction_panel("Age","b  Age-stratified estimates")
pS1c <- make_interaction_panel("Sex","c  Sex-stratified estimates")
figS1 <- pS1a / pS1b / pS1c + plot_layout(heights=c(.8,1,1))

save_pub <- function(plot,stem,width_mm=183,height_mm=145,dpi=600) {
  w<-width_mm/25.4;h<-height_mm/25.4
  # Some Windows Cairo devices cannot translate non-ASCII output paths under
  # the C locale. Render to an ASCII-safe temporary directory first and then
  # copy the completed files into the project output directory.
  tmp_dir <- file.path(tempdir(), "gca_publication_figures")
  dir.create(tmp_dir, recursive=TRUE, showWarnings=FALSE)
  tmp_svg <- file.path(tmp_dir,paste0(stem,".svg"))
  tmp_pdf <- file.path(tmp_dir,paste0(stem,".pdf"))
  tmp_png <- file.path(tmp_dir,paste0(stem,".png"))
  tmp_tiff <- file.path(tmp_dir,paste0(stem,".tiff"))
  svglite(tmp_svg,width=w,height=h);print(plot);dev.off()
  cairo_pdf(tmp_pdf,width=w,height=h,family="Arial");print(plot);dev.off()
  png(tmp_png,width=w,height=h,units="in",res=dpi,type="cairo");print(plot);dev.off()
  tiff(tmp_tiff,width=w,height=h,units="in",res=dpi,compression="lzw",type="cairo");print(plot);dev.off()
  for (src in c(tmp_svg,tmp_pdf,tmp_png,tmp_tiff)) {
    ok <- file.copy(src,file.path(fig_dir,basename(src)),overwrite=TRUE)
    if (!ok) stop("Failed to copy publication figure: ",src)
  }
}
save_pub(fig1,"Figure1_case_definition_and_coding_migration",183,165)
save_pub(fig2,"Figure2_survival_composition_and_adjustment",183,132)
save_pub(figS1,"Supplementary_Figure_S1_competing_risk_and_interactions",183,190)

cat("Formal next-phase analysis completed. Outputs:",out_dir,"\n")
