#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript run_methods_extension.R <input_csv> <output_dir>")
options(warn = as.integer(Sys.getenv("WARN_LEVEL", "0")))
input_csv <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
out <- args[[2]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

set.seed(20261011)
d <- read.csv(input_csv, check.names = FALSE, stringsAsFactors = FALSE,
              fileEncoding = "UTF-8-BOM")
get_series <- function(label) d[d$Definition == label, c("Year", "Rate", "Standard Error", "Count", "Population")]
icc0 <- get_series("ICC-observed")
icc1 <- get_series("ICC-corrected")
ecc0 <- get_series("ECC-observed")
ecc1 <- get_series("ECC-corrected")
ki <- get_series("K-intra")
# Independent disjoint components: A=ICC without K-intra, B=observed ECC, K=K-intra.
# Deterministic fractional weights, not Bernoulli reassignment of individual records.
stopifnot(identical(icc0$Year, ki$Year), identical(ecc0$Year, ki$Year))
audit <- data.frame(check=c("ICC rate additivity", "ECC rate additivity", "ICC variance additivity", "ECC variance additivity"),
 max_error=c(max(abs(icc0$Rate-icc1$Rate-ki$Rate)), max(abs(ecc1$Rate-ecc0$Rate-ki$Rate)),
 max(abs(icc0$`Standard Error`^2-icc1$`Standard Error`^2-ki$`Standard Error`^2)),
 max(abs(ecc1$`Standard Error`^2-ecc0$`Standard Error`^2-ki$`Standard Error`^2))))
audit$tolerance <- c(2e-6,2e-6,1e-5,1e-5)
audit$pass <- audit$max_error <= audit$tolerance
write.csv(audit,file.path(out,"partial_variance_component_audit.csv"),row.names=FALSE)
stopifnot(all(audit$pass))

fit_apc <- function(year, rate, se) {
  ok <- is.finite(rate) & rate > 0 & is.finite(se) & se > 0
  z <- data.frame(year = year[ok], rate = rate[ok], se = se[ok])
  m <- lm(log(rate) ~ I(year - min(year)), data = z, weights = rate^2 / se^2)
  b <- coef(m)[2]
  ci <- confint(m)[2, ]
  c(APC = 100 * (exp(b) - 1), LCL = 100 * (exp(ci[1]) - 1), UCL = 100 * (exp(ci[2]) - 1))
}

fractions <- c(0, .25, .50, .75, 1)
partial_long <- do.call(rbind, lapply(fractions, function(p) {
  do.call(rbind, lapply(c("ICC", "ECC"), function(site) {
    a <- if (site == "ICC") icc0 else ecc0
    b <- if (site == "ICC") icc1 else ecc1
    data.frame(site = site, fraction = p, year = a$Year,
               rate = a$Rate + p * (b$Rate - a$Rate),
               se = if(site=="ICC") sqrt(icc1$`Standard Error`^2+(1-p)^2*ki$`Standard Error`^2) else sqrt(ecc0$`Standard Error`^2+p^2*ki$`Standard Error`^2),
               count = a$Count + p * (b$Count - a$Count))
  }))
}))
partial_summary <- do.call(rbind, lapply(split(partial_long, interaction(partial_long$site, partial_long$fraction)), function(z) {
  est <- fit_apc(z$year, z$rate, z$se)
  data.frame(site = z$site[1], redistributed_percent = 100 * z$fraction[1],
             cases_1992_2023 = round(sum(z$count)), APC = est[1], LCL = est[2], UCL = est[3])
}))
partial_summary <- partial_summary[order(partial_summary$site, partial_summary$redistributed_percent), ]
write.csv(partial_long, file.path(out, "partial_redistribution_annual_series.csv"), row.names = FALSE)
write.csv(partial_summary, file.path(out, "partial_redistribution_loglinear_APC.csv"), row.names = FALSE)

# Scenario-based stress test. This uses a transparent BIC screen for large-scale
# mechanism exploration, not the NCI Joinpoint permutation procedure used for
# the empirical analyses.
years <- 2000:2023
tt <- years - 2000
candidate_breaks <- c(2004, 2008, 2012, 2016, 2020)
n_rep <- as.integer(Sys.getenv("SIM_REPS", "2000"))
if (!is.finite(n_rep) || n_rep < 1L) stop("SIM_REPS must be a positive integer")
nonnull_grid <- expand.grid(
  true_apc = c(0, 2, 5), mean_events_2000 = c(25, 100, 400),
  mechanism = c("over-inclusion", "under-ascertainment"),
  magnitude = c(.10, .25, .50), shape = c("linear", "step", "logistic"),
  transition_end = c(2010, 2018), stringsAsFactors = FALSE
)
null_grid <- expand.grid(
  true_apc = c(0, 2, 5), mean_events_2000 = c(25, 100, 400),
  stringsAsFactors = FALSE
)
null_grid$mechanism <- "none"
null_grid$magnitude <- 0
null_grid$shape <- "none"
null_grid$transition_end <- NA_integer_
grid <- rbind(nonnull_grid, null_grid)
grid$scenario_id <- sprintf("S%03d", seq_len(nrow(grid)))

profile <- function(level, shape, transition_end) {
  if (level == 0) return(rep(0, length(years)))
  if (shape == "step") return(ifelse(years < transition_end, level, 0))
  if (shape == "linear") {
    return(level * pmax(transition_end - years, 0) / (transition_end - min(years)))
  }
  if (shape == "logistic") {
    midpoint <- (min(years) + transition_end) / 2
    scale <- max((transition_end - min(years)) / 6, 1)
    raw <- plogis((midpoint - years) / scale)
    start_value <- plogis((midpoint - min(years)) / scale)
    end_value <- plogis((midpoint - transition_end) / scale)
    normalized <- (raw - end_value) / (start_value - end_value)
    return(level * pmin(1, pmax(0, normalized)))
  }
  stop("Unknown transition shape: ", shape)
}

simulate_counts <- function(apc, base, mechanism, magnitude, shape, transition_end) {
  truth <- base * (1 + apc / 100)^tt
  if (mechanism == "none") {
    mu <- truth
  } else {
    distortion <- profile(magnitude, shape, transition_end)
    mu <- if (mechanism == "over-inclusion") truth * (1 + distortion) else truth * (1 - distortion)
  }
  rpois(length(years), mu)
}

designs <- c(list(cbind(1,tt)),lapply(candidate_breaks,function(k) cbind(1,tt,pmax(0,years-k))))
qs <- lapply(designs,function(x) qr.Q(qr(x)))
slope_coef <- solve(crossprod(designs[[1]]),t(designs[[1]]))[2,]
screen_counts <- function(y, apc, transition_end = NA_integer_) {
  ly <- log(y + 0.5)
  rss <- vapply(qs,function(q) sum((ly-q%*%crossprod(q,ly))^2),numeric(1))
  bic <- length(y)*log(rss/length(y))+c(2,rep(4,length(candidate_breaks)))*log(length(y))
  bic0 <- bic[1]; best_bic <- min(bic[-1]); best_break <- candidate_breaks[which.min(bic[-1])]
  est <- unname(100 * (exp(sum(slope_coef*ly)) - 1))
  c(estimated_apc = est, bias = est - apc,
    break_selected = as.integer(best_bic < bic0),
    break_near_transition = as.integer(best_bic < bic0 && is.finite(transition_end) && abs(best_break - transition_end) <= 2),
    selected_break = ifelse(best_bic < bic0, best_break, NA))
}

one_rep <- function(apc, base, mechanism, magnitude, shape, transition_end) {
  y <- simulate_counts(apc, base, mechanism, magnitude, shape, transition_end)
  screen_counts(y, apc, transition_end)
}

rows <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  g <- grid[i, , drop = FALSE]
  rownames(g) <- NULL
  sims <- replicate(n_rep, one_rep(g$true_apc, g$mean_events_2000, g$mechanism,
                                   g$magnitude, g$shape, g$transition_end))
  sim_df <- as.data.frame(t(sims))
  rownames(sim_df) <- NULL
  scenario_df <- g[rep(1L, n_rep), , drop = FALSE]
  rownames(scenario_df) <- NULL
  rows[[i]] <- cbind(scenario_df, replicate = seq_len(n_rep), sim_df)
  if(i %% 30==0) cat("Completed scenario",i,"of",nrow(grid),"\n")
}
sim <- do.call(rbind, rows)
q <- function(x, p) unname(quantile(x, p, na.rm = TRUE, type = 8))
sim_summary <- do.call(rbind, lapply(split(sim, seq_len(nrow(grid))[rep(seq_len(nrow(grid)), each = n_rep)]), function(z) {
  data.frame(true_apc = z$true_apc[1], mean_events_2000 = z$mean_events_2000[1],
             mechanism = z$mechanism[1], distortion_percent = 100 * z$magnitude[1], shape = z$shape[1],
             transition_end = z$transition_end[1], simulations = nrow(z),
             median_estimated_apc = median(z$estimated_apc),
             apc_p025 = q(z$estimated_apc, .025), apc_p975 = q(z$estimated_apc, .975),
             median_bias = median(z$bias),
             breakpoint_selected_percent = 100 * mean(z$break_selected),
             breakpoint_near_transition_percent = 100 * mean(z$break_near_transition))
}))
write.csv(sim_summary, file.path(out, "expanded_simulation_summary.csv"), row.names = FALSE)
write.csv(sim, file.path(out, "expanded_simulation_replicates.csv"), row.names = FALSE)
write.csv(data.frame(parameter = c("calendar years", "true APC (%)", "mean events in 2000", "classification mechanism", "initial distortion", "transition shape", "transition end", "replicates", "candidate breakpoint years", "selection"),
                     value = c("2000-2023", "0, 2, 5", "25, 100, 400", "over-inclusion, under-ascertainment, or none", "0%, 10%, 25%, 50%", "linear, step, or logistic", "2010 or 2018", n_rep, "2004, 2008, 2012, 2016, 2020", "Gaussian log-count BIC screen: zero versus one breakpoint")),
          file.path(out, "expanded_simulation_parameters.csv"), row.names = FALSE)

# Recreate the identical calibration counts from their independent seed.
# Official GUI batch outputs for these 2000 series are retained separately.
calibration_scenarios <- expand.grid(
  true_apc = 5, mean_events_2000 = 400,
  mechanism = c("over-inclusion", "under-ascertainment"),
  magnitude = c(.10, .50), shape = c("linear", "logistic"),
  transition_end = 2018, stringsAsFactors = FALSE
)
calibration_scenarios$calibration_id <- sprintf("C%02d", seq_len(nrow(calibration_scenarios)))
calibration_reps <- as.integer(Sys.getenv("CALIBRATION_REPS", "250"))
set.seed(20261012)
calibration_rows <- vector("list", nrow(calibration_scenarios) * calibration_reps)
calibration_screen <- vector("list", length(calibration_rows))
idx <- 1L
for (i in seq_len(nrow(calibration_scenarios))) {
  g <- calibration_scenarios[i, , drop = FALSE]
  for (r in seq_len(calibration_reps)) {
    y <- simulate_counts(g$true_apc, g$mean_events_2000, g$mechanism,
                         g$magnitude, g$shape, g$transition_end)
    series_id <- paste0(g$calibration_id, "_R", sprintf("%03d", r))
    calibration_rows[[idx]] <- data.frame(
      series_id = series_id, calibration_id = g$calibration_id,
      replicate = r, year = years, count = y, log_count = log(y + 0.5)
    )
    ss <- screen_counts(y, g$true_apc, g$transition_end)
    calibration_screen[[idx]] <- data.frame(
      series_id = series_id, calibration_id = g$calibration_id,
      replicate = r, t(ss), check.names = FALSE
    )
    idx <- idx + 1L
  }
}
calibration_long <- do.call(rbind, calibration_rows)
calibration_bic <- do.call(rbind, calibration_screen)
write.csv(calibration_scenarios, file.path(out, "joinpoint_calibration_scenarios.csv"), row.names = FALSE)
write.csv(calibration_long, file.path(out, "joinpoint_calibration_input_long.csv"), row.names = FALSE)
write.csv(calibration_bic, file.path(out, "joinpoint_calibration_bic_results.csv"), row.names = FALSE)

# Figure 1: reusable four-stage framework and decision guide.
tiff(file.path(out, "Figure_1_Framework.tiff"), width = 4320, height = 3480, res = 600, compression = "lzw")
par(mar = c(0,0,0,0)); plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1))
cols <- c("#DCEAF7", "#E5F2EA", "#FFF0D5", "#F2E3F3")
titles <- c("1  Prespecify definitions", "2  Quantify sensitivity", "3  Audit migration if applicable", "4  Stress-test mechanisms")
body <- c("Contemporary primary definition\nHistorically plausible increments\nExact topography and morphology codes",
          "Counts, rates and code composition\nTrend magnitude and model structure\nPeriod-specific definition effects",
          "Only for cross-category migration\nPartial and complete reallocation\nAudit count-preserving identities",
          "Known-truth Poisson scenarios\nOver-inclusion and under-ascertainment\nBias and breakpoint selection")
xs <- c(.26,.74,.74,.26); ys <- c(.79,.79,.53,.53)
for(i in 1:4){
  rect(xs[i]-.20,ys[i]-.105,xs[i]+.20,ys[i]+.105,col=cols[i],border="#35556B",lwd=2)
  text(xs[i],ys[i]+.062,titles[i],font=2,cex=.70,col="#173A4F")
  text(xs[i],ys[i]-.022,body[i],cex=.52,col="#23343D")
}
arrows(.47,.79,.52,.79,length=.07,lwd=2,col="#557785")
arrows(.74,.68,.74,.64,length=.07,lwd=2,col="#557785")
arrows(.53,.53,.48,.53,length=.07,lwd=2,col="#557785")
text(.50,.375,"When is the full audit warranted?",font=2,cex=.72,col="#173A4F")
decision_x <- c(.20,.50,.80)
decision_text <- c("No plausible historical\ndefinition change\nCore disclosure may suffice",
                   "Coding update or disputed term\nApply Stages 1-2; add Stage 3\nonly for category migration",
                   "Conclusion depends on a\nbreakpoint or changing slope\nAdd Stage 4 stress tests")
decision_col <- c("#F5F7F8", "#FFF0D5", "#F2E3F3")
for(i in 1:3){
  rect(decision_x[i]-.13,.21,decision_x[i]+.13,.33,col=decision_col[i],border="#557785",lwd=1.7)
  text(decision_x[i],.27,decision_text[i],cex=.48,col="#23343D")
}
arrows(.34,.27,.36,.27,length=.06,lwd=1.7,col="#557785")
arrows(.64,.27,.66,.27,length=.06,lwd=1.7,col="#557785")
text(.50,.125,"Output: separate level, slope, structure and interpretation profiles",font=2,cex=.68,col="#173A4F")
text(.50,.085,"Boundary: evaluates dependence on operational definitions; does not reconstruct record-level truth.",cex=.52,col="#4A5960")
dev.off()

# Figure 3: partial redistribution.
tiff(file.path(out, "Figure_3_Partial_Redistribution.tiff"), width = 4320, height = 2700, res = 600, compression = "lzw")
par(mfrow=c(1,2), mar=c(5.2,5.1,2.6,1.5), mgp=c(3.0,.8,0), las=1)
site_cols <- c(ICC="#2B6F9C", ECC="#C4572D")
for(site in c("ICC","ECC")){
  z <- partial_summary[partial_summary$site==site,]
  plot(z$redistributed_percent,z$APC,type="b",pch=16,lwd=2.5,col=site_cols[site],
       xlab="Redistributed records (%)",ylab="Full-period log-linear APC (%)",
       main=ifelse(site=="ICC","A  Intrahepatic cholangiocarcinoma","B  Extrahepatic cholangiocarcinoma"),
       ylim=range(c(z$LCL,z$UCL)), cex.main=.72)
  arrows(z$redistributed_percent,z$LCL,z$redistributed_percent,z$UCL,angle=90,code=3,length=.04,col=site_cols[site])
  abline(h=0,lty=3,col="grey50")
}
dev.off()

# Figure 4: mechanism direction, transition shape, and information size.
tiff(file.path(out, "Figure_4_Simulation_Stress_Test.tiff"), width = 4320, height = 3000, res = 600, compression = "lzw")
par(mfrow=c(2,2), mar=c(4.4,4.7,2.2,1.0), mgp=c(2.7,.8,0), las=1)
shape_names <- c("linear","step","logistic")
shape_cols <- c("#2B6F9C","#C4572D","#6B5B95")
shape_pch <- c(16,17,15)
null400 <- sim_summary[sim_summary$true_apc==5 & sim_summary$mean_events_2000==400 & sim_summary$mechanism=="none",]
for(mech in c("over-inclusion","under-ascertainment")){
  z <- sim_summary[sim_summary$true_apc==5 & sim_summary$transition_end==2018 &
                   sim_summary$mean_events_2000==400 & sim_summary$mechanism==mech,]
  ylim <- range(c(null400$median_estimated_apc,z$median_estimated_apc)) + c(-.25,.25)
  panel <- ifelse(mech=="over-inclusion","A","B")
  ttl <- ifelse(mech=="over-inclusion","Historical over-inclusion","Improving ascertainment")
  plot(NA,xlim=c(0,50),ylim=ylim,xlab="Initial classification distortion (%)",
       ylab="Median full-period APC (%)",main=paste0(panel,"  ",ttl))
  for(j in seq_along(shape_names)){
    sh <- shape_names[j]; w <- z[z$shape==sh,]; o <- order(w$distortion_percent); w <- w[o,]
    xx <- c(0,w$distortion_percent); yy <- c(null400$median_estimated_apc,w$median_estimated_apc)
    lines(xx,yy,type="b",pch=shape_pch[j],lwd=2,col=shape_cols[j])
  }
  abline(h=5,lty=2,col="grey40")
  legend(ifelse(mech=="over-inclusion","bottomleft","topleft"),
         c("Linear","Step","Logistic"),col=shape_cols,pch=shape_pch,lwd=2,bty="n",cex=.76)
}
for(mech in c("over-inclusion","under-ascertainment")){
  z <- sim_summary[sim_summary$true_apc==5 & sim_summary$transition_end==2018 &
                   sim_summary$distortion_percent==50 & sim_summary$mechanism==mech,]
  panel <- ifelse(mech=="over-inclusion","C","D")
  ttl <- ifelse(mech=="over-inclusion","Over-inclusion","Under-ascertainment")
  plot(NA,xlim=c(1,3),ylim=c(0,100),xaxt="n",xlab="Mean events in 2000",
       ylab="Breakpoint selected (%)",main=paste0(panel,"  ",ttl," at 50%"))
  axis(1,at=1:3,labels=c(25,100,400))
  for(j in seq_along(shape_names)){
    sh <- shape_names[j]; w <- z[z$shape==sh,]; w <- w[match(c(25,100,400),w$mean_events_2000),]
    lines(1:3,w$breakpoint_selected_percent,type="b",pch=shape_pch[j],lwd=2,col=shape_cols[j])
  }
}
dev.off()

writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
cat("Completed partial redistribution, bidirectional expanded simulation, and callable calibration bundle.\n")
print(partial_summary, row.names=FALSE)
