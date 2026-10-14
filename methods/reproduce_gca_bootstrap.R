# Reproduce GCA Table 1 and Figure 2 incidence contrasts from aggregate cells.
# Preserve component-level Poisson draws shared across D1, D2 and D3.
args <- commandArgs(TRUE); stopifnot(length(args)==2)
inc_file<-args[1]; out_dir<-args[2];dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
options(stringsAsFactors=FALSE,scipen=999)
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
if (!all(c('histology','year','age','cases','population') %in% names(inc))) {
 inc$histology <- extract_hist(inc[["ICD-O-3 Hist/behav, malignant"]])
 inc$year <- suppressWarnings(as.integer(inc[["Year of diagnosis"]]))
 inc$age <- age_key(inc[["Age recode with <1 year olds and 90+"]])
 inc$cases <- suppressWarnings(as.numeric(inc[["Count"]]))
 inc$population <- suppressWarnings(as.numeric(inc[["Population"]]))
}
inc <- inc[inc$histology %in% 8243:8245 & inc$year %in% 2000:2023 & inc$age %in% names(std_w), ]
stopifnot(nrow(inc)==1440L,all(is.finite(inc$cases)),all(inc$population>0))
write_csv(inc[,c('histology','year','age','cases','population')],'GCA_Bootstrap_Annual_Cells.csv')

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
