args <- commandArgs(trailingOnly=TRUE)
data_dir <- if(length(args)) args[1] else 'Source_Data'
input_path <- function(name) {
 p <- file.path(data_dir,name)
 if(!file.exists(p)) stop('Missing input: ',p)
 p
}
read_export <- function(suffix) read.delim(file.path(data_dir, paste0('Joinpoint_Calibration_All.Export.', suffix, '.txt')), check.names=FALSE)
selected <- read_export('Final.Selected.Model')
aapc <- read_export('AAPC')
models <- read_export('Model.Estimates')
bic <- read.csv(input_path('joinpoint_calibration_bic_results.csv'))
if('estimated_apc' %in% names(bic)) bic$estimated_apc.tt <- bic$estimated_apc
bic$cohort_id <- (as.integer(sub('C','',bic$calibration_id))-1L)*250L+bic$replicate
stopifnot(nrow(selected)==2000L, length(unique(selected$cohort_id))==2000L, all(sort(selected$cohort_id)==1:2000))
stopifnot(nrow(aapc)==4000L, all(table(aapc$cohort_id)==2L))
zero <- aapc[aapc[['Joinpoint Model']]==0,c('cohort_id','AAPC')]
names(zero)[2] <- 'official_zero_apc'
chosen <- merge(selected,aapc,by='cohort_id')
chosen <- chosen[chosen$Model==chosen[['Joinpoint Model']], c('cohort_id','AAPC')]
names(chosen)[2] <- 'official_selected_aapc'
result <- merge(merge(merge(bic,selected,by='cohort_id'),zero,by='cohort_id'),chosen,by='cohort_id')
stopifnot(nrow(result)==2000L)
raw <- read.delim(file.path(data_dir,'Joinpoint_Calibration_Numeric_ID.txt'))
observed <- read_export('Data')
stopifnot(nrow(raw)==48000L, nrow(observed)==96000L, !anyDuplicated(raw[c('cohort_id','year')]))
check <- merge(raw,observed,by=c('cohort_id','year'))
stopifnot(nrow(check)==96000L, all(check$count==check[['Observed count']]), all(table(raw$cohort_id)==24L))
scenario <- read.csv(input_path('joinpoint_calibration_scenarios.csv'))
summary <- do.call(rbind,lapply(split(result,result$calibration_id),function(d) data.frame(
 calibration_id=d$calibration_id[1], n=nrow(d),
 bic_break_percent=100*mean(d$break_selected==1), official_break_percent=100*mean(d$Model==1),
 agreement_percent=100*mean(d$break_selected==d$Model),
 bic_median_apc=median(d$estimated_apc.tt), official_zero_median_apc=median(d$official_zero_apc),
 paired_median_difference=median(d$official_zero_apc-d$estimated_apc.tt),
 official_selected_median_aapc=median(d$official_selected_aapc))))
summary <- merge(scenario,summary,by='calibration_id')
for(method in c('bic','official')) {
 p <- summary[[paste0(method,'_break_percent')]]/100
 summary[[paste0(method,'_break_mcse_pp')]] <- 100*sqrt(p*(1-p)/summary$n)
 z <- qnorm(.975); n <- summary$n
 mid <- (p+z*z/(2*n))/(1+z*z/n)
 half <- z*sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
 summary[[paste0(method,'_break_wilson_lower_percent')]] <- 100*(mid-half)
 summary[[paste0(method,'_break_wilson_upper_percent')]] <- 100*(mid+half)
}
write.csv(result,file.path(data_dir,'Joinpoint_Calibration_Paired_Results.csv'),row.names=FALSE)
write.csv(summary,file.path(data_dir,'Joinpoint_Calibration_Comparison_Summary.csv'),row.names=FALSE)
print(summary,row.names=FALSE)
cat('Audit passed: 2000 cohorts, 48000 input records, 96000 exported observations identical, 4000 model-level AAPCs, complete paired mapping.\n')
