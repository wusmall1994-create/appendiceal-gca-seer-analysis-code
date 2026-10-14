args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1)
out <- args[1]; dir.create(out,showWarnings=FALSE,recursive=TRUE)
years <- 2000:2023; tt <- years-2000; candidates <- c(2004,2008,2012,2016,2020)
xs <- c(list(cbind(1,tt)),lapply(candidates,function(k)cbind(1,tt,pmax(0,years-k))))
qs <- lapply(xs,function(x)qr.Q(qr(x)))
slope <- solve(crossprod(xs[[1]]),t(xs[[1]]))[2,]
screen <- function(y) {
 z <- log(y+.5)
 rss <- vapply(qs,function(q)sum((z-q%*%crossprod(q,z))^2),numeric(1))
 b <- 24*log(rss/24)+c(2,rep(4,5))*log(24)
 k <- which.min(b)-1L
 c(apc=100*expm1(sum(slope*z)),selected=as.integer(k>0),break_year=if(k>0)candidates[k] else NA)
}
# Independent reference calculation checks QR projection and parameter counting.
check <- c(448,469,472,514,535,504,586,575,569,643,662,713,715,734,830,852,860,888,943,1075,1100,1108,1160,1200)
check_bic <- vapply(seq_along(xs),function(i){m<-lm.fit(xs[[i]],log(check+.5));24*log(sum(m$residuals^2)/24)+if(i==1)2*log(24) else 4*log(24)},numeric(1))
stopifnot(screen(check)['selected']==as.integer(which.min(check_bic)>1),abs(screen(check)['apc']-100*expm1(lm.fit(xs[[1]],log(check+.5))$coefficients[2]))<1e-10)
scenarios <- data.frame(control_id=sprintf('N%02d',1:8),
 condition=rep(c('no_drift','constant_20_percent','true_break_no_drift','true_break_with_drift'),each=2),
 baseline=rep(c(400,25),4),true_break=c(rep(NA,4),rep(2012,4)),
 apc_before=c(rep(5,4),rep(2,4)),apc_after=c(rep(5,4),rep(8,4)),
 R_replicates=2000,official_replicates_prepared=500)
set.seed(20261013); rows <- list(); inputs <- list(); mapping <- list(); pos <- 1L
for(i in 1:nrow(scenarios)) {
 g <- scenarios[i,]
 latent <- g$baseline*(1+g$apc_before/100)^pmin(tt,12)*(1+g$apc_after/100)^pmax(tt-12,0)
 distortion <- if(g$condition=='constant_20_percent')rep(.2,24) else if(g$condition=='true_break_with_drift') .5*pmax(2018-years,0)/18 else rep(0,24)
 target <- 100*expm1(sum(slope*log(latent)))
 for(r in 1:2000) {
  reference <- rpois(24,latent); extra <- rpois(24,latent*distortion); observed <- reference+extra
  s <- screen(observed); s0 <- screen(reference)
  rows[[pos]] <- data.frame(control_id=g$control_id,replicate=r,target_projection_apc=target,
    apc=s['apc'],bias=s['apc']-target,selected=s['selected'],break_year=s['break_year'],
    paired_apc_difference=s['apc']-s0['apc'],reference_selected=s0['selected'],
    near_true_break=as.integer(is.finite(g$true_break)&&is.finite(s['break_year'])&&abs(s['break_year']-g$true_break)<=2))
  if(r<=500) {
   id <- (i-1)*500+r
   inputs[[id]] <- data.frame(cohort_id=id,year=years,count=observed)
   mapping[[id]] <- data.frame(cohort_id=id,control_id=g$control_id,replicate=r)
  }
  pos <- pos+1L
 }
}
d <- do.call(rbind,rows); input <- do.call(rbind,inputs)
stopifnot(nrow(input)==96000,!anyDuplicated(input[c('cohort_id','year')]),all(input$count>0),all(table(input$cohort_id)==24))
write.table(input,file.path(out,'Joinpoint_Additional_Controls.txt'),sep='\t',row.names=FALSE,quote=FALSE,eol='\r\n')
write.csv(do.call(rbind,mapping),file.path(out,'Additional_Controls_Cohort_Map.csv'),row.names=FALSE)
write.csv(scenarios,file.path(out,'Additional_Controls_Scenarios.csv'),row.names=FALSE)
write.csv(d,file.path(out,'Additional_Controls_R_Replicates.csv'),row.names=FALSE)
s <- do.call(rbind,lapply(split(d,d$control_id),function(z){
 p<-mean(z$selected);n<-nrow(z);a<-qnorm(.975);mid<-(p+a*a/(2*n))/(1+a*a/n);half<-a*sqrt(p*(1-p)/n+a*a/(4*n*n))/(1+a*a/n)
 data.frame(control_id=z$control_id[1],n=n,target_projection_apc=z$target_projection_apc[1],
 median_apc=median(z$apc),median_bias=median(z$bias),median_paired_apc_difference=median(z$paired_apc_difference),
 selected_percent=100*p,mcse_pp=100*sqrt(p*(1-p)/n),wilson_lower=100*(mid-half),wilson_upper=100*(mid+half),
 near_true_break_percent=100*mean(z$near_true_break))}))
s<-merge(scenarios,s,by='control_id')
write.csv(s,file.path(out,'Additional_Controls_R_Summary.csv'),row.names=FALSE)
print(s,row.names=FALSE)
cat('Prepared 4000 official input cohorts only; official analyses NOT executed.\n')
