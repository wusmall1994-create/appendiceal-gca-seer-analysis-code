args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1)
out <- normalizePath(args[1],mustWork=TRUE)
set.seed(20261014)
years <- 2000:2023; tt <- years-2000; knots <- c(2004,2008,2012,2016,2020)
xx <- c(list(cbind(1,tt)),lapply(knots,function(k)cbind(1,tt,pmax(0,years-k))))
qq <- lapply(xx,function(x)qr.Q(qr(x)))
sl <- solve(crossprod(xx[[1]]),t(xx[[1]]))[2,]
screen <- function(y){
 z<-log(y+.5); rss<-vapply(qq,function(q)sum((z-q%*%crossprod(q,z))^2),numeric(1))
 b<-24*log(rss/24)+c(2,rep(4,5))*log(24); k<-which.min(b)-1L
 c(apc=100*expm1(sum(sl*z)),selected=as.integer(k>0),break_year=if(k>0)knots[k] else NA)
}
# Common uniforms make paired mechanism comparisons possible. NB marginals have
# exact mean mu and variance mu + mu^2/50. Gaussian copula rho is NOT count rho.
sc <- expand.grid(baseline=c(25,400),mechanism=c('none','over','under'),
                  process=c('Poisson','NB_iid','NB_AR1'),stringsAsFactors=FALSE)
sc$scenario_id<-sprintf('V%02d',seq_len(nrow(sc)));sc$replicates<-2000
sc$nb_size<-ifelse(sc$process=='Poisson',NA,50)
sc$latent_rho<-ifelse(sc$process=='NB_AR1',.5,0)
res<-vector('list',nrow(sc)*2000); inp<-list(); mp<-list(); pos<-1L; id<-0L
for(i in seq_len(nrow(sc))){
 g<-sc[i,]; mu0<-g$baseline*1.05^tt
 drift<-.5*pmax(2018-years,0)/18
 mu<-mu0*if(g$mechanism=='over')1+drift else if(g$mechanism=='under')1-drift else 1
 for(r in 1:2000){
  z<-rnorm(24)
  if(g$latent_rho>0)for(t in 2:24)z[t]<-.5*z[t-1]+sqrt(.75)*z[t]
  u<-pmin(1-1e-15,pmax(1e-15,pnorm(z)))
  y<-if(g$process=='Poisson')qpois(u,mu) else qnbinom(u,size=50,mu=mu)
  y0<-if(g$process=='Poisson')qpois(u,mu0) else qnbinom(u,size=50,mu=mu0)
  s<-screen(y);s0<-screen(y0)
  res[[pos]]<-data.frame(scenario_id=g$scenario_id,replicate=r,apc=s['apc'],
      selected=s['selected'],break_year=s['break_year'],paired_apc_difference=s['apc']-s0['apc'],zeros=sum(y==0))
  pos<-pos+1L
  # Representative official batch: all nine high-volume scenarios, first 250.
  if(g$baseline==400 && r<=250){
   id<-id+1L;inp[[id]]<-data.frame(cohort_id=id,year=years,count=y)
   mp[[id]]<-data.frame(cohort_id=id,scenario_id=g$scenario_id,replicate=r)
  }
 }
}
d<-do.call(rbind,res); input<-do.call(rbind,inp)
stopifnot(nrow(d)==36000,nrow(input)==54000,all(input$count>0),!anyDuplicated(input[c('cohort_id','year')]))
wilson<-function(k,n){p<-k/n;a<-qnorm(.975);c((p+a*a/(2*n)-a*sqrt(p*(1-p)/n+a*a/(4*n*n)))/(1+a*a/n),(p+a*a/(2*n)+a*sqrt(p*(1-p)/n+a*a/(4*n*n)))/(1+a*a/n))*100}
su<-do.call(rbind,lapply(split(d,d$scenario_id),function(z){p<-mean(z$selected);ci<-wilson(sum(z$selected),nrow(z));data.frame(scenario_id=z$scenario_id[1],median_apc=median(z$apc),lo_apc=unname(quantile(z$apc,.025)),hi_apc=unname(quantile(z$apc,.975)),median_paired_difference=median(z$paired_apc_difference),selection_pct=100*p,mcse_pp=100*sqrt(p*(1-p)/nrow(z)),wilson_lo=ci[1],wilson_hi=ci[2],zero_observations=sum(z$zeros))}))
su<-merge(sc,su,by='scenario_id',sort=TRUE)
write.csv(sc,file.path(out,'Robustness_Scenarios.csv'),row.names=FALSE)
write.csv(d,file.path(out,'Robustness_R_Replicates.csv'),row.names=FALSE)
write.csv(su,file.path(out,'Robustness_R_Summary.csv'),row.names=FALSE)
write.csv(do.call(rbind,mp),file.path(out,'Robustness_Cohort_Map.csv'),row.names=FALSE)
write.table(input,file.path(out,'Joinpoint_Distribution_Robustness.txt'),sep='\t',row.names=FALSE,quote=FALSE,eol='\r\n')
# Start-year check: a no-break WLS summary, not official segmented AAPC.
gca<-read.csv(file.path(out,'GCA_joinpoint_trends_source_data.csv'))
sr<-list();gi<-list();gm<-list();id<-0L
for(def in unique(gca$Definition))for(start in c(2000,2001)){
 a<-gca[gca$Definition==def & gca$year>=start,];a<-a[order(a$year),]
 stopifnot(all(diff(a$year)==1),all(a$se>0),all(a$age_adjusted_rate_per_100k>0))
 m<-lm(log(age_adjusted_rate_per_100k)~year,data=a,weights=(age_adjusted_rate_per_100k/se)^2)
 b<-coef(m)[2];ci<-confint(m)[2,];id<-id+1L
 sr[[id]]<-data.frame(cohort_id=id,definition=def,start=start,end=2023,n=nrow(a),apc=100*expm1(b),lower=100*expm1(ci[1]),upper=100*expm1(ci[2]))
 gi[[id]]<-data.frame(cohort_id=id,year=a$year,rate=a$age_adjusted_rate_per_100k,se=a$se)
 gm[[id]]<-data.frame(cohort_id=id,definition=def,start=start,end=2023)
}
write.csv(do.call(rbind,sr),file.path(out,'GCA_Start_Year_WLS_Summary.csv'),row.names=FALSE)
write.csv(do.call(rbind,gm),file.path(out,'GCA_Start_Year_Cohort_Map.csv'),row.names=FALSE)
write.table(do.call(rbind,gi),file.path(out,'Joinpoint_GCA_Start_Year.txt'),sep='\t',row.names=FALSE,quote=FALSE,eol='\r\n')
from2001<-do.call(rbind,gi)[do.call(rbind,gi)$cohort_id %in% c(2,4,6),]
from2001$cohort_id<-match(from2001$cohort_id,c(2,4,6))
write.table(from2001,file.path(out,'Joinpoint_GCA_From2001.txt'),sep='\t',row.names=FALSE,quote=FALSE,eol='\r\n')
print(su,row.names=FALSE);print(do.call(rbind,sr),row.names=FALSE)
cat('Joinpoint input files written.\n')
