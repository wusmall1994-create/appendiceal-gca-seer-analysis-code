out <- commandArgs(trailingOnly=TRUE)[1]
x<-read.delim(file.path(out,'Joinpoint_GCA_Start_Year.txt'))
a<-x[x$cohort_id %in% c(2,4,6),];a$cohort_id<-match(a$cohort_id,c(2,4,6))
stopifnot(nrow(a)==69,all(table(a$cohort_id)==23),all(a$year>=2001),!anyDuplicated(a[c('cohort_id','year')]))
write.table(a,file.path(out,'Joinpoint_GCA_From2001.txt'),sep='\t',row.names=FALSE,quote=FALSE,eol='\r\n')
y<-read.delim(file.path(out,'Joinpoint_Distribution_Robustness.txt'))
map<-read.csv(file.path(out,'Robustness_Cohort_Map.csv'))
r<-read.csv(file.path(out,'Robustness_R_Replicates.csv'))
stopifnot(nrow(y)==54000,length(unique(y$cohort_id))==2250,all(table(y$cohort_id)==24),all(y$count>0),all(table(map$scenario_id)==250),nrow(r)==36000)
# Independent lm.fit calculation on 30 deterministically selected official-input series.
errors<-c()
for(id in unique(round(seq(1,2250,length.out=30)))){
 b<-y[y$cohort_id==id,];m<-map[map$cohort_id==id,];z<-log(b$count+.5);tt<-b$year-2000
 X<-cbind(1,tt);fit<-lm.fit(X,z)
 apc<-100*expm1(fit$coefficients[2])
 bic<-c(24*log(sum(fit$residuals^2)/24)+2*log(24),sapply(c(2004,2008,2012,2016,2020),function(k){f<-lm.fit(cbind(X,pmax(0,b$year-k)),z);24*log(sum(f$residuals^2)/24)+4*log(24)}))
 target<-r[r$scenario_id==m$scenario_id & r$replicate==m$replicate,]
 stopifnot(nrow(target)==1,target$selected==as.integer(which.min(bic)>1))
 errors<-c(errors,abs(target$apc-apc))
}
stopifnot(max(errors)<1e-10)
writeLines(c('PASS: 36000 R replicates; 18 scenarios x 2000.',
 'PASS: official distribution inputs 2250 cohorts x 24 years, positive counts, 250 per selected scenario.',
 'PASS: GCA from-2001 inputs 3 cohorts x 23 years, with separate annual rate and SE fields.',
 paste('PASS: independent lm.fit / BIC check on 30 series, max APC error',format(max(errors),scientific=TRUE)),
 'NOT EXECUTED: new official Joinpoint distribution and start-year analyses.'),file.path(out,'Robustness_Input_QA.txt'))
cat(readLines(file.path(out,'Robustness_Input_QA.txt')),sep='\n')
