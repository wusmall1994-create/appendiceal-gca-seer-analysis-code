args <- commandArgs(trailingOnly=TRUE)
d <- if(length(args)) args[1] else 'Source_Data'
rd <- function(p,s) read.delim(file.path(d,paste0(p,'.Export.',s,'.txt')),check.names=FALSE)
wilson <- function(k,n) {z<-qnorm(.975);p<-k/n;c(100*(p+z*z/(2*n)-z*sqrt(p*(1-p)/n+z*z/(4*n*n)))/(1+z*z/n),100*(p+z*z/(2*n)+z*sqrt(p*(1-p)/n+z*z/(4*n*n)))/(1+z*z/n))}
p <- 'Joinpoint_Distribution_Robustness'
s<-rd(p,'Final.Selected.Model'); aa<-rd(p,'AAPC'); obs<-rd(p,'Data')
raw<-read.delim(file.path(d,paste0(p,'.txt')))
stopifnot(nrow(s)==2250,!anyDuplicated(s$cohort_id),setequal(s$cohort_id,1:2250),nrow(raw)==54000,nrow(obs)==108000,nrow(aa)==4500)
ch<-merge(raw,obs,by=c('cohort_id','year'));stopifnot(nrow(ch)==108000,all(ch$count==ch[['Observed count']]),all(table(raw$cohort_id)==24))
map<-read.csv(file.path(d,'Robustness_Cohort_Map.csv'));r<-read.csv(file.path(d,'Robustness_R_Replicates.csv'))
x<-merge(merge(map,r,by=c('scenario_id','replicate')),s,by='cohort_id')
z<-aa[aa[['Joinpoint Model']]==0,c('cohort_id','AAPC')];names(z)[2]<-'official_zero_apc';x<-merge(x,z,by='cohort_id')
z<-merge(s,aa,by='cohort_id');z<-z[z$Model==z[['Joinpoint Model']],c('cohort_id','AAPC')];names(z)[2]<-'official_selected_aapc';x<-merge(x,z,by='cohort_id')
stopifnot(nrow(x)==2250,all(table(x$scenario_id)==250),all(x$replicate<=250))
su<-do.call(rbind,lapply(split(x,x$scenario_id),function(q){ciR<-wilson(sum(q$selected),nrow(q));ciO<-wilson(sum(q$Model),nrow(q));data.frame(scenario_id=q$scenario_id[1],n=nrow(q),r_selection=100*mean(q$selected),r_low=ciR[1],r_high=ciR[2],official_selection=100*mean(q$Model),official_low=ciO[1],official_high=ciO[2],agreement=100*mean(q$selected==q$Model),r_apc=median(q$apc),official_zero_apc=median(q$official_zero_apc),paired_difference=median(q$official_zero_apc-q$apc),official_selected_aapc=median(q$official_selected_aapc))}))
write.csv(x,file.path(d,'Robustness_Official_Paired.csv'),row.names=FALSE);write.csv(su,file.path(d,'Robustness_Official_Summary.csv'),row.names=FALSE);print(su)
p<-'Joinpoint_GCA_From2001';s<-rd(p,'Final.Selected.Model');aa<-rd(p,'AAPC');obs<-rd(p,'Data');me<-rd(p,'Model.Estimates');ap<-rd(p,'APC')
raw<-read.delim(file.path(d,paste0(p,'.txt')))
stopifnot(nrow(s)==3,setequal(s$cohort_id,1:3),nrow(raw)==69,nrow(obs)==345,nrow(aa)==15)
ch<-merge(raw,obs,by=c('cohort_id','year'));stopifnot(nrow(ch)==345,all(abs(ch$rate-ch[['Observed rate']])<=.005001),all(abs(ch$se-ch[['Standard Error']])<=.00005001))
q<-merge(s,aa,by='cohort_id');q<-q[q$Model==q[['Joinpoint Model']],];q$definition<-paste0('D',q$cohort_id)
q$breaks<-vapply(q$cohort_id,function(id) paste(me$Joinpoint[me$cohort_id==id & me$Model==s$Model[s$cohort_id==id] & !is.na(me$Joinpoint)],collapse='; '),'')
write.csv(q,file.path(d,'GCA_From2001_Official_Summary.csv'),row.names=FALSE)
ap<-merge(s,ap,by=c('cohort_id','Model'));stopifnot(nrow(ap)==6);write.csv(ap,file.path(d,'GCA_From2001_Official_APC.csv'),row.names=FALSE)
print(q);print(ap)
for(prefix in c('Joinpoint_GCA_From2001','Joinpoint_Distribution_Robustness')) {
 rs<-readLines(file.path(d,paste0(prefix,'.Export.Run.Summary.txt')),warn=FALSE,skipNul=TRUE)
 stopifnot(any(grepl("Couldn't Be Executed\\s+0",rs)),any(grepl('Number of Excluded Cohorts\\s+0',rs)),any(grepl('Seed for Randomly Permuting Data\\s+7160',rs)))
}
writeLines(c('PASS: 2250 cohorts, 54000 input counts, 108000 exported observations; counts match exactly.',
'PASS: same first 250 R replicates per scenario; nine scenarios, no failures/exclusions.',
'PASS: three GCA cohorts, 69 input rows, 345 exported rows; rate/SE match export precision (0.005/0.00005).',
'PASS: permutation seed 7160; run summaries retained. GCA 0-4 EmpQ5001 final-only; robustness 0-1 parametric.',
'NOTE: R and official no-break slopes are distinct from selected-model AAPCs. Official comparisons use baseline 400.'),file.path(d,'Robustness_Official_Audit.txt'))
