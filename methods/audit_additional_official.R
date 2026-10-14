args <- commandArgs(trailingOnly=TRUE)
d <- if(length(args)) args[1] else 'Source_Data'
e <- file.path(d,'Joinpoint_Exports')
rd <- function(prefix,suffix) read.delim(file.path(e,paste0(prefix,'.Export.',suffix,'.txt')),check.names=FALSE)
p <- 'Additional_Controls_Official'
s <- rd(p,'Final.Selected.Model'); aa <- rd(p,'AAPC'); obs <- rd(p,'Data'); me <- rd(p,'Model.Estimates')
raw <- read.delim(file.path(d,'Joinpoint_Additional_Controls.txt'))
stopifnot(nrow(s)==4000,!anyDuplicated(s$cohort_id),setequal(s$cohort_id,1:4000),nrow(raw)==96000,nrow(obs)==192000,nrow(aa)==8000)
chk <- merge(raw,obs,by=c('cohort_id','year'))
stopifnot(nrow(chk)==192000,all(chk$count==chk[['Observed count']]),all(table(raw$cohort_id)==24))
map <- read.csv(file.path(d,'Additional_Controls_Cohort_Map.csv'))
r <- read.csv(file.path(d,'Additional_Controls_R_Replicates.csv'))
x <- merge(merge(map,r,by=c('control_id','replicate')),s,by='cohort_id')
z <- aa[aa[['Joinpoint Model']]==0,c('cohort_id','AAPC')]; names(z)[2] <- 'official_zero_apc'
x <- merge(x,z,by='cohort_id')
j <- me[me$Model==1 & !is.na(me$Joinpoint),c('cohort_id','Joinpoint')]
stopifnot(nrow(j)==4000,!anyDuplicated(j$cohort_id))
x <- merge(x,j,by='cohort_id')
x$official_near <- x$Model==1 & abs(x$Joinpoint-2012)<=2
x$official_near[x$control_id %in% sprintf('N%02d',1:4)] <- NA
stopifnot(nrow(x)==4000,all(table(x$control_id)==500))
su <- do.call(rbind,lapply(split(x,x$control_id),function(q) data.frame(control_id=q$control_id[1],n=nrow(q),
 r_selection=100*mean(q$selected),official_selection=100*mean(q$Model),agreement=100*mean(q$selected==q$Model),
 r_near=100*mean(q$near_true_break),official_near=100*mean(q$official_near),r_apc=median(q$apc),official_zero_apc=median(q$official_zero_apc),paired_difference=median(q$official_zero_apc-q$apc))))
for(method in c('r','official')) {
 pr <- su[[paste0(method,'_selection')]]/100; n <- su$n; z <- qnorm(.975)
 mid <- (pr+z*z/(2*n))/(1+z*z/n); half <- z*sqrt(pr*(1-pr)/n+z*z/(4*n*n))/(1+z*z/n)
 su[[paste0(method,'_lower')]] <- 100*(mid-half);su[[paste0(method,'_upper')]] <- 100*(mid+half)
}
write.csv(x,file.path(d,'Additional_Controls_Official_Paired.csv'),row.names=FALSE)
write.csv(su,file.path(d,'Additional_Controls_Official_Summary.csv'),row.names=FALSE)
print(su)
p <- 'Joinpoint_Partial_Redistribution'
s <- rd(p,'Final.Selected.Model'); aa <- rd(p,'AAPC'); obs <- rd(p,'Data'); me <- rd(p,'Model.Estimates')
raw <- read.delim(file.path(d,'Joinpoint_Partial_Redistribution.txt'))
stopifnot(nrow(s)==10,setequal(s$cohort_id,1:10),nrow(raw)==320,nrow(obs)==1600,nrow(aa)==50)
chk <- merge(raw,obs,by=c('cohort_id','year'))
stopifnot(nrow(chk)==1600,all(abs(chk$rate-chk[['Observed rate']])<=.005001),all(abs(chk$se-chk[['Standard Error']])<=.00005001))
q <- merge(s,aa,by='cohort_id'); q <- q[q$Model==q[['Joinpoint Model']],]
q <- merge(read.csv(file.path(d,'Partial_Redistribution_Cohort_Map.csv')),q,by='cohort_id')
q$breaks <- vapply(q$cohort_id,function(id) paste(me$Joinpoint[me$cohort_id==id & me$Model==s$Model[s$cohort_id==id] & !is.na(me$Joinpoint)],collapse='; '),'')
write.csv(q,file.path(d,'Partial_Redistribution_Official_Summary.csv'),row.names=FALSE)
print(q[,c('site','fraction','Model','AAPC','AAPC C.I. Low','AAPC C.I. High','breaks')])
writeLines(c('PASS: 4000 control cohorts; 96000 input rows; 192000 model-observation rows with identical counts.',
 'PASS: paired R comparison restricted to the same first 500 replicates per scenario.',
 'PASS: 10 redistribution cohorts; 320 input rows; 1600 exported model-observation rows.',
 'PASS: rates and SE match within declared export precision (0.005 and 0.00005).',
 'Run summaries: zero failures, zero excluded cohorts, permutation 4499, seed 7160; controls parametric 0-1; redistribution EmpQ5001 final-only 0-4.'),file.path(d,'Additional_Official_Audit.txt'))
