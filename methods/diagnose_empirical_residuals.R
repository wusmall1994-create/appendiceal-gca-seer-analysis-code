args<-commandArgs(TRUE);stopifnot(length(args)==2)
src<-normalizePath(args[1]);out<-normalizePath(args[2]);options(stringsAsFactors=FALSE)
rd<-function(x)read.delim(file.path(src,x),check.names=FALSE)
g<-read.csv(file.path(src,'GCA_joinpoint_trends_source_data.csv'))
pm<-read.csv(file.path(src,'Partial_Redistribution_Official_Summary.csv'),check.names=FALSE)
ps<-read.csv(file.path(src,'partial_redistribution_annual_series.csv'))
gm<-rd('GCA_Official_Provenance/GCA_D1_D2_D3_Joinpoint.Export.Model.Estimates.txt')
g1<-rd('Joinpoint_GCA_From2001.Export.Model.Estimates.txt')
g1map<-read.csv(file.path(src,'GCA_From2001_Official_Summary.csv'),check.names=FALSE)
pe<-rd('Joinpoint_Pending_Exports/Joinpoint_Partial_Redistribution.Export.Model.Estimates.txt')
inputs<-list();mods<-list();ii<-0L
add<-function(id,label,z,m,group){
 ii<<-ii+1L; z<-z[order(z$year),];m<-m[order(m$Segment),]
 z$series_id<-id;z$label<-label;z$group<-group
 inputs[[ii]]<<-z[,c('series_id','label','group','year','rate','se')];mods[[ii]]<<-m
}
for(i in 1:3){z<-g[substr(g$definition,1,2)==paste0('D',i),];z$rate<-z$age_adjusted_rate_per_100k
 add(i,paste0('D',i,' 2000-2023'),z,gm[gm$Definition==paste0('D',i),],'GCA primary')}
for(i in 1:nrow(pm)){m<-pm[i,];z<-ps[ps$site==m$site & ps$fraction==m$fraction,]
 add(3+i,paste0(m$site,' ',100*m$fraction,'%'),z,pe[pe$cohort_id==m$cohort_id & pe$Model==m$Model,],if(m$fraction%in%c(0,1))'Klatskin endpoints' else 'Klatskin intermediate')}
for(i in 1:3){z<-g[substr(g$definition,1,2)==paste0('D',i)&g$year>=2001,];z$rate<-z$age_adjusted_rate_per_100k
 m<-g1map[g1map$cohort_id==i,]
 add(13+i,paste0('D',i,' 2001-2023'),z,g1[g1$cohort_id==i & g1$Model==m$Model,],'GCA start-year')}
summ<-list();res<-list();acs<-list();audit<-list()
for(i in seq_along(inputs)){
 z<-inputs[[i]];m<-mods[[i]];k<-m$Joinpoint[is.finite(m$Joinpoint)];n<-nrow(z)
 x<-cbind(1,z$year-min(z$year));if(length(k))x<-cbind(x,sapply(k,function(v)pmax(0,z$year-v)))
 w<-(z$rate/z$se)^2;fit<-lm.wfit(x,log(z$rate),w)
 slopes<-cumsum(fit$coefficients[-1]);eta<-drop(x%*%fit$coefficients)
 e<-(log(z$rate)-eta)/(z$se/z$rate)
 # Unit-marginal-scale residuals; no leverage or residual-MSE rescaling.
 ac<-as.numeric(acf(e,plot=FALSE,lag.max=4)$acf)[-1]
 dw<-sum(diff(e)^2)/sum(e^2)
 ad<-max(abs(slopes-m[['Slope Estimate']]))
 sd<-abs(sum(e^2)-m$SSE[1])
 audit[[i]]<-data.frame(series_id=z$series_id[1],label=z$label[1],max_slope_difference=ad,SSE_difference=sd)
 summ[[i]]<-data.frame(series_id=z$series_id[1],label=z$label[1],group=z$group[1],n=n,breaks=paste(k,collapse='; '),DW=dw,ACF1=ac[1],ACF2=ac[2],ACF3=ac[3],ACF4=ac[4],reference_bound=1.96/sqrt(n),SSE=sum(e^2))
 z$fitted_rate<-exp(eta);z$log_residual<-log(z$rate)-eta;z$standardized_residual<-e
 res[[i]]<-z;acs[[i]]<-data.frame(series_id=z$series_id[1],label=z$label[1],lag=1:4,ACF=ac)
}
S<-do.call(rbind,summ);R<-do.call(rbind,res);A<-do.call(rbind,audit)
write.csv(S,file.path(out,'Empirical_Residual_Diagnostics.csv'),row.names=FALSE)
write.csv(R,file.path(out,'Empirical_Residual_Values.csv'),row.names=FALSE)
write.csv(do.call(rbind,acs),file.path(out,'Empirical_Residual_ACF.csv'),row.names=FALSE)
write.csv(A,file.path(out,'Empirical_Residual_Reconstruction_Audit.csv'),row.names=FALSE)
print(S[,c('label','DW','ACF1','ACF2','reference_bound')],row.names=FALSE)
print(A,row.names=FALSE)
# Reconstruction uses full-precision rates and SE; compare exported coefficients.
stopifnot(max(A$max_slope_difference)<2e-6,max(A$SSE_difference)<.001)
plotset<-function(ids,stem,height){
 draw<-function(){par(mfrow=c(length(ids),2),mar=c(3.3,4.0,2.1,1),mgp=c(1.9,.55,0),cex=.65)
 for(id in ids){z<-R[R$series_id==id,];s<-S[S$series_id==id,]
 plot(z$year,z$standardized_residual,type='b',pch=16,col='#2B6F9C',xlab='Year',ylab='Scaled log residual',main=z$label[1],cex.main=1)
 abline(h=0,lty=2,col='grey60');if(nchar(s$breaks))abline(v=as.numeric(strsplit(s$breaks,'; ')[[1]]),lty=3,col='grey70')
 plot(1:4,unlist(s[c('ACF1','ACF2','ACF3','ACF4')]),type='h',lwd=3,col='#2B6F9C',xlim=c(.5,4.5),ylim=c(-1,1),xlab='Lag (years)',ylab='ACF',main=paste0('DW = ',sprintf('%.2f',s$DW)),cex.main=1,xaxt='n')
 axis(1,at=1:4);abline(h=0,col='grey50');abline(h=c(-1,1)*s$reference_bound,lty=2,col='#B15B38')
 }}
 tiff(file.path(out,paste0(stem,'.tiff')),width=7,height=height,units='in',res=600,compression='lzw');draw();dev.off()
 png(file.path(out,paste0(stem,'_preview.png')),width=1050,height=height*150,res=150);draw();dev.off()
}
plotset(c(1:3,14:16),'Figure_S1_GCA_Residual_Diagnostics',10)
plotset(c(4,5,12,13),'Figure_S2_Klatskin_Residual_Diagnostics',7.5)
writeLines(c('PASS: reconstructed final-model WLS slopes and weighted SSE agree with official exports within declared precision.',
 'Residual diagnostics are conditional on selected breakpoints and are descriptive, not calibrated post-selection tests.',
 'ACF bounds +/-1.96/sqrt(n) are visual references only; no independence conclusion follows from staying inside them.',
 paste('Maximum slope discrepancy:',max(A$max_slope_difference)),paste('Maximum SSE discrepancy:',max(A$SSE_difference))),file.path(out,'Empirical_Residual_Audit.txt'))
