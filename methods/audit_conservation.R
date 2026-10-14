args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1)
root <- args[1]
d <- read.csv(file.path(root,'KT_all_definitions_clean_1992_2023.csv'),check.names=FALSE,fileEncoding='UTF-8-BOM')
s <- function(n) {z<-d[d$Definition==n,];z[order(z$Year),]}
ka<-s('K-all');ki<-s('K-intra');ke<-s('K-extra');i0<-s('ICC-observed');i1<-s('ICC-corrected');e0<-s('ECC-observed');e1<-s('ECC-corrected')
stopifnot(nrow(d)==224,all(vapply(split(d$Year,d$Definition),function(y)identical(sort(y),1992:2023),logical(1))),all(vapply(split(d$Population,d$Year),function(y)length(unique(y))==1,logical(1))))
rows<-list();ix<-1
for(v in c('Count','Rate')) {
 errors<-cbind(ka[[v]]-ki[[v]]-ke[[v]],i0[[v]]-i1[[v]]-ki[[v]],e1[[v]]-e0[[v]]-ki[[v]],i0[[v]]+e0[[v]]-i1[[v]]-e1[[v]])
 for(j in 1:4){rows[[ix]]<-data.frame(year=ka$Year,measure=v,identity=c('Kall=Kintra+Kextra','ICC0=ICC1+Kintra','ECC1=ECC0+Kintra','ICC0+ECC0=ICC1+ECC1')[j],error=errors[,j],tolerance=if(v=='Count')0 else 4e-6);ix<-ix+1}
}
a<-do.call(rbind,rows);a$pass<-abs(a$error)<=a$tolerance
stopifnot(all(a$pass))
write.csv(a,file.path(root,'KT_Annual_Conservation_Audit.csv'),row.names=FALSE)
cat('PASS: 128 annual count checks and 128 annual rate checks; common denominators; all seven series complete.\n')
