args <- commandArgs(trailingOnly = TRUE)
fp <- function(x) {fisher.test(matrix(x,nrow=2),alternative='great')$p.value}
fr <- function(x) {fisher.test(matrix(x,nrow=2),alternative='great')$estimate}
print("reading table")
dat <- read.table(args[1],sep='\t')
print("calculating p-values")
pvals <- apply(dat[3:6],1,fp)
oddsr <- apply(dat[3:6],1,fr)
logoddsr <- log(oddsr)
pdat <- cbind(dat,pvals,oddsr,logoddsr)
colnames(pdat) <- c('cluster','feature','c+f+','c+f-','c-f+','c-f-','p.value','OddsRatio','logOR')
print("writing table")
write.table(pdat,file=args[2],sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
