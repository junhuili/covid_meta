
############################################################################################################
# spearman correlation for each cohort with available severity and/or inflammatory biomarkers
############################################################################################################
s <-read.table('../data/species_relative.txt', head=T, sep='	', check.names = F,row.name=1)
# remove unused metadata
s$Age<-s$Gender<-s$BMI<-s$Time<-s$Condition<-s$Severity<-NULL

r <- cor(s,method ='spearman')
r[1:10,1:10]

#Computing the p-value of correlations
cor.mtest <- function(mat, method='spearman') {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], method='spearman')
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# matrix of the p-value of the correlation
p <- cor.mtest(r)

write.table(r[1:3,-c(1:3)], file='r.txt', sep=' ', row.names=F, quote=F)
write.table(p[1:3,-c(1:3)], file='p.txt', sep=' ', row.names=F, quote=F)


##########
# heatmap
##########
rm(list=ls())

r<-read.table('r.txt', header=T, sep='\t', check.names = F,row.names=1)
p<-read.table('p.txt', header=T, sep='\t', check.names = F,row.names=1)

myp <- ifelse(p < .001, "**", ifelse(p < .01, "*", " "))

col <- c("#004174","#3a8dc2","#97c9e1","#e4f0f6","#f5f8fa","#fff9f5","#feeade","#f4ac8c","#ce4c49","#7c0e29")
breaks <- c(-0.75, -0.5, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.5, 0.75)

require(gplots)
heatmap.2(t(r), col = col, na.color = "white", breaks = breaks, Rowv=F, density.info='none', srtCol=45, trace="none", scale="none", margins=c(28,28),
          cellnote=t(myp), notecol="green",notecex=4, key.xlab="", key.ylab="", key.title="", keysize=2, lhei=c(1,12), lwid=c(1,6), cexRow=2, cexCol=2)

