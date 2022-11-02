
library(vegan)
library(BiodiversityR)
library(MASS)

#PERMANOVA relative abundance matrix
############################################

r <-read.table('species_relative.txt', head=T, sep='	', check.names = F,row.name=1)
# remove unused metadata
r$Score<-r$'D-dimer'<-r$CRP<-NULL

# extract the first collected sample from each individual
r <- r[grepl('T1', r$'Time'), ]

# Bray-curtis: relative + bray
adonis2(r[,-c(1:6)] ~ r$Severity + r$Gender + r$Age + r$BMI + r$'other factors', permutations=999, method="bray", by="margin")



#PERMANOVA CLR transfermed absolute abundance matrix
############################################

#load the CLR transfermed absolute abundance matrix
a <-read.table('species_absolute.txt', head=T, sep='	', check.names = F,row.name=1)

# extract the first collected sample from each individual
a <- a[grepl('T1', a$'Time'), ]

# Aitchison: CLR transfermed absolute abundance + euclidean
adonis2(a[,-c(1:6)] ~ a$Severity + a$Gender + a$Age + a$BMI + a$'other factors', permutations=999, method="euclidean", by="margin")


