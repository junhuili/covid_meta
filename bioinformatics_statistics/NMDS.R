
####################################
#nonmetric multidimensional scaling
####################################

# load the required packages in the analysis
library (vegan)
library (ape)
library (ggplot2)

# load in the taxonomic data with first column as row names
m<-read.csv('merge_species_relative_meta.txt', header=T, sep='	', stringsAsFactors = FALSE, check.names = F,row.names = 1)
#m<-read.csv('merge_species_absolute_meta.txt', header=T, sep='	', stringsAsFactors = FALSE, check.names = F,row.names = 1)

m[1:5,1:9]
m$UNCLASSIFIED<-NULL
dim(m)

# load in the taxonomic data
n <-read.delim2('merge_species_relative_meta.txt', head=T, sep='	', check.names = F, stringsAsFactors = FALSE)
#n<-read.delim2('merge_species_absolute_meta.txt', header=T, sep='	', stringsAsFactors = FALSE, check.names = F)

n$UNCLASSIFIED<-NULL
n[1:4,1:10]

# subtset metadata
m0<-n[,c(1:9)]
head(m0)
dim(m0) #

# subset the tax data
m[1:5,1:9]
m1<-m[,-c(1:8)]
dim(m1) #
m1[1:5,1:6]


# Bray-curtis distance using relative abundance
set.seed(123)
NMDS<- metaMDS(m1, dist = "bray")

# Aitchison distance using clr-transformed absolute abudance
#set.seed(123)
#NMDS<- metaMDS(m1, dist = "euclidean")


m3<-scores(NMDS,display=c("sites"))
head(m3)

# combining NMDS score and metadata
m4<-cbind(m3,m0)
m4$Severity <- factor(m4$Severity, levels=c("Control","Mild","Moderate","Severe","Critical","Fatal","COVID"))


require(cowplot)

p<-ggplot(m4, aes(NMDS1,NMDS2,colour=Severity))

p<-p+geom_point(data=m4,aes(NMDS1,NMDS2,fill=Severity,colour=Severity),size=1.5,alpha=0.5,stroke=0.8,shape=21)
p<-p + scale_fill_manual(values=c("Control"="#008080","Asymptomatic"="#87CEEB","Mild"="blue","Moderate"="orange","Severe"="#E34234","Critical"="#953553","Fatal"="black","COVID"="yellow","non-ICU"="#8a75ba","ICU"="#e334b7"))
p<-p + scale_color_manual(values=c("Control"="#008080","Asymptomatic"="#87CEEB","Mild"="blue","Moderate"="orange","Severe"="#E34234","Critical"="#953553","Fatal"="black","COVID"="yellow","non-ICU"="#8a75ba","ICU"="#e334b7"))
p<-p+theme_bw()
p<-p+theme(panel.background = element_rect(fill = NULL),legend.position="right",
           panel.grid.major = element_line(colour = NULL))
p<-p + theme(axis.text.x=element_text(colour="black",angle=0,hjust=0.5,vjust=0.5,face="plain",size=9),
             axis.text.y=element_text(colour="black",size=10),
             axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=0),
             axis.title.y = element_text(colour="black",size=14)) #,face="bold"
p<-p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=1),strip.text.x = element_text(size=12,angle=0,face="bold"), strip.text.y=element_text(size=16,face="bold"),strip.background=element_rect(colour=NA, fill = NA))
p<-p+theme(legend.position="NULL")
p


# without unclassified COVID
#################################################################################

exclude<-'COVID'
m4 <- m4[!grepl(exclude, m4$'Severity'), ]
dim(m4)

m4$Severity <- factor(m4$Severity, levels=c("Control","Mild","Moderate","Severe","Critical","Fatal"))

require(cowplot)
tail(m4)
p<-ggplot(m4, aes(NMDS1,NMDS2,colour=Severity))#

p<-p+geom_point(data=m4,aes(NMDS1,NMDS2,fill=Severity,colour=Severity),size=1.5,alpha=0.5,stroke=0.8,shape=21)
p<-p + scale_fill_manual(values=c("Control"="#008080","Asymptomatic"="#87CEEB","Mild"="blue","Moderate"="orange","Severe"="#E34234","Critical"="#953553","Fatal"="black","COVID"="yellow","non-ICU"="#8a75ba","ICU"="#e334b7"))
p<-p + scale_color_manual(values=c("Control"="#008080","Asymptomatic"="#87CEEB","Mild"="blue","Moderate"="orange","Severe"="#E34234","Critical"="#953553","Fatal"="black","COVID"="yellow","non-ICU"="#8a75ba","ICU"="#e334b7"))

require("ggrepel")
set.seed(42)
p<-p+theme_bw()
p<-p+theme(panel.background = element_rect(fill = NULL),legend.position="right",
           panel.grid.major = element_line(colour = NULL))

p<-p + theme(axis.text.x=element_text(colour="black",angle=0,hjust=0.5,vjust=0.5,face="plain",size=9),
             axis.text.y=element_text(colour="black",size=10),
             axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=0),
             axis.title.y = element_text(colour="black",size=14)) #,face="bold"
p<-p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=1),strip.text.x = element_text(size=12,angle=0,face="bold"), strip.text.y=element_text(size=16,face="bold"),strip.background=element_rect(colour=NA, fill = NA))
p<-p+theme(legend.position="NULL")
p<-p+theme(legend.position="right")
p



#nmds_bray-curtis_log2_day
#################################################################################

library (vegan)
library (ape)
library (ggplot2)

m5 <- m4[!grepl('Control', m4$'Day'), ]
m5$`Day`<-as.numeric(m5$`Day`)
summary(m5$Day)
m5$`log2(Day)`<-log2(m5$`Day`+1)

m4$Condition <- factor(m4$Condition, levels=c("Control","COVID"))

centroids <- aggregate(cbind(NMDS1,NMDS2) ~ Condition, m4, mean)
segs <- merge(m4, setNames(centroids, c('Condition','oNMDS1','oNMDS2')), by = 'Condition', sort = FALSE)


p<-ggplot(m5, aes(NMDS1,NMDS2))

p<-p+geom_segment(data=segs, mapping=aes(xend=oNMDS1, yend=oNMDS2,alpha=Condition),colour="#008080",size=0.2)

p<-p+geom_point(data=m4,aes(NMDS1,NMDS2,alpha=Condition),colour="#008080",fill="#008080",size=1.5,stroke=0.3,shape=21)
p<-p + scale_alpha_manual(values=c("Control"=0.9,"COVID"=0))
p<-p+geom_point(aes(fill=`log2(Day)`),size=1.5,alpha=0.9,stroke=0.3,shape=21,colour='black')
p<-p+scale_color_gradient(low="#feeade", high="#ff6000")
p<-p+scale_fill_gradient(low="#feeade", high="#ff6000")
p<-p+geom_density_2d(aes(color = ..level..),alpha=0.9,size=0.5)

p<-p+geom_point(data=centroids,aes(alpha=Condition),size=3,fill="#008080",colour="#008080",shape=21)
p<-p+theme_bw()
p<-p+theme(panel.background = element_rect(fill = NULL),legend.position="right",
           panel.grid.major = element_line(colour = NULL))

p<-p + theme(axis.text.x=element_text(colour="black",angle=0,hjust=0.5,vjust=0.5,face="plain",size=9),
             axis.text.y=element_text(colour="black",size=10),
             axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=0),
             axis.title.y = element_text(colour="black",size=14)) #,face="bold"
p<-p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=1),strip.text.x = element_text(size=12,angle=0,face="bold"), strip.text.y=element_text(size=16,face="bold"),strip.background=element_rect(colour=NA, fill = NA))
p<-p+theme(legend.position="NULL")
p


#################################################################################
#longtitudinal analyses of PRJNA714459 and PRJNA703303
#################################################################################
# load the required packages in the analysis
library (vegan)
library (ape)
library (ggplot2)

# load in the taxonomic data with first column as row names
m <-read.table('species_relative_meta.txt', head=T, sep='	', check.names = F,row.name=1)
m[1:5,1:7]
m$UNCLASSIFIED<-NULL

# load in the taxonomic data
n <-read.delim2('species_relative_meta.txt', head=T, sep='	', check.names = F, stringsAsFactors = FALSE)
n$UNCLASSIFIED<-NULL
n[1:4,1:8]

# subset metadata
m0<-n[,c(1:7)]
head(m0)

# subset the tax data
m[1:5,1:7]
m1<-m[,-c(1:6)]
m1[1:5,1:6]

set.seed(123)
NMDS<- metaMDS(m1, dist = "bray")

m3<-scores(NMDS,display=c("sites"))
head(m3)

# combining NMDS score and metadata
m4<-cbind(m3,m0)
head(m4)
dim(m4)

# Phase
unique(m$Phase) #Control M0 M1 M6 M9
m4$Phase <- factor(m4$Phase, levels=c("Control","M0","M1","M6","M9"))
centroids <- aggregate(cbind(NMDS1,NMDS2) ~ Phase, m4, mean)
segs <- merge(m4, setNames(centroids, c('Phase','oNMDS1','oNMDS2')),
              by = 'Phase', sort = FALSE)

p<-ggplot(m4, aes(NMDS1,NMDS2,colour=Phase))#
p<-p+stat_ellipse(data=m4,aes(NMDS1,NMDS2,colour=Phase),alpha=0.4,type='t',size=0.8, linetype = 1) 
p<-p+geom_point(data=m4,aes(NMDS1,NMDS2,fill=Phase,colour=Phase),size=1.5,alpha=0.7,stroke=0.5) 
p<-p+geom_segment(data=segs, mapping=aes(xend=oNMDS1, yend=oNMDS2, colour=Phase),alpha=0.7,size=0.3)
p<-p+geom_point(data=centroids,size=3,aes(colour=Phase),alpha=0.9)
p<-p + scale_fill_manual(values=c("Control"="#008080","M0"="#953553","M1"="#E34234","M6"="orange","M9"="#87CEEB"))
p<-p + scale_color_manual(values=c("Control"="#008080","M0"="#953553","M1"="#E34234","M6"="orange","M9"="#87CEEB"))
p<-p + scale_shape_manual(values=c(21,24,0))
p<-p+theme_bw()
p<-p+theme(panel.background = element_rect(fill = NULL),legend.position="right",panel.grid.major = element_line(colour = NULL))
p<-p + theme(axis.text.x=element_text(colour="black",angle=0,hjust=0.5,vjust=0.5,face="plain",size=9),
             axis.text.y=element_text(colour="black",size=10),
             axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=0),
             axis.title.y = element_text(colour="black",size=14))
p<-p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=1),strip.text.x = element_text(size=12,angle=0,face="bold"), strip.text.y=element_text(size=16,face="bold"),strip.background=element_rect(colour=NA, fill = NA))
p<-p+theme(legend.position="NULL")
p



######################################################################################################
#PERMANOVA between COVID samples at each month versus healthy controls (PRJNA714459 and PRJNA703303)
######################################################################################################
library(vegan)
library(BiodiversityR)
library(MASS)

#######################
#Control M0 M1 M6 M9
#######################
c <- m[grepl('Control', m$'Phase'), ]
c[1:2,1:5]
dim(c)

c1 <- m[grepl('M0', m$'Phase'), ]
c1[1:2,1:5]

c2 <- m[grepl('M1', m$'Phase'), ]
c2[1:2,1:5]
c2$Phase

c3 <- m[grepl('M6', m$'Phase'), ]
c3[1:2,1:5]
c3$Phase

c4 <- m[grepl('M9', m$'Phase'), ]
c4[1:2,1:5]
c4$Phase

d1<-rbind(c1,c)
set.seed(36) 
adonis2(d1[,-c(1:6)] ~ d1$Phase, permutations=999, method="bray")

d2<-rbind(c2,c)
set.seed(36) 
adonis2(d2[,-c(1:6)] ~ d2$Phase, permutations=999, method="bray")

d3<-rbind(c3,c)
set.seed(36) 
adonis2(d3[,-c(1:6)] ~ d3$Phase, permutations=999, method="bray")

d4<-rbind(c4,c)
set.seed(36)
adonis2(d4[,-c(1:6)] ~ d4$Phase, permutations=999, method="bray")

#################################################################################


