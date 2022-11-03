#######################
# dysbiosis score
#######################
library(vegan)
library(BiodiversityR)
library(MASS)
library(ggplot2)
library(ggbreak)

m<-read.csv('../data/merge_species_rel_meta.txt', header=T, sep='	', stringsAsFactors = FALSE, check.names = F,row.names = 1)
m$UNCLASSIFIED<-NULL
dim(m)

m1 <- as.matrix(as.data.frame(lapply(m[,-c(1:8)], as.numeric)))
m1[1:5,1:5]

# bray-curtis distances between every two samples
d <- vegdist(m1, method = "bray")

d1 <- data.frame(t(combn(rownames(m),2)), as.numeric(d))
names(d1) <- c("Row.names", "c2", "bray.rel")
d1[1:5,]
"   Row.names         c2  bray.rel
1 ERR5445742 ERR5445743 0.9303872
2 ERR5445742 ERR5445744 0.9196088
3 ERR5445742 ERR5445745 0.8352233
4 ERR5445742 ERR5445746 0.9314631
5 ERR5445742 ERR5445747 0.7911285"

n<-read.csv('up/merge_species_rel_meta.txt', header=T, sep='	', stringsAsFactors = FALSE, check.names = F)
n[1:2,1:4]
"   Row.names Subject Condition Severity
1 ERR5445742    M001     COVID     Mild
2 ERR5445743    M002     COVID     Mild"

n1<-subset(n,select=c('Row.names','Severity'))
head(n1)
dim(n1)

# to extract distance matrix of each COVID sample to each healthy control
##########################################################################
n2 <- n1[grepl('Control', n1$'Severity'), ]
dim(n2) 
head(n2)

d1[1:5,]
d2a<-merge(n2,d1,by=c('Row.names'),all.x=T)
head(d2a)
dim(d2a)

head(n2)
d2b<-merge(n2,d1,by.x=c('Row.names'),by.y=c('c2'),all.x=T)
head(d2b)
dim(d2b)
d2b<-d2b[complete.cases(d2b[ ,3]),]
dim(d2b)
names(d2b)<-c('Row.names','Severity','c2','bray.rel')

d3<-rbind(d2a,d2b)
dim(d3)
write.table(d3, file='bray.rel.txt', sep='	', row.names=F, quote=F)
##########################################################################


# median Bray-Curtis distance of a COVID microbiome to that of every healthy control
##########################################################################
head(d3)
d4<-aggregate(bray.rel ~ c2, d3, median)
head(d4)
dim(d4)
names(d4)<-c('Row.names','bray.rel')

df<-merge(n[,c(1:9)], d4, by=c('Row.names'),all=T)
head(df)
dim(df)
write.table(df, file='bray.rel.median.txt', sep='	', row.names=F, quote=F)

library(plyr)
mp <- ddply(df, "Severity", summarise, percent90=quantile(bray.rel, probs = c(0.9)))
mp <- mp[grepl('Control', mp$'Severity'), ]
percent90 <- mp[[2]]

p<-ggplot(df, aes(x=bray.rel,color=Severity,fill=Severity)) #
p<-p+ geom_density(alpha=0.1)
p<-p+ theme_classic() 
p<-p+ geom_vline(data=percent90, color="black", linetype="dashed") 
p<-p + theme_classic() + labs(x = "Dysbiosis score") + labs(y = "Density")
p<-p + scale_fill_manual(values=c("Control"="#008080","Asymptomatic"="#87CEEB","Mild"="blue","Moderate"="orange","Severe"="#E34234","Critical"="#953553","Fatal"="black","COVID"="yellow"))
p<-p + scale_color_manual(values=c("Control"="#008080","Asymptomatic"="#87CEEB","Mild"="blue","Moderate"="orange","Severe"="#E34234","Critical"="#953553","Fatal"="black","COVID"="yellow"))
p<-p+theme(legend.position="NULL")
p <- p + scale_x_continuous(limits = c(0.6,1), expand = c(0.005,0.005), breaks=c(0.6,0.7,0.8,0.9,1))
p <- p + scale_y_continuous(limits = c(0,5), expand = c(0,0), breaks=c(0,1,2,3,4,5))
p <- p + annotate("text", x=0.905, y=4, label= "Dysbiosis")
p <- p + theme(text = element_text(size=15), axis.text.x = element_text(angle=0, hjust=0.5))
p




#######################
# dysbiotic proportion
#######################
library(plyr)
mu <- ddply(df, "Severity", summarise, grp.median=median(bray.rel))
head(mu)

mp <- ddply(df, "Severity", summarise, percent90=quantile(bray.rel, probs = c(0.9)))
mp <- mp[grepl('Control', mp$'Severity'), ]
head(mp)

df$hit<-1

df1 <-df[ which(df$bray.rel > mp[1,2]),]
df2<-aggregate(hit~Severity,df1,sum)
names(df2)<-c('Severity','hit90')
df3<-aggregate(hit~Severity,df,sum)
df4<-merge(df2,df3)
df4$percent<-df4$hit90/df4$hit
write.table(df4, file='bray.rel.90p.count.txt', sep='	', row.names=F, quote=F)


# manually curated "dysbiosis.txt" based on "bray.rel.90p.count.txt" for figure generation
r <-read.table('dysbiosis.txt', head=T, sep="\t", check.names = F,stringsAsFactors = FALSE) #
r$Seq <- factor(r$Seq, levels=c("shotgun metagenome","16S amplicon"))
r$Severity <- factor(r$Severity, levels=c("Control","Asymptomatic","Mild","Moderate","Severe","Critical","Fatal"))

library(ggplot2)
p <- ggplot(r, aes(Score, as.numeric(as.character(Dysbiosis)) ))
p <- p + geom_point(alpha=1,size=3,aes(shape=Seq,fill=NULL,color=Severity), stroke=1.5)
p <- p + geom_point(alpha=0.2,size=3,aes(shape=Seq,fill=Severity,color=Severity), stroke=1.5)
p <- p + scale_shape_manual(values=c(21,24))
p <- p + ggpubr::stat_cor(method = "spearman", aes(colour=Seq), label.x = 1, color='black', size=3) 
p <- p + geom_smooth(method="lm", formula = y ~ x, aes(colour=Seq,linetype=Seq), se=T, size=1, alpha=0.2)
p <- p + labs(y = "Dysbiosis frequency (%)") + labs(x = "Severity score") + theme_classic() + theme(legend.justification = "bottom")
p <- p + theme(text = element_text(size=13), axis.text.x = element_text(angle=0, hjust=0.5))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black"))
p <- p + expand_limits(x=c(-0.3, 10.3))
p <- p + scale_x_continuous(expand = c(0,0), breaks=c(0,1,2,4,6,8,10))
p <- p + expand_limits(y=c(0, 100))
p <- p + scale_y_continuous(expand = c(0.01,0.01), breaks=c(0,20,40,60,80,100))
p <- p + theme(legend.position="bottom", legend.box = "horizontal") + guides(colour = guide_legend(nrow = 2))
p<-p + scale_fill_manual(values=c("shotgun metagenome"="black","16S amplicon"="black","Control"="#008080","Asymptomatic"="#87CEEB","Mild"="blue","Moderate"="orange","Severe"="#E34234","Critical"="#953553","Fatal"="black"))
p<-p + scale_color_manual(values=c("shotgun metagenome"="black","16S amplicon"="black","Mild"="blue","Control"="#008080","Asymptomatic"="#87CEEB","Moderate"="orange","Severe"="#E34234","Critical"="#953553","Fatal"="black"))
p <- p + theme(legend.position="NULL")
p <- p + annotate("text", x=6.5, y=105, label= "Shotgun metagenome",size=4)
p <- p + annotate("text", x=4, y=12, label= "16S amplicon",size=4)
p


r1 <- r[which(r$Seq == 'shotgun metagenome'),]
fit1 <- lm(Dysbiosis~Score, data=r1)
summary(fit1)

r2 <- r[which(r$Seq == '16S amplicon'),]
fit2 <- lm(Dysbiosis~Score, data=r2)
summary(fit2)


library("ggpubr")
res1<-cor.test(r1$Dysbiosis,r1$Score, method ='spearman')
res1$p.value
res1$estimate

res2<-cor.test(r2$Dysbiosis,r2$Score, method ='spearman')
res2$p.value
res2$estimate
#########################################################################################################
#########################################################################################################



################################################################################
# Within-subject Bray-Curtis distance of gut microbiome in COVID-19 individuals
################################################################################
a <-read.delim2('../data/shift.txt', head=T, sep='	', check.names = F, stringsAsFactors = FALSE)
a1<-a[,c(2,10)]
dim(unique(a1)) 

a$Severity <- factor(a$Severity, levels=c("Mild","Moderate","Severe","Critical","COVID"))
a$Day<-as.numeric(as.character(a$Day))+1
a$bray.rel<-as.numeric(as.character(a$bray.rel))

library(ggplot2)
p <- ggplot(a, aes(as.numeric(as.character(Day)), bray.rel))
p <- p + geom_point(alpha=0.5,size=3,shape=21,aes(fill=Severity),color='black')
p <- p + geom_smooth(method="lm", formula = y ~ x, aes(colour=Severity,fill=Severity), se=F, size=1, alpha=0.2)
p <- p + geom_smooth(method="lm", formula = y ~ x,  colour='green', se=T, fill='green', size=1.5, alpha=0.1)
p <- p + labs(y = "Within-subject Bray-Curtis Disssimilarity") + labs(x = "Day") + theme_classic() + theme(legend.justification = "bottom")
p <- p + theme(text = element_text(size=13), axis.text.x = element_text(angle=0, hjust=0.5))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black"))
p <- p + expand_limits(y=c(0, 1.02))
p <- p + scale_y_continuous(expand = c(0.01,0), breaks=c(0,0.2,0.4,0.6,0.8,1))
p <- p + expand_limits(x=c(-10, 280))
p <- p + scale_x_continuous(expand = c(0.05,0.05), breaks=c(0,7,14,30,60,180,270))
p <- p + scale_fill_discrete(name="Day", breaks=a$Day) + labs(colour="Severity")
p <- p + theme(legend.position="bottom", legend.box = "horizontal") + guides(colour = guide_legend(nrow = 2))
p<-p + scale_fill_manual(values=c("Mild"="blue","Moderate"="orange","Severe"="#E34234","Critical"="#953553","COVID"="yellow")) #"Control"="#008080","Asymptomatic"="#87CEEB","Fatal"="black",
p<-p + scale_color_manual(values=c("Mild"="blue","Moderate"="orange","Severe"="#E34234","Critical"="#953553","COVID"="yellow")) #"Control"="#008080","Asymptomatic"="#87CEEB","Fatal"="black",

require(scales)
p<-p + scale_x_continuous(limits = c(-10, 280), expand = c(0.05,0.05), trans = log2_trans(), breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(2^.x))) #expand = c(0,0), 
p

library("ggExtra")
ggMarginal(p, type = "density",colour = "green", fill = "green",margins = "y",alpha=0.05)

require(cowplot)
ybox <- axis_canvas(p, axis = "y") + 
  geom_boxplot(data=a, aes(y=bray.rel, x=factor(Severity), fill=factor(Severity), color=factor(Severity)),alpha=0.1) +
  scale_color_manual(values=c("Control"="#008080","Asymptomatic"="#87CEEB","Mild"="blue","Moderate"="orange","Severe"="#E34234","Critical"="#953553","Fatal"="black","COVID"="yellow","non-ICU"="#8a75ba","ICU"="#e334b7")) +
  scale_fill_manual(values=c("Control"="#008080","Asymptomatic"="#87CEEB","Mild"="blue","Moderate"="orange","Severe"="#E34234","Critical"="#953553","Fatal"="black","COVID"="yellow","non-ICU"="#8a75ba","ICU"="#e334b7")) +
  scale_x_discrete()

p1 <- insert_yaxis_grob(p, ybox, grid::unit(0.7, "in"), position = "right")
ggdraw(p1)


#################################################################################################
#linear fit to each disease severity: lm(within-subject Bray-Curtis dissimilarities ~ log(Day, 2)
#################################################################################################
a$Severity <- factor(a$Severity, levels=c("Mild","Moderate","Severe","Critical","COVID"))

fit <- lm(bray.rel~log(Day,2), data=a)
summary(fit)

a1 <- a[which(a$Severity == 'COVID'),]
fit1 <- lm(bray.rel~log(Day,2), data=a1)
summary(fit1)

a2 <- a[which(a$Severity == 'Mild'),]
fit2 <- lm(bray.rel~log(Day,2), data=a2)
summary(fit2)

a3 <- a[which(a$Severity == 'Moderate'),]
fit3 <- lm(bray.rel~log(Day,2), data=a3)
summary(fit3)

a4 <- a[which(a$Severity == 'Severe'),]
fit4 <- lm(bray.rel~log(Day,2), data=a4)
summary(fit4)

a5 <- a[which(a$Severity == 'Critical'),]
fit5 <- lm(bray.rel~log(Day,2), data=a5)
summary(fit5)




###################################
#  longitudinal dysbiosis scores
###################################

# load dysbiosis score: median Bray-Curtis distance of a COVID microbiome to that of every healthy control
a <-read.table('bray.rel.median.txt', head=T, sep='\t', check.names = F) 

# subset individuals with multiple samples only
a <- a[grepl('Yes', a$'Longitudinal'), ]

a$Day<-as.numeric(as.character(a$Day))+1
a$bray.rel<-as.numeric(as.character(a$bray.rel))

p <- ggplot(a, aes(Day, bray.rel))
p <- p + geom_point(alpha=0.2,size=2)
p <- p + ggpubr::stat_cor(method = "spearman", label.x = 0.1, color='black', size=3) 
p <- p + geom_smooth(method="lm", se=F, colour = 'red', size=1, alpha=0.2) 
p <- p + geom_smooth(method="loess", se=F, colour = 'blue', size=1, alpha=0.2, span=1)
p <- p + labs(y = "Dysbiosis score") + labs(x = "Day") + theme_classic() + theme(legend.justification = "bottom")
p <- p + theme(text = element_text(size=13), axis.text.x = element_text(angle=0, hjust=0.5))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black"))
p<-p + facet_wrap(. ~ Severity, scales="fixed",ncol = 6)
p <- p + expand_limits(y=c(0.6, 1.05))
p <- p + scale_y_continuous(expand = c(0.01,0), breaks=c(0.6,0.7,0.8,0.9,1.0))
p <- p + scale_fill_discrete(name="Day", breaks=a$Day)
p <- p + theme(legend.position="bottom", legend.box = "horizontal") + guides(colour = guide_legend(nrow = 2))
require(scales)
p<-p + scale_x_continuous(expand = c(0.1,0.1), trans = log2_trans(), breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(2^.x)))
p
#########################################################################################################
#########################################################################################################


