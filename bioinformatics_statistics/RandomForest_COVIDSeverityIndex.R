# load necessary packages
require(pROC)
require(randomForest)
require(ggplot2)

df<-read.csv('../data/merge_species_rel_meta.txt', header=T, sep='	', stringsAsFactors = FALSE, check.names = F,row.names = 1)

# subset biomarkers
d<-subset(select=c('Condition','Cohort','Clostridium_innocuum','Enterococcus_faecalis','Enterococcus_faecium','Staphylococcus_epidermidis','Faecalibacterium_prausnitzii','Actinomyces_odontolyticus','Adlercreutzia_equolifaciens','Anaerostipes_hadrus','Bacteroides_eggerthii','Bacteroides_ovatus','Bacteroides_stercoris','Bacteroides_xylanisolvens','Barnesiella_intestinihominis','Bifidobacterium_adolescentis','Bifidobacterium_bifidum','Bifidobacterium_longum','Bifidobacterium_pseudocatenulatum','Blautia_wexlerae','Collinsella_aerofaciens','Coprococcus_catus','Coprococcus_comes','Coprococcus_eutactus','Dorea_formicigenerans','Dorea_longicatena','Eubacterium_eligens','Eubacterium_hallii','Eubacterium_ramulus','Eubacterium_rectale','Eubacterium_ventriosum','Fusicatenibacter_saccharivorans','Lachnospira_pectinoschiza','Paraprevotella_clara','Parasutterella_excrementihominis','Roseburia_faecis','Roseburia_hominis','Roseburia_intestinalis','Roseburia_inulinivorans','Ruminococcus_bromii','Ruminococcus_lactaris','Ruminococcus_torques'), df)
d[1:5,1:9]

# train data
##############
df_train <- d[which(d$Cohort == 'cohort1'),]
df_train$Condition
df_train[1:5,1:5]
df_train <- transform(df_train, Condition=as.factor(Condition))
str(df_train[,1:5])

set.seed(123) # for replicated results
model <- randomForest(formula = Condition ~ ., data = df_train[,-c(2:3)], importance=TRUE)

# test data
###################
df_test<-d[which(d$Cohort == 'cohort2'),]
df_test[1:5,1:5]
dim(df_test)
df_test <- transform(df_test, Condition=as.factor(Condition))
str(df_test[,1:5])

set.seed(123)
predictions <-predict(model, df_test[,-c(2:3)], type = "prob")[, 'COVID']

#model performace: AUC value
roc<-roc(df_test$Condition ~ predictions, plot = F, print.auc = TRUE)
sprintf("%0.4f", auc(roc)) 

##########
#ROC plot
##########
roc1 <- roc(df_test$Condition, predictions, plot=TRUE, auc.polygon=F, max.auc.polygon=F, grid=F, xlab="False Positive Rate", ylab="True Postive Rate", print.auc=TRUE, show.thres=T)
ci1 <- ci.se(roc1, specificities=seq(0, 1, len=20))
plot(ci1, type="shape", col="lightblue", xlim = c(0, 1))
ci(roc1) #95% CI:xx-xx
text(0.27, 0.4, '95% CI:xx-xx')


##################################################
#probability of being COVID: COVID Severity Index
##################################################
data <- data.frame(Index = predictions, Group = df_test$Condition)
data1<-merge(df_test[,2:3],data,by=c('row.names'),all=T)
data1$Severity <- factor(data1$Severity, levels=c("Control","Mild","Moderate","Severe","Critical","Fatal","COVID"))

#COVID Severity Index plot
ggplot(data1, aes(x = Severity, y = predictions, col = Severity)) +
  geom_boxplot() + theme_classic() + labs(x = "") + labs(y = "COVID Severity Index") +
  geom_jitter(width = 0.2, size = 1.5, alpha=0.5, fill='lightgrey', shape=21) +
  scale_fill_manual(values=c("Control"="#008080","Asymptomatic"="#87CEEB","Mild"="blue","Moderate"="orange","Severe"="#E34234","Critical"="#953553","Fatal"="black","COVID"="yellow")) +
  scale_color_manual(values=c("Control"="#008080","Asymptomatic"="#87CEEB","Mild"="blue","Moderate"="orange","Severe"="#E34234","Critical"="#953553","Fatal"="black","COVID"="yellow")) +
  scale_y_continuous(limits = c(0,1.02), expand = c(0,0), breaks=c(0,0.2,0.4,0.6,0.8,1)) +
  theme(text = element_text(size=14), axis.text.x = element_text(angle=90, hjust=0.5)) +
  theme(legend.position="NULL")

