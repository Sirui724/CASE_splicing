a <- read.table('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\some_file\\up_vs_down.txt',sep='\t',head=T)

a$cancer = factor(a$cancer, levels=c('STAD','LUAD*','LUAD','HNSC','ESCC*','BLCA','PRAD','BRCA','ESCA','GBM','UCEC','READ','KIRP','KICH','LUSC','KIRC','LIHC','COAD','PAAD','THCA'))
a$type = factor(a$type,levels=c('short','long'))
p <- ggplot(a,aes(x =cancer, y = up_vs_down))+geom_point(aes(colour = factor(type),shape=factor(type)),size=3)+ theme_bw()+scale_color_brewer(palette = "Accent")+theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
p+theme(panel.grid.major=element_line(colour=NA))



######number
a <- read.table('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\some_file\\number.txt',sep='\t',head=TRUE)
ggplot(data=a,mapping=aes(x=as.factor(a$number),y=a$short))+geom_bar(stat="identity",width=0.7)+theme_classic()+xlab("number of cancer type")+ylab("number of AS events")
ggplot(data=a,mapping=aes(x=as.factor(a$number),y=a$all))+geom_bar(stat="identity",width=0.7)+theme_classic()+xlab("number of cancer type")+ylab("number of AS events")

library(reshape2)
a$long=a$all-a$short
b <- as.data.frame(a[,-3])
b <- melt(b,id.vars = c("number"))
b$variable <- factor(b$variable,levels=c('long','short'))
ggplot(b, aes(x =as.factor(b$number), y = b$value, fill =as.factor(b$variable)) )+geom_bar(stat = "identity", position = "stack",width=0.6)+theme_bw()+theme(panel.grid=element_blank()) +scale_fill_manual(values=c('grey37','red3'))+ guides(fill=guide_legend(title=NULL))+xlab("Number of cancer type")+ylab("Number of ASE")+theme_classic()
b$variable <- factor(b$variable,levels=c('short','long'))
ggplot(b, aes(x =as.factor(b$number), y = b$value, fill =as.factor(b$variable)) )+geom_bar(stat = "identity", position = "stack",width=0.6)+theme_bw()+theme(panel.grid=element_blank()) +scale_fill_manual(values=c('red3','grey37'))+ guides(fill=guide_legend(title=NULL))+xlab("Number of cancer type")+ylab("Number of ASE")+theme_classic()


a <- read.table('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\some_file\\number.txt',sep='\t',head=TRUE)
aa <- a[,c(1,4,6)]
b <- melt(aa,id.vars = c("number"))
b$variable <- factor(b$variable,levels=c('short_','long_'))
ggplot(b, aes(x =as.factor(b$number), y = b$value, fill =as.factor(b$variable)) )+geom_bar(stat = "identity", position = "stack",width=0.6)+theme_bw()+theme(panel.grid=element_blank()) +scale_fill_manual(values=c('red3','grey37'))+ guides(fill=guide_legend(title=NULL))+xlab("Number of cancer type")+ylab("Number of ASE")+theme_classic()



###########merge all shortexon and all samples

COAD <- read.table('D:\\wanglab\\smallexon\\TCGA_smallexon\\AC_data\\COAD_AS.txt',head=TRUE)
COAD <- filter(COAD,COAD$exon_legth<=60)
COAD <- cbind(COAD[,2],COAD[,7:333])
rownames(COAD) <- COAD[,1]
COAD <- COAD[,2:326]
COAD$name <- rownames(COAD)
A <-COAD

a <- c('BLCA','BRCA','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LIHC','LUAD','LUSC')
a <- c('PAAD','PRAD','READ','STAD','THCA','UCEC')
for (i in a){
  z <- paste(i,sep="")
  z <-read.table(paste('D:\\wanglab\\smallexon\\TCGA_smallexon\\AC_data\\',i,'_AS.txt',sep=""),head=TRUE)
  z <- filter(z,z$exon_legth<=60)
  l <-length(z)-5
  z <- cbind(z[,2],z[,7:l])
  rownames(z) <- z[,1]
  ll=l-5
  z <- z[,2:ll]
  z$name<- rownames(z)
  A<- merge(A,z,by='name',all=TRUE)
}

write.table(t(A),file='D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\shortexon_merge_t.txt',sep='\t',row.names=TRUE)
write.table(A,file='D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\shortexon_merge.txt',sep='\t',row.names=TRUE)

##########

A <- read.table('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\shortexon_merge.txt',sep='\t',head=TRUE,row.names=1)
B<-t(A)
B <- as.data.frame(B)
data_AS <- read.table('D:\\wanglab\\smallexon\\TCGA_smallexon\\clinical_merge.txt',sep='\t',row.names=1,head=TRUE)
z <- row.names(B)
library(stringr)
zz <- paste(str_sub(z,start=1L,end=4L),str_sub(z,start=6L,end=7L),str_sub(z,start=9L,end=12L),str_sub(z,start=14L,end=15L),sep='-')
B$name=zz
data_AS$name <- rownames(data_AS)
A <- merge(data_AS,B,by='name')
write.table(A,file='D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\shortexon_merge_clinical.txt',sep='\t',row.names=TRUE)

#####################
D <- read.table('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\some_file\\shortexon_merge.txt',sep='\t',head=TRUE)
f_na<-function(x) sum(is.na(x))
na <- apply(D,1,f_na)
D$NA_number <- na
library(stringr)
tumor <- c()
normal <- c()

for (i in seq(from=1,to=7509)){
  if (str_sub(colnames(D)[i],start=14L,end=15L)=='01'){
    tumor<- c(tumor,i)
  }
  else if (str_sub(colnames(D)[i],start=14L,end=15L)=='11'){
    normal<-c(normal,i)
  }
}
length(tumor)
length(normal)

f_na<-function(x) sum(is.na(x))
na <- apply(D[,tumor],1,f_na)
D$NA_number_t <- na
na <- apply(D[,normal],1,f_na)
D$NA_number_n <- na

library(dplyr)
D_f <- filter(D,D$NA_number_t<=5500)
D_f <- filter(D_f,D_f$NA_number_n<=450)


delt_psi <- apply(D_f,1,function(x) mean(as.numeric(x[tumor]),na.rm = T)-mean(as.numeric(x[normal]),na.rm = T))
D_f$delt_psi <- delt_psi
D_f2 <- filter(D_f,abs(delt_psi)>0.1)
p_value <- apply(D_f2,1,function(x) t.test(as.numeric(x[tumor]),as.numeric(x[normal]))$p.value)
D_f2$p_value <- p_value
D_f2 <- filter(D_f2,p_value<0.001)

TEMP <-c()
a <- apply(D_f2,2,f_na)
for (i in seq(from=1,to=7509)){
  if (a[i]<=20){
    TEMP<- c(TEMP,i)
  }
}

D_f3 <- D_f2[,TEMP]
tumor <- c()
normal <- c()

for (i in seq(from=1,to=length(D_f3))){
  if (str_sub(colnames(D_f3)[i],start=14L,end=15L)=='01'){
    tumor<- c(tumor,i)
  }
  else if (str_sub(colnames(D_f3)[i],start=14L,end=15L)=='11'){
    normal<-c(normal,i)
  }
}
length(tumor)
length(normal)


#D_f2[is.na(D_f2)]=mean(as.matrix(D_f2[,2:7529]),na.rm=T)
D_f3[is.na(D_f3)]=mean(as.matrix(D_f3[,2:length(D_f3)]),na.rm=T)

library(mixOmics)
trans.pca <- pca(t(cbind(as.matrix(D_f3[,tumor]),as.matrix(D_f3[,normal]))),ncomp = 7, center = TRUE, scale = TRUE)
plot(trans.pca)
class <- c(rep('tumor',length(tumor)),rep('normal',length(normal)))
plotIndiv(trans.pca, comp = c(1, 2), ind.names = F, group = class,legend = F, title = 'pca',col=c('red','steelblue1'),ellipse = F)
plotIndiv(trans.pca, comp = c(1, 2), ind.names = F, group = class,legend = T, title = 'PCA',col=c('steelblue1','darkorange'),ellipse =F,pch=16,cex=0.5,style = 'lattice',ellipse.level = 0.9,legend.title = "")
plotIndiv(trans.pca, comp = c(1, 2), ind.names = F, group = class,legend = T, title = 'PCA',col=c('steelblue','darkorange2'),ellipse =T,pch=16,cex=0.4,ellipse.level = 0.8,legend.title = "")

trans.plsda <- plsda(t(cbind(as.matrix(D_f3[,tumor]),as.matrix(D_f3[,normal]))),class,ncomp=7)
plotIndiv(trans.plsda, comp = c(1, 2), ind.names = F, group = class,legend = T, title = 'PLS-DA',col=c('steelblue1','darkorange'),ellipse =F,pch=16,cex=0.5,style = 'lattice',ellipse.level = 0.9,legend.title = "")
plotIndiv(trans.plsda, comp = c(1, 2), ind.names = F, group = class,legend = T, title = 'PLS-DA',col=c('steelblue','darkorange2'),ellipse =T,pch=16,cex=0.4,ellipse.level = 0.8,legend.title = "")

#trans.plsda <- plsda(d_rf_2,as.factor(d_rf$type),ncomp=7)
#data_pre <- LUAD_t
#d_rf_2 <- d_rf_2[,order(names(d_rf_2))]
#data_pre <- as.data.frame(lapply(data_pre,as.numeric))
#data_pre <-data_pre[,order(names(data_pre))]
#pred <- predict(trans.plsda,data_pre)


library(randomForest)
library(pROC)
d_rf <- D_f3
d_rf <- as.data.frame(t(d_rf),stringsAsFactors=FALSE)
colnames(d_rf) <- d_rf[1,]
d_rf <- d_rf[-1,]
d_rf$type <- str_sub(rownames(d_rf),start=14L,end=15L)

d_rf[is.na(d_rf)]=mean(as.numeric(as.matrix(d_rf[,1:47])),na.rm=T)

n <-c()
for (i in seq(from=1,to=46)){
  n <- c(n,str_c(str_split_fixed(colnames(d_rf)[i],':',9)[1,5],str_split_fixed(colnames(d_rf)[i],':',9)[1,6],sep='.'))
}
n <- c(n,'type')
colnames(d_rf) <- n


d_rf <- as.data.frame(lapply(d_rf,as.numeric))
d_rf <- filter(d_rf,d_rf$type=='1' | d_rf$type=='11')
rf <- randomForest(as.factor(d_rf$type) ~., data=d_rf, importance=TRUE,proximity=TRUE,mtry=3,ntree=1300)
plot(rf)
ntree=1300

set.seed(99)
n <-length(names(d_rf))
for (i in 1:(n-1)){
  mtry_fit <- randomForest(as.factor(d_rf$type) ~., data=d_rf,mtry=i)
  err <- mean(mtry_fit$err.rate)
  print(i)
  print(err)
}

rf <- randomForest(as.factor(d_rf$type) ~., data=d_rf, importance=TRUE,proximity=TRUE,mtry=22,ntree=1300)
plot(rf)

importance(rf)
varImpPlot(rf,n.var=10)


train_sub = sample(nrow(d_rf),7/10*nrow(d_rf))
train_data = d_rf[train_sub,]
test_data = d_rf[-train_sub,]
rf <- randomForest(as.factor(train_data$type) ~., data=train_data, importance=TRUE,proximity=TRUE,mtry=22,ntree=1300)
plot(rf)


pre_ran <- predict(rf,newdata=test_data)
obs_p_ran = data.frame(prob=pre_ran,obs=test_data$type)
table(test_data$type,pre_ran,dnn=c("真实值","预测值"))
ran_roc <- roc(test_data$type,as.numeric(pre_ran))
plot(ran_roc, print.auc=T, auc.polygon=TRUE,grid=c(1,1),grid.col='grey', max.auc.polygon=T,auc.polygon.col="pink", print.thres=TRUE)

test_all <- c()
pre_all <-c()
auc_all <- c()
for (i in 1:100){
  set.seed(i)
  train_sub = sample(nrow(d_rf),7/10*nrow(d_rf))
  train_data = d_rf[train_sub,]
  test_data = d_rf[-train_sub,]
  rf <- randomForest(as.factor(train_data$type) ~., data=train_data, importance=TRUE,proximity=TRUE,mtry=22,ntree=1300)
  pre_ran <- predict(rf,newdata=test_data)
 # obs_p_ran = data.frame(prob=pre_ran,obs=test_data$type)
 # table(test_data$type,pre_ran,dnn=c("真实值","预测值"))
  ran_roc <- roc(test_data$type,as.numeric(pre_ran))
 # print(ran_roc$auc)
  test_all <- c(test_all,test_data$type)
  pre_all <- c(pre_all,pre_ran)
  auc_all <- c(auc_all,as.numeric(ran_roc$auc))
}
ran_roc <- roc(test_all,pre_all)
plot(ran_roc, print.auc=T, auc.polygon=TRUE,grid=c(1,1),grid.col='grey', max.auc.polygon=T,auc.polygon.col="pink", print.thres=TRUE)


#################



####################test by LUAD data from wyb
LUAD <- read.table("D:\\wanglab\\project_wyb\\clinical\\clinical_miso_analysis\\miso_out_NA_gene.txt",head=TRUE)
LUAD_ <- filter(LUAD,LUAD$event_name%in%D_f2$name)
data(ROCR.simple)
df <- data.frame(ROCR.simple)
pred <- prediction(df$predictions, df$labels)
perf <- performance(pred,"tpr","fpr")
plot(perf,colorize=TRUE)
LUAD_ <- LUAD_[,2:18]
n <-c()
for (i in seq(from=1,to=26)){
  n <- c(n,str_c('X',str_c(str_split_fixed(LUAD_$event_name[i],':',9)[1,5],str_split_fixed(LUAD_$event_name[i],':',9)[1,6],sep='.'),sep=''))
}
LUAD_$n <- n

LUAD_t <-as.data.frame(t(LUAD_),stringsAsFactors=FALSE)
colnames(LUAD_t) <- LUAD_t[18,]
LUAD_t <-LUAD_t[2:17,]


d_rf_2 <-d_rf[,which(names(d_rf)%in%colnames(LUAD_t))]
rf <- randomForest(as.factor(d_rf$type) ~., data=d_rf_2, importance=TRUE,proximity=TRUE)
plot(rf)
rf

set.seed(99)
n <-length(names(d_rf_2))
for (i in 1:(n-1)){
  mtry_fit <- randomForest(as.factor(d_rf$type) ~., data=d_rf,mtry=i)
  err <- mean(mtry_fit$err.rate)
  print(i)
  print(err)
}

LUAD_t[is.na(LUAD_t)]=mean(as.numeric(as.matrix(LUAD_t)),na.rm=T)

pre_ran <- predict(rf,newdata=LUAD_t)
obs_p_ran = data.frame(prob=pre_ran,obs=c(1,1,1,11,11,11,11,11,1,11,11,1,1,11,1,1))
table(c(1,1,1,11,11,11,11,11,1,11,11,1,1,11,1,1),pre_ran,dnn=c("真实值","预测值"))
ran_roc <- roc(c(1,1,1,11,11,11,11,11,1,11,11,1,1,11,1,1),as.numeric(pre_ran))
plot(ran_roc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),grid.col=c("green", "red"), max.auc.polygon=TRUE,auc.polygon.col="skyblue", print.thres=TRUE)


###############################################################
########test by ESCC data
ESCC <- read.table("D:\\wanglab\\食管癌组学数据分析\\analysis\\a.txt",head=TRUE)
ESCC_ <- filter(ESCC,ESCC$event_name%in%D_f2$name)
ESCC_ <- filter(ESCC_,ESCC_$NA_number<=20)
ESCC_ <- ESCC_[,2:180]
n <-c()
for (i in seq(from=1,to=32)){
  n <- c(n,str_c('X',str_c(str_split_fixed(ESCC_$event_name[i],':',9)[1,5],str_split_fixed(ESCC_$event_name[i],':',9)[1,6],sep='.'),sep=''))
}
ESCC_$n <- n


ESCC_t <-as.data.frame(t(ESCC_),stringsAsFactors=FALSE)
colnames(ESCC_t) <- ESCC_t[180,]
ESCC_t <-ESCC_t[2:179,]


d_rf_3 <-d_rf[,which(names(d_rf)%in%colnames(ESCC_t))]
rf <- randomForest(as.factor(d_rf$type) ~., data=d_rf_3, importance=TRUE,proximity=TRUE,mtry=9,ntree=1500)
plot(rf)
rf


ESCC_t[is.na(ESCC_t)]=mean(as.numeric(as.matrix(ESCC_t)),na.rm=T)

pre_ran <- predict(rf,newdata=ESCC_t)
t <- rep(c(11,1),89)
obs_p_ran = data.frame(prob=pre_ran,obs=t)
table(t,pre_ran,dnn=c("真实值","预测值"))
ran_roc <- roc(t,as.numeric(pre_ran))
plot(ran_roc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),grid.col=c("green", "red"), max.auc.polygon=TRUE,auc.polygon.col="skyblue", print.thres=TRUE)



############################################################
##########ESCC only
d_ESCC <- read.table('D:\\wanglab\\smallexon\\TCGA_smallexon\\AC_data\\ESCA_AS.txt',head=TRUE)
d_ESCC <-  filter(d_ESCC,d_ESCC$NA_number_t<=30)
d_ESCC <-  filter(d_ESCC,d_ESCC$NA_number_n<=3)
d_ESCC <- filter(d_ESCC,abs(d_ESCC$delt_psi)>=0.2)
d_E <- cbind(d_ESCC$name,d_ESCC[,7:202])
d_E2 <- as.data.frame(t(d_E),stringsAsFactors=FALSE)
colnames(d_E2) <- d_E2[1,]
d_E2 <- d_E2[-1,]
d_E2$type <- str_sub(rownames(d_E2),start=14L,end=15L)
d_E2[is.na(d_E2)]=mean(as.numeric(as.matrix(d_E2[,1:36])),na.rm=T)


n <-c()
for (i in seq(from=1,to=36)){
  n <- c(n,str_c('X',str_c(str_split_fixed(colnames(d_E2)[i],':',9)[1,5],str_split_fixed(colnames(d_E2)[i],':',9)[1,6],sep='.'),sep=''))
}
n <- c(n,'type')
colnames(d_E2) <- n

d_E2 <- as.data.frame(lapply(d_E2,as.numeric)) 
d_E2 <- filter(d_E2,d_E2$type=='1' | d_E2$type=='11')
rf <- randomForest(as.factor(d_E2$type) ~., data=d_E2, importance=TRUE,proximity=TRUE,mtry=18,ntree=500)
plot(rf)
rf
importance(rf)
varImpPlot(rf,n.var=10)

t <- rep(c(11,1),89)
ESCC_t2 <-ESCC_t[,which(names(ESCC_t)%in%names(d_E2))]
rf <- randomForest(as.factor(t) ~., data=ESCC_t2, importance=TRUE,proximity=TRUE,mtry=3,ntree=400)



pre_ran <- predict(rf,newdata=d_E2)
#t <- rep(c(11,1),89)
obs_p_ran = data.frame(prob=pre_ran,obs=d_E2$type)
table(d_E2$type,pre_ran,dnn=c("真实值","预测值"))
ran_roc <- roc(d_E2$type,as.numeric(pre_ran))
plot(ran_roc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),grid.col=c("green", "red"), max.auc.polygon=TRUE,auc.polygon.col="skyblue", print.thres=TRUE)

#############################


data_clinical <- read.table('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\some_file\\shortexon_merge_clinical.txt',sep='\t',row.names=1,head=TRUE)

library(survival)
library(survminer)

A=c()
C=c()
for (i in seq(from=10,to=526)){
  d <- data_clinical[order(data_clinical[,i]),]
  #print(colnames(d)[i])
  d <- filter(d,as.numeric(d[,i])>=0)
  aa <- round(length(d[,1])/4)
  class <- c(rep('low',aa),rep(NA,length(d[,1])-2*aa),rep('high',aa))
  d$psi_ <- class
  fit <- coxph(Surv(OS_month, OSS_) ~ class, data = d)
  A <- cbind(i,colnames(d)[i],as.numeric(coef(summary(fit))[,5]),as.numeric(round(exp(confint(fit)), 2)[1]),as.numeric(round(exp(confint(fit)), 2)[2]),as.numeric(round(exp(coef(fit)), 2)))
  C<-rbind(C,A)
}

C <- as.data.frame(C,stringsAsFactors = F)
colnames(C) <- c('i','name','p','low','high','HR')
C$p <- as.numeric(C$p)
C$HR<- as.numeric(C$HR)
C$low<- as.numeric(C$low)
C$high<- as.numeric(C$high)
library(forestplot)
C <- C[order(C[,6]),]
#C <- filter(C, C$p<0.05)
C$logp <- -log10(C$p+0.000001)

write.table(C,file='D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\hr_clinical_2.txt',sep='\t',row.names=FALSE)
CC <- read.table('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\hr_clinical_2.txt',sep='\t',head=TRUE,stringsAsFactors = F)

name <- read.table('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\some_file\\517_shortexon.txt',head=TRUE)

n <-c()
for (i in seq(from=1,to=517)){
  n <- c(n,str_c(str_split_fixed(name$name[i],':',9)[1,1],str_c(str_split_fixed(name$name[i],':',9)[1,5], str_split_fixed(name$name[i],':',9)[1,6],sep='-'),sep=':'))
}
name <- as.data.frame(name[,3])
name$name_2 <- n

name <- unique(name)
CC_ <-merge(CC,name,by="name_2")
colnames(CC_)[12] <- 'gene_name'
CC_$gene_name <- as.character(CC_$gene_name)

CC_ <- CC_[order(CC_$HR),]
CC_ <- filter(CC_,CC_$p<=0.05)

CC_ <- filter(CC_,CC_$p<=0.01)
CC_ <- filter(CC_,CC_$p<=0.001)
forestplot(cbind(CC_$gene_name,CC_$name_2,CC_$HR_,CC_$p),mean=log(CC_$HR),upper=log(CC_$high),lower=log(CC_$low),boxsize=0.3,zero=0,graph.pos=3,graphwidth = unit(70,"mm"),lineheight="auto",col = fpColors(lines="gray45", box="gray7"),txt_gp=fpTxtGp(label=gpar(cex=0.6), ticks=gpar(cex=0.7)),ci.vertices=TRUE, ci.vertices.height = 0.2)
CC_2 <- filter(CC_,CC_$p<=0.0000000001)
forestplot(cbind(CC_2$gene_name,CC_2$name_2,CC_2$HR_,CC_2$p),mean=log(CC_2$HR),upper=log(CC_2$high),lower=log(CC_2$low),boxsize=0.3,zero=0,graph.pos=3,graphwidth = unit(70,"mm"),lineheight="auto",col = fpColors(lines="gray45", box="gray7"),txt_gp=fpTxtGp(label=gpar(cex=0.6), ticks=gpar(cex=0.7)),ci.vertices=TRUE, ci.vertices.height = 0.2)
CC_2 <- filter(CC_,CC_$p<=0.00001)
forestplot(cbind(CC_2$gene_name,CC_2$HR_,CC_2$p),mean=log(CC_2$HR),upper=log(CC_2$high),lower=log(CC_2$low),boxsize=0.3,zero=0,graph.pos=2,graphwidth = unit(40,"mm"),lineheight="auto",col = fpColors(lines="gray45", box="red3"),txt_gp=fpTxtGp(label=gpar(cex=0.7), ticks=gpar(cex=0.65)),ci.vertices=TRUE, ci.vertices.height = 0.2)

CC_$log_HR <- log(CC_$HR)
CC_$log_low <- log(CC_$low)
CC_$log_high <- log(CC_$high)
CC_2 <- filter(CC_,CC_$p<=0.0000000001)



plot(sort(CC$logP),pch=16,cex=0.4,xlab="",ylab='log(P)')
abline(h=0)
abline(h=1.3)
abline(h=-1.3)


for (i in seq(from=10,to=526)){
  d <- data_clinical[order(data_clinical[,i]),]
  print(colnames(d)[i])
  d <- filter(d,as.numeric(d[,i])>=0)
  aa <- round(length(d[,1])/4)
  class <- c(rep('low',aa),rep(NA,length(d[,1])-2*aa),rep('high',aa))
  d$psi_ <- class
  fit <-survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
  ggsurvplot(fit,
             pval = TRUE, conf.int = FALSE,
             risk.table = FALSE, # Add risk table
             risk.table.col = "strata", # Change risk table color by groups
             # linetype = "strata", # Change line type by groups
             #  surv.median.line = "hv", # Specify median survival
             # ggtheme = theme_bw(legend.position = "right"), # Change ggplot2 theme
             palette = c("#E7B800", "#2E9FDF"))
  ggsave(paste('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\km_plot\\',colnames(d)[i],'.pdf',sep=""))
}



i=191
d <- data_clinical[order(data_clinical[,i]),]
print(colnames(d)[i])
d <- filter(d,as.numeric(d[,i])>=0)
aa <- round(length(d[,1])/4)
class <- c(rep('low',aa),rep(NA,length(d[,1])-2*aa),rep('high',aa))
d$psi_ <- class
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,
           pval = TRUE, conf.int = FALSE,
           risk.table = FALSE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           # linetype = "strata", # Change line type by groups
           #  surv.median.line = "hv", # Specify median survival
           # ggtheme = theme_bw(legend.position = "right"), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))
an <- ggsurvplot(fit,palette = c("#E7B800", "#2E9FDF"))
an$plot + theme (legend.position = "right")

pdf("D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\km_plot\\myplot.pdf")
myplot <- ggsurvplot(fit,
                     pval = TRUE, conf.int = FALSE,
                     risk.table = FALSE, # Add risk table
                     risk.table.col = "strata", # Change risk table color by groups
                     # linetype = "strata", # Change line type by groups
                     #  surv.median.line = "hv", # Specify median survival
                     # ggtheme = theme_bw(legend.position = "right"), # Change ggplot2 theme
                     palette = c("#E7B800", "#2E9FDF"))
print(myplot)
dev.off()

#####################long exon 

COAD <- read.table('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\file\\COAD_AS.txt',head=TRUE)
COAD <- filter(COAD,COAD$exon_legth>60)
COAD <- cbind(COAD[,2],COAD[,7:333])
rownames(COAD) <- COAD[,1]
COAD <- COAD[,2:326]
COAD$name <- rownames(COAD)
A <-COAD

a <- c('BLCA','BRCA','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LIHC','LUAD','LUSC')
a <- c('PAAD','PRAD','READ','STAD','THCA','UCEC')
for (i in a){
  z <- paste(i,sep="")
  z <-read.table(paste('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\file\\',i,'_AS.txt',sep=""),head=TRUE)
  z <- filter(z,z$exon_legth>60)
  l <-length(z)-5
  z <- cbind(z[,2],z[,7:l])
  rownames(z) <- z[,1]
  ll=l-5
  z <- z[,2:ll]
  z$name<- rownames(z)
  A<- merge(A,z,by='name',all=TRUE)
}

write.table(t(A),file='D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\file\\longexon_merge_t.txt',sep='\t',row.names=TRUE)
write.table(A,file='D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\file\\longexon_merge.txt',sep='\t',row.names=TRUE)


A <- read.table('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\file\\longexon_merge.txt',sep='\t',head=TRUE,row.names=1)
B<-t(A)
B <- as.data.frame(B)
data_AS <- read.table('D:\\wanglab\\smallexon\\TCGA_smallexon\\clinical_merge.txt',sep='\t',row.names=1,head=TRUE)
z <- row.names(B)
library(stringr)
zz <- paste(str_sub(z,start=1L,end=4L),str_sub(z,start=6L,end=7L),str_sub(z,start=9L,end=12L),str_sub(z,start=14L,end=15L),sep='-')
B$name=zz
data_AS$name <- rownames(data_AS)
A <- merge(data_AS,B,by='name')
write.table(A,file='D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\file\\longexon_merge_clinical.txt',sep='\t',row.names=TRUE)


data_clinical <- read.table('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\file\\longexon_merge_clinical.txt',sep='\t',row.names=1,head=TRUE)

library(survival)
library(survminer)

A=c()
C=c()
for (i in seq(from=10,to=1623)){
  d <- data_clinical[order(data_clinical[,i]),]
  #print(colnames(d)[i])
  d <- filter(d,as.numeric(d[,i])>=0)
  aa <- round(length(d[,1])/4)
  class <- c(rep('low',aa),rep(NA,length(d[,1])-2*aa),rep('high',aa))
  d$psi_ <- class
  fit <- coxph(Surv(OS_month, OSS_) ~ class, data = d)
  A <- cbind(i,colnames(d)[i],as.numeric(round(coef(summary(fit))[,5],10)),as.numeric(round(exp(confint(fit)), 2)[1]),as.numeric(round(exp(confint(fit)), 2)[2]),as.numeric(round(exp(coef(fit)), 2)))
  C<-rbind(C,A)
}

C <- as.data.frame(C,stringsAsFactors = F)
colnames(C) <- c('i','name','p','low','high','HR')
C$p <- as.numeric(C$p)
C$HR<- as.numeric(C$HR)
C$low<- as.numeric(C$low)
C$high<- as.numeric(C$high)
C$logp <- -log(C$p+0.000001)

write.table(C,file='D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\file\\hr_clinical_long.txt',sep='\t',row.names=FALSE)
CC <- read.table('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\some_file\\hr_clinical_long.txt',sep='\t',head=TRUE,stringsAsFactors = F)

CC_ <- filter(CC,CC$p<=0.05)

CC_ <- filter(CC_,CC_$p<=0.01)
CC_ <- filter(CC_,CC_$p<=0.001)

plot(sort(CC$logP),pch=16,cex=0.4,xlab="",ylab='log(P)')
abline(h=-1.30103)
abline(h=1.30103)
###############


data_cli <- read.table('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\some_file\\clinical_8.txt',sep='\t',row.names=1,head=TRUE)
m <- apply(data_cli,1,function(x) mean(as.numeric(x[c(17,18,19,20,21,22,23,24)]),na.rm=T))
data_cli$m8 <-m


i=26
d <- data_cli[order(data_cli[,i]),]
print(colnames(d)[i])
d <- filter(d,as.numeric(d[,i])>=0)
d <- filter(d,as.numeric(d$OS_month)>=0)
aa <- round(length(d[,1])/4)
bb <- round(length(d[,1])/2)
class <- c(rep('low',aa),rep(NA,length(d[,1])-2*aa),rep('high',aa)) 
class2<-c(rep('group1',aa),rep('group2',aa),rep('group3',aa),rep('group4',length(d[,1])-3*aa)) 
class3<-c(rep('low',bb),rep('high',length(d[,1])-bb)) 

d$psi_ <- class
d$psi_2 <- class2
d$psi_3<- class3
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,
           pval = TRUE, conf.int = FALSE,
           risk.table = FALSE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           # linetype = "strata", # Change line type by groups
           #  surv.median.line = "hv", # Specify median survival
           # ggtheme = theme_bw(legend.position = "right"), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))
an <- ggsurvplot(fit,palette = c("#E7B800", "#2E9FDF"))
an$plot + theme (legend.position = "right")

fit2 <- coxph(Surv(OS_month, OSS_) ~ class3, data = d)
cbind(as.numeric(round(coef(summary(fit2))[,5],10)),as.numeric(round(exp(confint(fit2)), 2)[1]),as.numeric(round(exp(confint(fit2)), 2)[2]),as.numeric(round(exp(coef(fit2)), 2)))



i=26
d<-filter(data_cli,data_cli$cancer_type=='blca_tcga_pan_can_atlas_2018')
d <- d[order(d[,i]),]
print(colnames(d)[i])
d <- filter(d,as.numeric(d[,i])>=0)
d <- filter(d,as.numeric(d$OS_month)>=0)
aa <- round(length(d[,1])/4)
class <- c(rep('low',aa),rep(NA,length(d[,1])-2*aa),rep('high',aa))
d$psi_ <- class
bb <- round(length(d[,1])/2)
class3<-c(rep('low',bb),rep('high',length(d[,1])-bb)) 
d$psi_3<- class3
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_3, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))


d<-filter(data_cli,data_cli$cancer_type=='brca_tcga')
d <- d[order(d[,i]),]
print(colnames(d)[i])
d <- filter(d,as.numeric(d[,i])>=0)
d <- filter(d,as.numeric(d$OS_month)>=0)
aa <- round(length(d[,1])/4)
class <- c(rep('low',aa),rep(NA,length(d[,1])-2*aa),rep('high',aa))
d$psi_ <- class
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))


d<-filter(data_cli,data_cli$cancer_type=='coadread_tcga')
d <- d[order(d[,i]),]
print(colnames(d)[i])
d <- filter(d,as.numeric(d[,i])>=0)
d <- filter(d,as.numeric(d$OS_month)>=0)
aa <- round(length(d[,1])/4)
class <- c(rep('low',aa),rep(NA,length(d[,1])-2*aa),rep('high',aa))
d$psi_ <- class
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))


d<-filter(data_cli,data_cli$cancer_type=='esca_tcga')
d <- d[order(d[,i]),]
print(colnames(d)[i])
d <- filter(d,as.numeric(d[,i])>=0)
d <- filter(d,as.numeric(d$OS_month)>=0)
aa <- round(length(d[,1])/4)
class <- c(rep('low',aa),rep(NA,length(d[,1])-2*aa),rep('high',aa))
d$psi_ <- class
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))


d<-filter(data_cli,data_cli$cancer_type=='gbm_tcga')
d <- d[order(d[,i]),]
print(colnames(d)[i])
d <- filter(d,as.numeric(d[,i])>=0)
d <- filter(d,as.numeric(d$OS_month)>=0)
aa <- round(length(d[,1])/4)
class <- c(rep('low',aa),rep(NA,length(d[,1])-2*aa),rep('high',aa))
d$psi_ <- class
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))


d<-filter(data_cli,data_cli$cancer_type=='hnsc_tcga')
d <- d[order(d[,i]),]
print(colnames(d)[i])
d <- filter(d,as.numeric(d[,i])>=0)
d <- filter(d,as.numeric(d$OS_month)>=0)
aa <- round(length(d[,1])/4)
class <- c(rep('low',aa),rep(NA,length(d[,1])-2*aa),rep('high',aa))
d$psi_ <- class
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))



d<-filter(data_cli,data_cli$cancer_type=='kich_tcga')
d <- d[order(d[,i]),]
print(colnames(d)[i])
d <- filter(d,as.numeric(d[,i])>=0)
d <- filter(d,as.numeric(d$OS_month)>=0)
aa <- round(length(d[,1])/4)
class <- c(rep('low',aa),rep(NA,length(d[,1])-2*aa),rep('high',aa))
aa <- round(length(d[,1])/2)
class <- c(rep('low',aa),rep('high',length(d[,1])-aa))
d$psi_ <- class
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))



d<-filter(data_cli,data_cli$cancer_type=='kirc_tcga')
d <- d[order(d[,i]),]
print(colnames(d)[i])
d <- filter(d,as.numeric(d[,i])>=0)
d <- filter(d,as.numeric(d$OS_month)>=0)
aa <- round(length(d[,1])/4)
class <- c(rep('low',aa),rep(NA,length(d[,1])-2*aa),rep('high',aa))
d$psi_ <- class
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))


d<-filter(data_cli,data_cli$cancer_type=='kirp_tcga')
d <- d[order(d[,i]),]
print(colnames(d)[i])
d <- filter(d,as.numeric(d[,i])>=0)
d <- filter(d,as.numeric(d$OS_month)>=0)
aa <- round(length(d[,1])/4)
class <- c(rep('low',aa),rep(NA,length(d[,1])-2*aa),rep('high',aa))
d$psi_ <- class
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))


d<-filter(data_cli,data_cli$cancer_type=='lihc_tcga')
d <- d[order(d[,i]),]
print(colnames(d)[i])
d <- filter(d,as.numeric(d[,i])>=0)
d <- filter(d,as.numeric(d$OS_month)>=0)
aa <- round(length(d[,1])/4)
class <- c(rep('low',aa),rep(NA,length(d[,1])-2*aa),rep('high',aa))
d$psi_ <- class
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))



d<-filter(data_cli,data_cli$cancer_type=='luad_tcga')
d <- d[order(d[,i]),]
print(colnames(d)[i])
d <- filter(d,as.numeric(d[,i])>=0)
d <- filter(d,as.numeric(d$OS_month)>=0)
aa <- round(length(d[,1])/4)
class <- c(rep('low',aa),rep(NA,length(d[,1])-2*aa),rep('high',aa))
d$psi_ <- class
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))


d<-filter(data_cli,data_cli$cancer_type=='lusc_tcga')
d <- d[order(d[,i]),]
print(colnames(d)[i])
d <- filter(d,as.numeric(d[,i])>=0)
d <- filter(d,as.numeric(d$OS_month)>=0)
aa <- round(length(d[,1])/4)
class <- c(rep('low',aa),rep(NA,length(d[,1])-2*aa),rep('high',aa))
d$psi_ <- class
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))


d<-filter(data_cli,data_cli$cancer_type=='paad_tcga')
d <- d[order(d[,i]),]
print(colnames(d)[i])
d <- filter(d,as.numeric(d[,i])>=0)
d <- filter(d,as.numeric(d$OS_month)>=0)
aa <- round(length(d[,1])/4)
class <- c(rep('low',aa),rep(NA,length(d[,1])-2*aa),rep('high',aa))
d$psi_ <- class
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))


d<-filter(data_cli,data_cli$cancer_type=='prad_tcga')
d <- d[order(d[,i]),]
print(colnames(d)[i])
d <- filter(d,as.numeric(d[,i])>=0)
d <- filter(d,as.numeric(d$OS_month)>=0)
aa <- round(length(d[,1])/4)
class <- c(rep('low',aa),rep(NA,length(d[,1])-2*aa),rep('high',aa))
d$psi_ <- class
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))


d<-filter(data_cli,data_cli$cancer_type=='stad_tcga')
d <- d[order(d[,i]),]
print(colnames(d)[i])
d <- filter(d,as.numeric(d[,i])>=0)
d <- filter(d,as.numeric(d$OS_month)>=0)
aa <- round(length(d[,1])/4)
class <- c(rep('low',aa),rep(NA,length(d[,1])-2*aa),rep('high',aa))
d$psi_ <- class
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))



d<-filter(data_cli,data_cli$cancer_type=='ucec_tcga')
d <- d[order(d[,i]),]
print(colnames(d)[i])
d <- filter(d,as.numeric(d[,i])>=0)
d <- filter(d,as.numeric(d$OS_month)>=0)
aa <- round(length(d[,1])/4)
class <- c(rep('low',aa),rep(NA,length(d[,1])-2*aa),rep('high',aa))
d$psi_ <- class
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))
######################

i=26
data_cli <- data_cli[order(data_cli[,i]),]
print(colnames(data_cli[i]))
data_cli <- filter(data_cli,as.numeric(data_cli[,i])>=0)
data_cli <- filter(data_cli,as.numeric(data_cli$OS_month)>=0)
aa <- round(length(data_cli[,1])/4)
class <- c(rep('low',aa),rep(NA,length(data_cli[,1])-2*aa),rep('high',aa))
data_cli$psi_ <- class
bb <- round(length(data_cli[,1])/2)
class3<-c(rep('low',bb),rep('high',length(data_cli[,1])-bb)) 
data_cli$psi_3<- class3

d<-filter(data_cli,data_cli$cancer_type=='blca_tcga_pan_can_atlas_2018')
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))

d<-filter(data_cli,data_cli$cancer_type=='brca_tcga')
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))

d<-filter(data_cli,data_cli$cancer_type=='coadread_tcga')
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))

d<-filter(data_cli,data_cli$cancer_type=='esca_tcga')
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))

d<-filter(data_cli,data_cli$cancer_type=='gbm_tcga')
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))

d<-filter(data_cli,data_cli$cancer_type=='hnsc_tcga')
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))

d<-filter(data_cli,data_cli$cancer_type=='kich_tcga')
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))


d<-filter(data_cli,data_cli$cancer_type=='kirc_tcga')
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))

d<-filter(data_cli,data_cli$cancer_type=='kirp_tcga')
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))


d<-filter(data_cli,data_cli$cancer_type=='lihc_tcga')
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))

d<-filter(data_cli,data_cli$cancer_type=='luad_tcga')
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))



d<-filter(data_cli,data_cli$cancer_type=='lusc_tcga')
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))



d<-filter(data_cli,data_cli$cancer_type=='paad_tcga')
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))



d<-filter(data_cli,data_cli$cancer_type=='prad_tcga')
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))



d<-filter(data_cli,data_cli$cancer_type=='stad_tcga')
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))



d<-filter(data_cli,data_cli$cancer_type=='ucec_tcga')
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,pval = TRUE, conf.int = FALSE,risk.table = FALSE,risk.table.col = "strata", palette = c("#E7B800", "#2E9FDF"))


##############
BLCA <- read.table('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\GBM_AS.txt',head=TRUE)
BLCA <- filter(BLCA,BLCA$exon_legth<=60)

library(randomForest)
library(pROC)

d_rf <- as.data.frame(t(BLCA[,7:(length(BLCA)-4)]),stringsAsFactors=FALSE)
colnames(d_rf) <- BLCA$name
d_rf[is.na(d_rf)]=mean(as.numeric(as.matrix(d_rf)),na.rm=T)
d_rf$type <- str_sub(rownames(d_rf),start=14L,end=15L)

n <-c()
for (i in seq(from=1,to=(length(d_rf)-1))){
  n <- c(n,str_c(str_split_fixed(colnames(d_rf)[i],':',9)[1,5],str_split_fixed(colnames(d_rf)[i],':',9)[1,6],sep='.'))
}
n <- c(n,'type')
colnames(d_rf) <- n


d_rf <- as.data.frame(lapply(d_rf,as.numeric))
d_rf <- filter(d_rf,d_rf$type=='1' | d_rf$type=='11')
rf <- randomForest(as.factor(d_rf$type) ~., data=d_rf, importance=TRUE,proximity=TRUE,mtry=3,ntree=1300)
plot(rf)
ntree=1300

set.seed(6)

train_sub = sample(nrow(d_rf),7/10*nrow(d_rf))
train_data = d_rf[train_sub,]
test_data = d_rf[-train_sub,]
rf <- randomForest(as.factor(train_data$type) ~., data=train_data, importance=TRUE,proximity=TRUE,mtry=22,ntree=1300)
plot(rf)


pre_ran <- predict(rf,newdata=test_data)
obs_p_ran = data.frame(prob=pre_ran,obs=test_data$type)
table(test_data$type,pre_ran,dnn=c("真实值","预测值"))
ran_roc <- roc(test_data$type,as.numeric(pre_ran))
plot(ran_roc, print.auc=T, auc.polygon=TRUE,grid=c(1,1),grid.col='grey', max.auc.polygon=T,auc.polygon.col="pink", print.thres=TRUE)

#每种cancer
#library(randomForest)
#library(pROC)
BLCA <- read.table('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\PRAD_AS.txt',head=TRUE)
BLCA <- filter(BLCA,BLCA$exon_legth<=60)

d_rf <- as.data.frame(t(BLCA[,7:(length(BLCA)-4)]),stringsAsFactors=FALSE)
colnames(d_rf) <- BLCA$name
d_rf[is.na(d_rf)]=mean(as.numeric(as.matrix(d_rf)),na.rm=T)
d_rf$type <- str_sub(rownames(d_rf),start=14L,end=15L)

n <-c()
for (i in seq(from=1,to=(length(d_rf)-1))){
  n <- c(n,str_c(str_split_fixed(colnames(d_rf)[i],':',9)[1,5],str_split_fixed(colnames(d_rf)[i],':',9)[1,6],sep='.'))
}
n <- c(n,'type')
colnames(d_rf) <- n

d_rf <- as.data.frame(lapply(d_rf,as.numeric))
d_rf <- filter(d_rf,d_rf$type=='1' | d_rf$type=='11')
rf <- randomForest(as.factor(d_rf$type) ~., data=d_rf, importance=TRUE,proximity=TRUE,mtry=3,ntree=1300)

test_all_ <- c()
pre_all_ <-c()
auc_all_ <- c()
for (i in 1:100){
  set.seed(i)
  train_sub = sample(nrow(d_rf),7/10*nrow(d_rf))
  train_data = d_rf[train_sub,]
  test_data = d_rf[-train_sub,]
  rf <- randomForest(as.factor(train_data$type) ~., data=train_data, importance=TRUE,proximity=TRUE,mtry=22,ntree=1300)
  pre_ran <- predict(rf,newdata=test_data)
  # obs_p_ran = data.frame(prob=pre_ran,obs=test_data$type)
  # table(test_data$type,pre_ran,dnn=c("真实值","预测值"))
  ran_roc <- roc(test_data$type,as.numeric(pre_ran))
  # print(ran_roc$auc)
  test_all_ <- c(test_all_,test_data$type)
  pre_all_ <- c(pre_all_,pre_ran)
  auc_all_ <- c(auc_all_,as.numeric(ran_roc$auc))
}
ran_roc_ <- roc(test_all_,pre_all_)
plot(ran_roc_, print.auc=T, auc.polygon=TRUE,grid=c(1,1),grid.col='grey', max.auc.polygon=T,auc.polygon.col="pink", print.thres=TRUE)
print(mean(auc_all_))


#######################
dd <- read.table('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\some_file\\AUC.txt',head=T,row.names = 1)
ggradar(dd)
dd <- dd[-1,]


#######
data_immue <- read.csv('D:\\wanglab\\smallexon\\new_60\\04immue\\infiltration_estimation_for_tcga.csv')
data_m8 <- as.data.frame(cbind(rownames(d),d$m8,d$psi_))
colnames(data_m8) <- c('sample_name','m8','group')
colnames(data_immue)[1] <- 'sample_name'
dd <- merge(data_immue,data_m8,by='sample_name')
library(dplyr)
dd$group <- as.character(dd$group)
dd_ <- filter(dd,dd$group=='high' | dd$group=='low')
row.names(dd_)=dd_$sample_name
dd_ <- dd_[,-1]
dd_timer <- cbind(dd_[,1:6],dd_$group)
library(reshape2)
dd_timer$`dd_$group`<-as.character(dd_timer$`dd_$group`)
dd_timer <- melt(dd_timer)
ggplot(data=dd_timer,aes(y=value,x=variable,fill=`dd_$group`)) + geom_boxplot()+stat_compare_means(aes(group=`dd_$group`), label = "p.signif",method='t.test')+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggplot(data=dd_timer,aes(y=value,x=variable,fill=`dd_$group`)) + geom_boxplot()+stat_compare_means(aes(group=`dd_$group`), label = "p.signif",method='t.test')+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+xlab("")+ylab("")+scale_fill_manual(values=c("#E7B800", "#2E9FDF"))+theme(panel.grid=element_blank())

dd_CIBERSORT <- cbind(dd_[,7:28],dd_$group)
dd_CIBERSORT$`dd_$group`<-as.character(dd_CIBERSORT$`dd_$group`)
dd_CIBERSORT <- melt(dd_CIBERSORT)

dd_XCELL <- cbind(dd_[,73:111],dd_$group)
dd_XCELL$`dd_$group`<-as.character(dd_XCELL$`dd_$group`)
dd_XCELL <- melt(dd_XCELL)
ggplot(data=dd_XCELL,aes(y=value,x=variable,fill=`dd_$group`)) + geom_boxplot()+stat_compare_means(aes(group=`dd_$group`), label = "p.signif",method='t.test')+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+xlab("")+ylab("")+scale_fill_manual(values=c("#E7B800", "#2E9FDF"))+theme(panel.grid=element_blank())

library(stringr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
dd_CIBERSORT<-separate(dd_CIBERSORT,variable, sep = "_", remove = TRUE,into = c("variable", "n"))
dd_CIBERSORT$`dd_$group`<- factor(dd_CIBERSORT$`dd_$group`,levels=c('low','high'))
dd_CIBERSORT$variable <- factor(dd_CIBERSORT$variable,levels=c('B.cell.naive','B.cell.memory','B.cell.plasma','T.cell.CD8.','T.cell.CD4..naive','T.cell.CD4..memory.resting','T.cell.CD4..memory.activated','T.cell.follicular.helper','T.cell.regulatory..Tregs.','T.cell.gamma.delta','NK.cell.resting','NK.cell.activated','Monocyte','Macrophage.M0','Macrophage.M1','Macrophage.M2','Myeloid.dendritic.cell.resting','Myeloid.dendritic.cell.activated','Mast.cell.resting','Mast.cell.activated','Eosinophil','Neutrophil'))
a1 <- filter(dd_CIBERSORT,dd_CIBERSORT$variable=='B.cell.naive' | dd_CIBERSORT$variable=='B.cell.memory' |dd_CIBERSORT$variable=='B.cell.plasma')
a2 <- filter(dd_CIBERSORT,dd_CIBERSORT$variable=='T.cell.CD4..memory.resting' | dd_CIBERSORT$variable=='T.cell.CD4..memory.activated' )
a3 <- filter(dd_CIBERSORT,dd_CIBERSORT$variable=='Mast.cell.resting' | dd_CIBERSORT$variable=='Mast.cell.activated')
a4 <- filter(dd_CIBERSORT,dd_CIBERSORT$variable=='NK.cell.resting' | dd_CIBERSORT$variable=='NK.cell.activated')
a5 <- filter(dd_CIBERSORT,dd_CIBERSORT$variable=='T.cell.follicular.helper' | dd_CIBERSORT$variable=='T.cell.regulatory..Tregs.')
a6 <-  filter(dd_CIBERSORT,dd_CIBERSORT$variable=='B.cell.naive' | dd_CIBERSORT$variable=='B.cell.memory' |dd_CIBERSORT$variable=='B.cell.plasma' |dd_CIBERSORT$variable=='T.cell.follicular.helper' | dd_CIBERSORT$variable=='T.cell.regulatory..Tregs.')



b1 <-ggplot(data=a1,aes(y=log(value),x=variable,fill=`dd_$group`)) + geom_boxplot()+stat_compare_means(aes(group=`dd_$group`), label = "p.signif",method='t.test')+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+xlab("")+ylab("")+scale_fill_manual(values=c("#E7B800", "#2E9FDF"))+theme(panel.grid=element_blank())+ guides(fill=FALSE)
b2 <-ggplot(data=a2,aes(y=log(value),x=variable,fill=`dd_$group`)) + geom_boxplot()+stat_compare_means(aes(group=`dd_$group`), label = "p.signif",method='t.test')+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+xlab("")+ylab("")+scale_fill_manual(values=c("#E7B800", "#2E9FDF"))+theme(panel.grid=element_blank())+ guides(fill=FALSE)
b3 <-ggplot(data=a5,aes(y=log(value),x=variable,fill=`dd_$group`)) + geom_boxplot()+stat_compare_means(aes(group=`dd_$group`), label = "p.signif",method='t.test')+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+xlab("")+ylab("")+scale_fill_manual(values=c("#E7B800", "#2E9FDF"))+theme(panel.grid=element_blank())+ guides(fill=FALSE)
b4 <-ggplot(data=a3,aes(y=log(value),x=variable,fill=`dd_$group`)) + geom_boxplot()+stat_compare_means(aes(group=`dd_$group`), label = "p.signif",method='t.test')+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+xlab("")+ylab("")+scale_fill_manual(values=c("#E7B800", "#2E9FDF"))+theme(panel.grid=element_blank())+ guides(fill=FALSE)
b5 <-ggplot(data=a4,aes(y=log(value),x=variable,fill=`dd_$group`)) + geom_boxplot()+stat_compare_means(aes(group=`dd_$group`), label = "p.signif",method='t.test')+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+xlab("")+ylab("")+scale_fill_manual(values=c("#E7B800", "#2E9FDF"))+theme(panel.grid=element_blank())+ guides(fill=FALSE)
b6 <-ggplot(data=a6,aes(y=log(value),x=variable,fill=`dd_$group`)) + geom_boxplot()+stat_compare_means(aes(group=`dd_$group`), label = "p.signif",method='t.test')+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+xlab("")+ylab("")+scale_fill_manual(values=c("#E7B800", "#2E9FDF"))+theme(panel.grid=element_blank())+ guides(fill=FALSE)

ggarrange(b1,b2,b3,b4,b5, nrow =1, ncol = 5,align = "h")

dd_CIBERSORT$B.cell.memory_naive <-as.numeric(dd_CIBERSORT$B.cell.memory_CIBERSORT)/(as.numeric(dd_CIBERSORT$B.cell.naive_CIBERSORT)+0.001)
aa1 <- ggplot(data=dd_CIBERSORT,aes(y=log(B.cell.memory_naive),x=dd_$group,fill=`dd_$group`))+stat_compare_means( label = "p.signif",method='t.test') + geom_boxplot()+theme_bw()+xlab("")+ylab("B.cell.memory_naive")+scale_fill_manual(values=c("#E7B800", "#2E9FDF"))+theme(panel.grid=element_blank())+ guides(fill=FALSE)

dd_CIBERSORT$B.cell.plasma_naive <- as.numeric(dd_CIBERSORT$B.cell.plasma_CIBERSORT)/(as.numeric(dd_CIBERSORT$B.cell.naive_CIBERSORT)+0.001)
aa2 <- ggplot(data=dd_CIBERSORT,aes(y=log(B.cell.plasma_naive),x=dd_$group,fill=`dd_$group`)) +stat_compare_means(aes(group=`dd_$group`), label = "p.signif",method='t.test')+ geom_boxplot()+theme_bw()+xlab("")+ylab("B.cell.plasma_naive")+scale_fill_manual(values=c("#E7B800", "#2E9FDF"))+theme(panel.grid=element_blank())+ guides(fill=FALSE)

dd_CIBERSORT$T.cell.CD4..memory <- as.numeric(dd_CIBERSORT$T.cell.CD4..memory.activated_CIBERSORT)/(as.numeric(dd_CIBERSORT$T.cell.CD4..memory.resting_CIBERSORT)+0.001)
aa3 <- ggplot(data=dd_CIBERSORT,aes(y=log(T.cell.CD4..memory),x=dd_$group,fill=`dd_$group`)) +stat_compare_means(aes(group=`dd_$group`), label = "p.signif",method='t.test')+ geom_boxplot()+theme_bw()+xlab("")+ylab("T.cell.CD4..memory")+scale_fill_manual(values=c("#E7B800", "#2E9FDF"))+theme(panel.grid=element_blank())+ guides(fill=FALSE)

dd_CIBERSORT$NK.cell <- as.numeric(dd_CIBERSORT$NK.cell.activated_CIBERSORT)/(as.numeric(dd_CIBERSORT$NK.cell.resting_CIBERSORT)+0.001)
aa4 <- ggplot(data=dd_CIBERSORT,aes(y=log(NK.cell),x=dd_$group,fill=`dd_$group`)) +stat_compare_means(aes(group=`dd_$group`), label = "p.signif",method='t.test')+ geom_boxplot()+theme_bw()+xlab("")+ylab("NK.cell")+scale_fill_manual(values=c("#E7B800", "#2E9FDF"))+theme(panel.grid=element_blank())+ guides(fill=FALSE)

dd_CIBERSORT$Mast.cell <- as.numeric(dd_CIBERSORT$Mast.cell.activated_CIBERSORT)/(as.numeric(dd_CIBERSORT$Mast.cell.resting_CIBERSORT)+0.001)
aa5 <- ggplot(data=dd_CIBERSORT,aes(y=log(Mast.cell),x=dd_$group,fill=`dd_$group`)) +stat_compare_means(aes(group=`dd_$group`), label = "p.signif",method='t.test')+ geom_boxplot()+theme_bw()+xlab("")+ylab("Mast.cell")+scale_fill_manual(values=c("#E7B800", "#2E9FDF"))+theme(panel.grid=element_blank())+ guides(fill=FALSE)



ggplot(data=dd_CIBERSORT,aes(y=value,x=variable,fill=`dd_$group`)) + geom_boxplot()+stat_compare_means(aes(group=`dd_$group`), label = "p.signif",method='t.test')+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+xlab("")+ylab("")+scale_fill_manual(values=c("#E7B800", "#2E9FDF"))+theme(panel.grid=element_blank())

ggplot(data=a1,aes(y=value,x=variable,fill=`dd_$group`)) + geom_boxplot()+stat_compare_means(aes(group=`dd_$group`), label = "p.signif",method='t.test')+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+xlab("")+ylab("")+scale_fill_manual(values=c("#E7B800", "#2E9FDF"))+theme(panel.grid=element_blank())

dd_QUANTISEQ <- cbind(dd_[,51:61],dd_$group)
dd_QUANTISEQ$`dd_$group`<-as.character(dd_QUANTISEQ$`dd_$group`)
dd_QUANTISEQ <- melt(dd_QUANTISEQ)
ggplot(data=dd_QUANTISEQ,aes(y=value,x=variable,fill=`dd_$group`)) + geom_boxplot()+stat_compare_means(aes(group=`dd_$group`), label = "p.signif",method='t.test')+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+xlab("")+ylab("")+scale_fill_manual(values=c("#E7B800", "#2E9FDF"))+theme(panel.grid=element_blank())



##################

data_clinical <- read.table('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\some_file\\shortexon_merge_clinical.txt',sep='\t',row.names=1,head=TRUE)

library(survival)
library(survminer)

da <- filter(data_clinical,data_clinical$cancer_type=='esca_tcga')
A=c()
C=c()
for (i in seq(from=10,to=526)){
  d <- da[order(da[,i]),]
  #print(colnames(d)[i])
  d <- filter(d,as.numeric(d[,i])>=0)
  if (length(d[,1])>=10){
    aa <- round(length(d[,1])/4)
    class <- c(rep('low',aa),rep(NA,length(d[,1])-2*aa),rep('high',aa))
    d$psi_ <- class
    fit <- coxph(Surv(OS_month, OSS_) ~ class, data = d)
    A <- cbind(i,colnames(d)[i],as.numeric(coef(summary(fit))[,5]),as.numeric(round(exp(confint(fit)), 2)[1]),as.numeric(round(exp(confint(fit)), 2)[2]),as.numeric(round(exp(coef(fit)), 2)))
    C<-rbind(C,A)
  }
}

C <- as.data.frame(C,stringsAsFactors = F)
colnames(C) <- c('i','name','p','low','high','HR')
C$p <- as.numeric(C$p)
C$HR<- as.numeric(C$HR)
C$low<- as.numeric(C$low)
C$high<- as.numeric(C$high)
C <- filter(C,C$p<=0.05)

i <- c(38,294,442,481)
da <- da[order(da[,38]),]
temp <- sum(is.na(da[,38])==FALSE)
aa <- round(temp/4)
da$num334 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,294],decreasing = T),]
temp <- sum(is.na(da[,294])==FALSE)
aa <- round(temp/4)
da$num190 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,442],decreasing = T),]
temp <- sum(is.na(da[,442])==FALSE)
aa <- round(temp/4)
da$num326 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,481],decreasing = T),]
temp <- sum(is.na(da[,481])==FALSE)
aa <- round(temp/4)
da$num170 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
m <- apply(da,1,function(x) mean(as.numeric(x[c(527,528,529,530)]),na.rm=T))
da$m4 <-m

i=531
d <- da[order(da[,i]),]
print(colnames(d)[i])
d <- filter(d,as.numeric(d[,i])>=0)
d <- filter(d,as.numeric(d$OS_month)>=0)
aa <- round(length(d[,1])/4)
bb <- round(length(d[,1])/2)
class <- c(rep('low',aa),rep(NA,length(d[,1])-2*aa),rep('high',aa)) 
class3<-c(rep('low',bb),rep('high',length(d[,1])-bb)) 

d$psi_ <- class
d$psi_3<- class3
fit <- survfit(Surv(d$OS_month, OSS) ~ psi_, data = d)
print(fit)
summary(fit)$table
ggsurvplot(fit,
           pval = TRUE, conf.int = FALSE,
           risk.table = FALSE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           # linetype = "strata", # Change line type by groups
           #  surv.median.line = "hv", # Specify median survival
           # ggtheme = theme_bw(legend.position = "right"), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))
an <- ggsurvplot(fit,palette = c("#E7B800", "#2E9FDF"))
an$plot + theme (legend.position = "right")
fit2 <- coxph(Surv(OS_month, OSS_) ~ class3, data = d)
cbind(as.numeric(round(coef(summary(fit2))[,5],10)),as.numeric(round(exp(confint(fit2)), 2)[1]),as.numeric(round(exp(confint(fit2)), 2)[2]),as.numeric(round(exp(coef(fit2)), 2)))


###########
ddd1 <- read.table('D:\\wanglab\\smallexon\\new_60\\06临床样本验证\\MYH11.txt',head=TRUE)
#ggplot(ddd,aes(x=ddd$type,y=PSI,fill=type))+stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+ geom_boxplot(size=0.5,outlier.fill="white",outlier.color="white")+geom_jitter(aes(fill=type),width =0.2,shape = 21,size=2.5)
ggplot(ddd,aes(x=ddd$type,y=PSI,fill=type))+stat_boxplot(geom = "errorbar",width=0.15,color="black")+ geom_boxplot(width=0.5,size=0.5,outlier.fill="white",outlier.color="white",fill=c('darkorange','steelblue1'))+geom_jitter(aes(fill=type),width =0.2,shape = 21,size=2.5)+theme(panel.grid=element_blank()) +xlab('')+scale_fill_manual(values = c('darkorange2','steelblue2'))+theme_classic()
ddd2 <- read.table('D:\\wanglab\\smallexon\\new_60\\06临床样本验证\\RPS24.txt',head=TRUE)
ddd3 <- read.table('D:\\wanglab\\smallexon\\new_60\\06临床样本验证\\TEAD2.txt',head=TRUE)
ddd4 <- read.table('D:\\wanglab\\smallexon\\new_60\\06临床样本验证\\INSR.txt',head=TRUE)
ddd5 <- read.table('D:\\wanglab\\smallexon\\new_60\\06临床样本验证\\ZFAND1.txt',head=TRUE)
ddd6 <- read.table('D:\\wanglab\\smallexon\\new_60\\06临床样本验证\\MBNL1.txt',head=TRUE)

A <- ggplot(ddd1,aes(x=ddd$type,y=PSI,fill=type))+stat_boxplot(geom = "errorbar",width=0.15,color="black")+ geom_boxplot(width=0.5,size=0.5,outlier.fill="white",outlier.color="white",fill=c('darkorange','steelblue1'))+geom_jitter(aes(fill=type),width =0.2,shape = 21,size=2.5)+theme(panel.grid=element_blank()) +xlab('MYH11')+scale_fill_manual(values = c('darkorange2','steelblue2'))+theme_classic()
B <- ggplot(ddd2,aes(x=ddd2$type,y=PSI,fill=type))+stat_boxplot(geom = "errorbar",width=0.15,color="black")+ geom_boxplot(width=0.5,size=0.5,outlier.fill="white",outlier.color="white",fill=c('darkorange','steelblue1'))+geom_jitter(aes(fill=type),width =0.2,shape = 21,size=2.5)+theme(panel.grid=element_blank()) +xlab('RPS24')+scale_fill_manual(values = c('darkorange2','steelblue2'))+theme_classic()
C <- ggplot(ddd3,aes(x=ddd$type,y=PSI,fill=type))+stat_boxplot(geom = "errorbar",width=0.15,color="black")+ geom_boxplot(width=0.5,size=0.5,outlier.fill="white",outlier.color="white",fill=c('darkorange','steelblue1'))+geom_jitter(aes(fill=type),width =0.2,shape = 21,size=2.5)+theme(panel.grid=element_blank()) +xlab('TEAD2')+scale_fill_manual(values = c('darkorange2','steelblue2'))+theme_classic()+ylim(c(0.3,0.6))
D <- ggplot(ddd4,aes(x=ddd$type,y=PSI,fill=type))+stat_boxplot(geom = "errorbar",width=0.15,color="black")+ geom_boxplot(width=0.5,size=0.5,outlier.fill="white",outlier.color="white",fill=c('darkorange','steelblue1'))+geom_jitter(aes(fill=type),width =0.2,shape = 21,size=2.5)+theme(panel.grid=element_blank()) +xlab('INSR')+scale_fill_manual(values = c('darkorange2','steelblue2'))+theme_classic()+ylim(c(0.4,0.8))
E <- ggplot(ddd5,aes(x=ddd$type,y=PSI,fill=type))+stat_boxplot(geom = "errorbar",width=0.15,color="black")+ geom_boxplot(width=0.5,size=0.5,outlier.fill="white",outlier.color="white",fill=c('darkorange','steelblue1'))+geom_jitter(aes(fill=type),width =0.2,shape = 21,size=2.5)+theme(panel.grid=element_blank()) +xlab('ZFAND1')+scale_fill_manual(values = c('darkorange2','steelblue2'))+theme_classic()+ylim(c(0.2,0.8))
F <- ggplot(ddd6,aes(x=ddd6$type,y=PSI,fill=type))+stat_boxplot(geom = "errorbar",width=0.15,color="black")+ geom_boxplot(width=0.5,size=0.5,outlier.fill="white",outlier.color="white",fill=c('darkorange','steelblue1'))+geom_jitter(aes(fill=type),width =0.2,shape = 21,size=2.5)+theme(panel.grid=element_blank()) +xlab('MBNL1')+scale_fill_manual(values = c('darkorange2','steelblue2'))+theme_classic()+ylim(c(0.4,0.8))

ggpubr::ggarrange(A,B,C,D,E,F, nrow =3, ncol = 2)


ddd7 <- read.table('D:\\wanglab\\smallexon\\new_60\\06临床样本验证\\FLNA.txt',head=TRUE)
ddd7 <- ddd7[-c(7,8),]
A <-ggplot(ddd7,aes(x=type,y=PSI,fill=type))+stat_boxplot(geom = "errorbar",width=0.15,color="black")+ geom_boxplot(width=0.5,size=0.5,outlier.fill="white",outlier.color="white",fill=c('darkorange','steelblue1'))+geom_jitter(aes(fill=type),width =0.2,shape = 21,size=2.5)+theme(panel.grid=element_blank()) +xlab('FLNA')+scale_fill_manual(values = c('darkorange2','steelblue2'))+theme_classic()+ylim(c(0.4,0.8))+stat_compare_means(aes(group=type), label = "p.format",method='t.test')
ddd8 <- read.table('D:\\wanglab\\smallexon\\new_60\\06临床样本验证\\ARHGAP17.txt',head=TRUE)
ddd8 <- ddd8[-c(7,8),]
B<-ggplot(ddd8,aes(x=type,y=PSI,fill=type))+stat_boxplot(geom = "errorbar",width=0.15,color="black")+ geom_boxplot(width=0.5,size=0.5,outlier.fill="white",outlier.color="white",fill=c('darkorange','steelblue1'))+geom_jitter(aes(fill=type),width =0.2,shape = 21,size=2.5)+theme(panel.grid=element_blank()) +xlab('ARHGAP17')+scale_fill_manual(values = c('darkorange2','steelblue2'))+theme_classic()+ylim(c(0.4,0.6))+stat_compare_means(aes(group=type), label = "p.format",method='t.test')
ddd9 <- read.table('D:\\wanglab\\smallexon\\new_60\\06临床样本验证\\DST.txt',head=TRUE)
ddd9 <- ddd9[-c(7,8),]
C<- ggplot(ddd9,aes(x=type,y=PSI,fill=type))+stat_boxplot(geom = "errorbar",width=0.15,color="black")+ geom_boxplot(width=0.5,size=0.5,outlier.fill="white",outlier.color="white",fill=c('darkorange','steelblue1'))+geom_jitter(aes(fill=type),width =0.2,shape = 21,size=2.5)+theme(panel.grid=element_blank()) +xlab('DST')+scale_fill_manual(values = c('darkorange2','steelblue2'))+theme_classic()+ylim(c(0.8,1))+stat_compare_means(aes(group=type), label = "p.format",method='t.test')
ddd10 <- read.table('D:\\wanglab\\smallexon\\new_60\\06临床样本验证\\PPP3CC.txt',head=TRUE)
ddd10 <- ddd10[-c(7,8),]
D <- ggplot(ddd10,aes(x=type,y=PSI,fill=type))+stat_boxplot(geom = "errorbar",width=0.15,color="black")+ geom_boxplot(width=0.5,size=0.5,outlier.fill="white",outlier.color="white",fill=c('darkorange','steelblue1'))+geom_jitter(aes(fill=type),width =0.2,shape = 21,size=2.5)+theme(panel.grid=element_blank()) +xlab('PPP3CC')+scale_fill_manual(values = c('darkorange2','steelblue2'))+theme_classic()+ylim(c(0.2,0.7))+stat_compare_means(aes(group=type), label = "p.format",method='t.test')
ddd11 <- read.table('D:\\wanglab\\smallexon\\new_60\\06临床样本验证\\MBNL1_new.txt',head=TRUE)
ddd11 <- ddd11[-c(7,8),]
E <- ggplot(ddd11,aes(x=type,y=PSI,fill=type))+stat_boxplot(geom = "errorbar",width=0.15,color="black")+ geom_boxplot(width=0.5,size=0.5,outlier.fill="white",outlier.color="white",fill=c('darkorange','steelblue1'))+geom_jitter(aes(fill=type),width =0.2,shape = 21,size=2.5)+theme(panel.grid=element_blank()) +xlab('MBNL1')+scale_fill_manual(values = c('darkorange2','steelblue2'))+theme_classic()+ylim(c(0,0.8))+stat_compare_means(aes(group=type), label = "p.format",method='t.test')
ggpubr::ggarrange(A,B,C,D,E,F,nrow =1, ncol = 6)
ddd12 <- read.table('D:\\wanglab\\smallexon\\new_60\\06临床样本验证\\MARK3.txt',head=TRUE)
F <- ggplot(ddd12,aes(x=type,y=PSI,fill=type))+stat_boxplot(geom = "errorbar",width=0.15,color="black")+ geom_boxplot(width=0.5,size=0.5,outlier.fill="white",outlier.color="white",fill=c('darkorange','steelblue1'))+geom_jitter(aes(fill=type),width =0.2,shape = 21,size=2.5)+theme(panel.grid=element_blank()) +xlab('MARK3')+scale_fill_manual(values = c('darkorange2','steelblue2'))+theme_classic()+ylim(c(0.3,0.8))+stat_compare_means(aes(group=type), label = "p.format",method='t.test')

ddd13 <- read.table('D:\\wanglab\\smallexon\\new_60\\06临床样本验证\\TPD52L2_.txt',head=TRUE)
G <- ggplot(ddd13,aes(x=ddd13$type,y=PSI,fill=type))+stat_boxplot(geom = "errorbar",width=0.15,color="black")+ geom_boxplot(width=0.5,size=0.5,outlier.fill="white",outlier.color="white",fill=c('darkorange','steelblue1'))+geom_jitter(aes(fill=type),width =0.2,shape = 21,size=2.5)+theme(panel.grid=element_blank()) +xlab('TPD52L2')+scale_fill_manual(values = c('darkorange2','steelblue2'))+theme_classic()+ylim(c(0.5,1))+stat_compare_means(aes(group=type), label = "p.format",method='t.test')

library(ggpubr)
A <- ggplot(ddd,aes(x=ddd$type,y=PSI,fill=type))+stat_boxplot(geom = "errorbar",width=0.15,color="black")+ geom_boxplot(width=0.5,size=0.5,outlier.fill="white",outlier.color="white",fill=c('darkorange','steelblue1'))+geom_jitter(aes(fill=type),width =0.2,shape = 21,size=2.5)+theme(panel.grid=element_blank()) +xlab('MYH11')+scale_fill_manual(values = c('darkorange2','steelblue2'))+theme_classic()
normal <-c(0.51,0.53,0.51,0.51,0.59,0.52,0.56)
cancer <-c(0.41,0.54,0.54,0.53,0.49,0.48,0.51)
d <- data.frame(normal = normal, cancer = cancer)
ggpaired(d, cond1 = "cancer", cond2 = "normal",
         fill = "condition", palette = "jco")


#####
d <- read.table('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\some_file\\tcga_all.txt',head=TRUE)
d[is.na(d)] <- 0
bk <- c(seq(-1,-0.01,by=0.01),seq(0,1,by=0.01))
pheatmap(d[,7:24],show_rownames = F,color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)))
d1 <- filter(d,d$exon_legth<=60)
d2 <- filter(d,d$exon_legth>60)
d1_ <- d1[,7:24]
d1_[d1_==0]<-NA
d1_$name <- d1$name
d1_2 <- melt(d1_, id.vars = c("name"))
colnames(d1_2)<-c('gene_name','sample',"trv_type")
d1_2$trv_type[d1_2$trv_type!='NA']<-'3_prime_flanking_region'
#genes_to_plot <- c("EGFR","KRAS","MAP2K1","MET","ERBB2","BRAF","RB1","RBM10","TP53","ROS1","ALK","RET")
#waterfall(mut_2, fileType="Custom",plotMutBurden = F, plotGenes=genes_to_plot2,mainPalette=col,rmvSilent = F,geneOrder = genes_to_plot2,variant_class_order=mut_type,mainXlabel=FALSE,mainGrid=F,mainDropMut=F,mainRecurCutoff = 0)
waterfall(d1_2, fileType="MGI",mainGrid = FALSE,plotMutBurden = FALSE,mainXlabel = F)


d2_ <- d2[,7:24]
d2_[d2_==0]<-NA
d2_$name <- d2$name
d2_2 <- melt(d2_, id.vars = c("name"))
colnames(d2_2)<-c('gene','sample',"variant_class")
d2_2$trv_type[d2_2$trv_type!='NA']<-'3_prime_flanking_region'
d2_2 <- filter(d2_2,d2_2$trv_type=='3_prime_flanking_region')
waterfall(d2_2, fileType="Custom",mainGrid = FALSE,plotMutBurden = FALSE,mainXlabel = F)

mut_type <- c("nonsense","frameshift","splice_site","missense","inframe_indel","fusion_gene","3_prime_flanking_region")

sampleOrder=c('GBM','LUSC','KICH','BRCA','UCEC','READ','STAD','ESCA','BLCA','LUAD','KIRC','COAD','HNSC','KIRP','PRAD','LIHC','THCA','PAAD')



###########

dd1 <- read.table('D:\\wanglab\\smallexon\\new_60\\06临床样本验证\\肠癌\\MBNL1.txt',head=TRUE)
dd1$type <-factor(dd1$type,levels=c('T','N'))
#dd1 <- dd1[1:6,]
pp1 <- ggplot(dd1,aes(x=type,y=PSI,fill=type))+stat_boxplot(geom = "errorbar",width=0.15,color="black")+ geom_boxplot(width=0.5,size=0.5,outlier.fill="white",outlier.color="white",fill=c('darkorange','steelblue1'))+geom_jitter(aes(fill=type),width =0.2,shape = 21,size=2.5)+theme(panel.grid=element_blank()) +xlab('MBNL1')+scale_fill_manual(values = c('darkorange2','steelblue2'))+theme_classic()+ylim(c(0.3,0.8))+stat_compare_means(aes(group=type), label = "p.format",method='t.test')

dd2 <- read.table('D:\\wanglab\\smallexon\\new_60\\06临床样本验证\\肠癌\\ARHG.txt',head=TRUE)
dd2$type <-factor(dd2$type,levels=c('T','N'))
pp2 <-ggplot(dd2,aes(x=type,y=PSI,fill=type))+stat_boxplot(geom = "errorbar",width=0.15,color="black")+ geom_boxplot(width=0.5,size=0.5,outlier.fill="white",outlier.color="white",fill=c('darkorange','steelblue1'))+geom_jitter(aes(fill=type),width =0.2,shape = 21,size=2.5)+theme(panel.grid=element_blank()) +xlab('ARHGAP17')+scale_fill_manual(values = c('darkorange2','steelblue2'))+theme_classic()+ylim(c(0.2,0.5))+stat_compare_means(aes(group=type), label = "p.format",method='t.test')

dd3 <- read.table('D:\\wanglab\\smallexon\\new_60\\06临床样本验证\\肠癌\\PPP3CC.txt',head=TRUE)
dd3$type <-factor(dd3$type,levels=c('T','N'))
pp3 <-ggplot(dd3,aes(x=type,y=PSI,fill=type))+stat_boxplot(geom = "errorbar",width=0.15,color="black")+ geom_boxplot(width=0.5,size=0.5,outlier.fill="white",outlier.color="white",fill=c('darkorange','steelblue1'))+geom_jitter(aes(fill=type),width =0.2,shape = 21,size=2.5)+theme(panel.grid=element_blank()) +xlab('PPP3CC')+scale_fill_manual(values = c('darkorange2','steelblue2'))+theme_classic()+ylim(c(0.1,0.6))+stat_compare_means(aes(group=type), label = "p.format",method='t.test')

dd4 <- read.table('D:\\wanglab\\smallexon\\new_60\\06临床样本验证\\肠癌\\TRD.txt',head=TRUE)
dd4$type <-factor(dd4$type,levels=c('T','N'))
pp4 <-ggplot(dd4,aes(x=type,y=PSI,fill=type))+stat_boxplot(geom = "errorbar",width=0.15,color="black")+ geom_boxplot(width=0.5,size=0.5,outlier.fill="white",outlier.color="white",fill=c('darkorange','steelblue1'))+geom_jitter(aes(fill=type),width =0.2,shape = 21,size=2.5)+theme(panel.grid=element_blank()) +xlab('TPD52L2')+scale_fill_manual(values = c('darkorange2','steelblue2'))+theme_classic()+ylim(c(0.2,1))+stat_compare_means(aes(group=type), label = "p.format",method='t.test')

dd5 <- read.table('D:\\wanglab\\smallexon\\new_60\\06临床样本验证\\肠癌\\DST.txt',head=TRUE)
dd5$type <-factor(dd5$type,levels=c('T','N'))
pp5 <-ggplot(dd5,aes(x=type,y=PSI,fill=type))+stat_boxplot(geom = "errorbar",width=0.15,color="black")+ geom_boxplot(width=0.5,size=0.5,outlier.fill="white",outlier.color="white",fill=c('darkorange','steelblue1'))+geom_jitter(aes(fill=type),width =0.2,shape = 21,size=2.5)+theme(panel.grid=element_blank()) +xlab('TPD52L2')+scale_fill_manual(values = c('darkorange2','steelblue2'))+theme_classic()+ylim(c(0.75,0.9))+stat_compare_means(aes(group=type), label = "p.format",method='t.test')

###
a <- c(28.7,20.2,31,22.8,71.3,79.8,69,77.2)
b <- c('short_cancer','shoter_other','long_cancer','long_other','short_cancer','shoter_other','long_cancer','long_other')
c<- c(rep('pol2_bind',4),rep('no_bind',4))
ab <- as.data.frame(cbind(a,b,c))
ab$a <-as.numeric(as.character(ab$a))
ab$b <-factor(ab$b,levels=c('long_other','long_cancer','shoter_other','short_cancer'))



###
d=read.table('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\some_file\\cli_count.txt',head=TRUE)
ggplot(d,aes(x=factor(type, levels=levels(d$type)[c(5,4,3,2,1,6)]),y=number, fill=factor(class,levels=c("short","long"))))+geom_bar(stat='identity',  position=position_dodge())
ggplot(d,aes(x=factor(type, levels=levels(d$type)[c(5,4,3,2,1,6)]),y=number, fill=factor(class,levels=c("short","long"))))+geom_bar(stat='identity',  position=position_dodge(),width=0.7)+scale_fill_manual(values=c("#EC7357","forestgreen"))

d <- filter(d,d$type=='P<0.01' | d$type=='P<0.0001')
