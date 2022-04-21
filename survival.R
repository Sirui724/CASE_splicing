
data_clinical <- read.table('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\some_file\\shortexon_merge_clinical.txt',sep='\t',row.names=1,head=TRUE)

library(survival)
library(survminer)

da <- filter(data_clinical,data_clinical$cancer_type=='blca_tcga_pan_can_atlas_2018')
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

i <- c(334,190,326,170)
da <- da[order(da[,334]),]
temp <- sum(is.na(da[,334])==FALSE)
aa <- round(temp/4)
da$num334 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,190]),]
temp <- sum(is.na(da[,190])==FALSE)
aa <- round(temp/4)
da$num190 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,326]),]
temp <- sum(is.na(da[,326])==FALSE)
aa <- round(temp/4)
da$num326 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,170],decreasing = T),]
temp <- sum(is.na(da[,170])==FALSE)
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

da <- filter(data_clinical,data_clinical$cancer_type=='brca_tcga')
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
da$num1 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,294],decreasing = T),]
temp <- sum(is.na(da[,294])==FALSE)
aa <- round(temp/4)
da$num2 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,442],decreasing = T),]
temp <- sum(is.na(da[,442])==FALSE)
aa <- round(temp/4)
da$num3 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,481],decreasing = T),]
temp <- sum(is.na(da[,481])==FALSE)
aa <- round(temp/4)
da$num4 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
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
fit2 <- coxph(Surv(OS_month, OSS_) ~ class, data = d)
cbind(as.numeric(round(coef(summary(fit2))[,5],10)),as.numeric(round(exp(confint(fit2)), 2)[1]),as.numeric(round(exp(confint(fit2)), 2)[2]),as.numeric(round(exp(coef(fit2)), 2)))


###

da <- filter(data_clinical,data_clinical$cancer_type=='coadread_tcga')
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

i <- c(123,349,276,145)
da <- da[order(da[,123],decreasing = T),]
temp <- sum(is.na(da[,123])==FALSE)
aa <- round(temp/4)
da$num1 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,349],decreasing = F),]
temp <- sum(is.na(da[,349])==FALSE)
aa <- round(temp/4)
da$num2 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,276],decreasing = F),]
temp <- sum(is.na(da[,276])==FALSE)
aa <- round(temp/4)
da$num3 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,145],decreasing = T),]
temp <- sum(is.na(da[,145])==FALSE)
aa <- round(temp/4)
da$num4 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
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
           palette = c( "#2E9FDF","#E7B800"))
an <- ggsurvplot(fit,palette = c("#E7B800", "#2E9FDF"))
an$plot + theme (legend.position = "right")
fit2 <- coxph(Surv(OS_month, OSS_) ~ class, data = d)
cbind(as.numeric(round(coef(summary(fit2))[,5],10)),as.numeric(round(exp(confint(fit2)), 2)[1]),as.numeric(round(exp(confint(fit2)), 2)[2]),as.numeric(round(exp(coef(fit2)), 2)))


####
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

i <- c(17,170,295,516)
da <- da[order(da[,17],decreasing = T),]
temp <- sum(is.na(da[,17])==FALSE)
aa <- round(temp/4)
da$num1 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,170],decreasing = T),]
temp <- sum(is.na(da[,170])==FALSE)
aa <- round(temp/4)
da$num2 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,295],decreasing = F),]
temp <- sum(is.na(da[,295])==FALSE)
aa <- round(temp/4)
da$num3 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,516],decreasing = T),]
temp <- sum(is.na(da[,516])==FALSE)
aa <- round(temp/4)
da$num4 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
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
           palette = c( "#E7B800","#2E9FDF"))
an <- ggsurvplot(fit,palette = c("#E7B800", "#2E9FDF"))
an$plot + theme (legend.position = "right")
fit2 <- coxph(Surv(OS_month, OSS_) ~ class, data = d)
cbind(as.numeric(round(coef(summary(fit2))[,5],10)),as.numeric(round(exp(confint(fit2)), 2)[1]),as.numeric(round(exp(confint(fit2)), 2)[2]),as.numeric(round(exp(coef(fit2)), 2)))


###
da <- filter(data_clinical,data_clinical$cancer_type=='gbm_tcga')
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

i <- c(128,143,449,425)
da <- da[order(da[,128],decreasing = T),]
temp <- sum(is.na(da[,128])==FALSE)
aa <- round(temp/4)
da$num1 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,143],decreasing = F),]
temp <- sum(is.na(da[,143])==FALSE)
aa <- round(temp/4)
da$num2 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,449],decreasing = F),]
temp <- sum(is.na(da[,449])==FALSE)
aa <- round(temp/4)
da$num3 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,425],decreasing = T),]
temp <- sum(is.na(da[,425])==FALSE)
aa <- round(temp/4)
da$num4 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
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
           palette = c( "#E7B800","#2E9FDF"))
an <- ggsurvplot(fit,palette = c("#E7B800", "#2E9FDF"))
an$plot + theme (legend.position = "right")
fit2 <- coxph(Surv(OS_month, OSS_) ~ class, data = d)
cbind(as.numeric(round(coef(summary(fit2))[,5],10)),as.numeric(round(exp(confint(fit2)), 2)[1]),as.numeric(round(exp(confint(fit2)), 2)[2]),as.numeric(round(exp(coef(fit2)), 2)))


###
da <- filter(data_clinical,data_clinical$cancer_type=='hnsc_tcga')
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

i <- c(185,329,438,494)
da <- da[order(da[,185],decreasing = F),]
temp <- sum(is.na(da[,185])==FALSE)
aa <- round(temp/4)
da$num1 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,329],decreasing = F),]
temp <- sum(is.na(da[,329])==FALSE)
aa <- round(temp/4)
da$num2 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,438],decreasing = T),]
temp <- sum(is.na(da[,438])==FALSE)
aa <- round(temp/4)
da$num3 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,494],decreasing = T),]
temp <- sum(is.na(da[,494])==FALSE)
aa <- round(temp/4)
da$num4 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
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
           palette = c( "#E7B800","#2E9FDF"))
an <- ggsurvplot(fit,palette = c("#E7B800", "#2E9FDF"))
an$plot + theme (legend.position = "right")
fit2 <- coxph(Surv(OS_month, OSS_) ~ class, data = d)
cbind(as.numeric(round(coef(summary(fit2))[,5],10)),as.numeric(round(exp(confint(fit2)), 2)[1]),as.numeric(round(exp(confint(fit2)), 2)[2]),as.numeric(round(exp(coef(fit2)), 2)))


###
da <- filter(data_clinical,data_clinical$cancer_type=='kich_tcga')
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

i <- c(194,204,508,498)
da <- da[order(da[,194],decreasing = T),]
temp <- sum(is.na(da[,194])==FALSE)
aa <- round(temp/4)
da$num1 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,204],decreasing = T),]
temp <- sum(is.na(da[,204])==FALSE)
aa <- round(temp/4)
da$num2 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,508],decreasing = T),]
temp <- sum(is.na(da[,508])==FALSE)
aa <- round(temp/4)
da$num3 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,498],decreasing = F),]
temp <- sum(is.na(da[,498])==FALSE)
aa <- round(temp/4)
da$num4 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
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
           palette = c( "#E7B800","#2E9FDF"))
an <- ggsurvplot(fit,palette = c("#E7B800", "#2E9FDF"))
an$plot + theme (legend.position = "right")
fit2 <- coxph(Surv(OS_month, OSS_) ~ class, data = d)
cbind(as.numeric(round(coef(summary(fit2))[,5],10)),as.numeric(round(exp(confint(fit2)), 2)[1]),as.numeric(round(exp(confint(fit2)), 2)[2]),as.numeric(round(exp(coef(fit2)), 2)))

###





da <- filter(data_clinical,data_clinical$cancer_type=='kirc_tcga')
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

i <- c(382,17,12,214)
da <- da[order(da[,382],decreasing = T),]
temp <- sum(is.na(da[,382])==FALSE)
aa <- round(temp/4)
da$num1 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,17],decreasing = T),]
temp <- sum(is.na(da[,17])==FALSE)
aa <- round(temp/4)
da$num2 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,12],decreasing = F),]
temp <- sum(is.na(da[,12])==FALSE)
aa <- round(temp/4)
da$num3 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,214],decreasing = T),]
temp <- sum(is.na(da[,214])==FALSE)
aa <- round(temp/4)
da$num4 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
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
           palette = c( "#E7B800","#2E9FDF"))
an <- ggsurvplot(fit,palette = c("#E7B800", "#2E9FDF"))
an$plot + theme (legend.position = "right")
fit2 <- coxph(Surv(OS_month, OSS_) ~ class, data = d)
cbind(as.numeric(round(coef(summary(fit2))[,5],10)),as.numeric(round(exp(confint(fit2)), 2)[1]),as.numeric(round(exp(confint(fit2)), 2)[2]),as.numeric(round(exp(coef(fit2)), 2)))


###
da <- filter(data_clinical,data_clinical$cancer_type=='kirp_tcga')
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

i <- c(381,88,373,431)
da <- da[order(da[,381],decreasing = T),]
temp <- sum(is.na(da[,381])==FALSE)
aa <- round(temp/4)
da$num1 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,88],decreasing = T),]
temp <- sum(is.na(da[,88])==FALSE)
aa <- round(temp/4)
da$num2 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,373],decreasing = T),]
temp <- sum(is.na(da[,373])==FALSE)
aa <- round(temp/4)
da$num3 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,431],decreasing = T),]
temp <- sum(is.na(da[,431])==FALSE)
aa <- round(temp/4)
da$num4 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
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
           palette = c( "#E7B800","#2E9FDF"))
an <- ggsurvplot(fit,palette = c("#E7B800", "#2E9FDF"))
an$plot + theme (legend.position = "right")
fit2 <- coxph(Surv(OS_month, OSS_) ~ class, data = d)
cbind(as.numeric(round(coef(summary(fit2))[,5],10)),as.numeric(round(exp(confint(fit2)), 2)[1]),as.numeric(round(exp(confint(fit2)), 2)[2]),as.numeric(round(exp(coef(fit2)), 2)))


####
da <- filter(data_clinical,data_clinical$cancer_type=='lihc_tcga')
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

i <- c(419,179,62,234)
da <- da[order(da[,419],decreasing = T),]
temp <- sum(is.na(da[,419])==FALSE)
aa <- round(temp/4)
da$num1 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,179],decreasing = T),]
temp <- sum(is.na(da[,179])==FALSE)
aa <- round(temp/4)
da$num2 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,62],decreasing = T),]
temp <- sum(is.na(da[,62])==FALSE)
aa <- round(temp/4)
da$num3 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,234],decreasing = T),]
temp <- sum(is.na(da[,234])==FALSE)
aa <- round(temp/4)
da$num4 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
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
           palette = c( "#E7B800","#2E9FDF"))
an <- ggsurvplot(fit,palette = c("#E7B800", "#2E9FDF"))
an$plot + theme (legend.position = "right")
fit2 <- coxph(Surv(OS_month, OSS_) ~ class, data = d)
cbind(as.numeric(round(coef(summary(fit2))[,5],10)),as.numeric(round(exp(confint(fit2)), 2)[1]),as.numeric(round(exp(confint(fit2)), 2)[2]),as.numeric(round(exp(coef(fit2)), 2)))



###
da <- filter(data_clinical,data_clinical$cancer_type=='luad_tcga')
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

i <- c(211,461,72,349)
da <- da[order(da[,211],decreasing = T),]
temp <- sum(is.na(da[,211])==FALSE)
aa <- round(temp/4)
da$num1 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,461],decreasing = F),]
temp <- sum(is.na(da[,461])==FALSE)
aa <- round(temp/4)
da$num2 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,72],decreasing = T),]
temp <- sum(is.na(da[,72])==FALSE)
aa <- round(temp/4)
da$num3 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,349],decreasing = F),]
temp <- sum(is.na(da[,349])==FALSE)
aa <- round(temp/4)
da$num4 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
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
           palette = c( "#E7B800","#2E9FDF"))
an <- ggsurvplot(fit,palette = c("#E7B800", "#2E9FDF"))
an$plot + theme (legend.position = "right")
fit2 <- coxph(Surv(OS_month, OSS_) ~ class, data = d)
cbind(as.numeric(round(coef(summary(fit2))[,5],10)),as.numeric(round(exp(confint(fit2)), 2)[1]),as.numeric(round(exp(confint(fit2)), 2)[2]),as.numeric(round(exp(coef(fit2)), 2)))

##
da <- filter(data_clinical,data_clinical$cancer_type=='lusc_tcga')
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

i <- c(74,302,36,167)
da <- da[order(da[,74],decreasing = F),]
temp <- sum(is.na(da[,74])==FALSE)
aa <- round(temp/4)
da$num1 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,302],decreasing = F),]
temp <- sum(is.na(da[,302])==FALSE)
aa <- round(temp/4)
da$num2 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,36],decreasing = T),]
temp <- sum(is.na(da[,36])==FALSE)
aa <- round(temp/4)
da$num3 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,167],decreasing = T),]
temp <- sum(is.na(da[,167])==FALSE)
aa <- round(temp/4)
da$num4 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
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
           palette = c( "#E7B800","#2E9FDF"))
an <- ggsurvplot(fit,palette = c("#E7B800", "#2E9FDF"))
an$plot + theme (legend.position = "right")
fit2 <- coxph(Surv(OS_month, OSS_) ~ class, data = d)
cbind(as.numeric(round(coef(summary(fit2))[,5],10)),as.numeric(round(exp(confint(fit2)), 2)[1]),as.numeric(round(exp(confint(fit2)), 2)[2]),as.numeric(round(exp(coef(fit2)), 2)))


###
da <- filter(data_clinical,data_clinical$cancer_type=='paad_tcga')
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

i <- c(285,430,349,236)
da <- da[order(da[,285],decreasing = T),]
temp <- sum(is.na(da[,285])==FALSE)
aa <- round(temp/4)
da$num1 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,430],decreasing = T),]
temp <- sum(is.na(da[,430])==FALSE)
aa <- round(temp/4)
da$num2 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,349],decreasing = F),]
temp <- sum(is.na(da[,349])==FALSE)
aa <- round(temp/4)
da$num3 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,236],decreasing = F),]
temp <- sum(is.na(da[,236])==FALSE)
aa <- round(temp/4)
da$num4 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
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
           palette = c( "#E7B800","#2E9FDF"))
an <- ggsurvplot(fit,palette = c("#E7B800", "#2E9FDF"))
an$plot + theme (legend.position = "right")
fit2 <- coxph(Surv(OS_month, OSS_) ~ class, data = d)
cbind(as.numeric(round(coef(summary(fit2))[,5],10)),as.numeric(round(exp(confint(fit2)), 2)[1]),as.numeric(round(exp(confint(fit2)), 2)[2]),as.numeric(round(exp(coef(fit2)), 2)))


##
da <- filter(data_clinical,data_clinical$cancer_type=='stad_tcga')
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

i <- c(89,237,127,430)
da <- da[order(da[,89],decreasing = F),]
temp <- sum(is.na(da[,89])==FALSE)
aa <- round(temp/4)
da$num1 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,237],decreasing = F),]
temp <- sum(is.na(da[,237])==FALSE)
aa <- round(temp/4)
da$num2 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,127],decreasing = T),]
temp <- sum(is.na(da[,127])==FALSE)
aa <- round(temp/4)
da$num3 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,430],decreasing = F),]
temp <- sum(is.na(da[,430])==FALSE)
aa <- round(temp/4)
da$num4 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
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
           palette = c( "#E7B800","#2E9FDF"))
an <- ggsurvplot(fit,palette = c("#E7B800", "#2E9FDF"))
an$plot + theme (legend.position = "right")
fit2 <- coxph(Surv(OS_month, OSS_) ~ class, data = d)
cbind(as.numeric(round(coef(summary(fit2))[,5],10)),as.numeric(round(exp(confint(fit2)), 2)[1]),as.numeric(round(exp(confint(fit2)), 2)[2]),as.numeric(round(exp(coef(fit2)), 2)))

###
da <- filter(data_clinical,data_clinical$cancer_type=='thca_tcga')
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

i <- c(238,127,345,391)
da <- da[order(da[,238],decreasing = T),]
temp <- sum(is.na(da[,238])==FALSE)
aa <- round(temp/4)
da$num1 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,127],decreasing = F),]
temp <- sum(is.na(da[,127])==FALSE)
aa <- round(temp/4)
da$num2 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,345],decreasing = T),]
temp <- sum(is.na(da[,345])==FALSE)
aa <- round(temp/4)
da$num3 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,391],decreasing = F),]
temp <- sum(is.na(da[,391])==FALSE)
aa <- round(temp/4)
da$num4 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
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
           palette = c( "#E7B800","#2E9FDF"))
an <- ggsurvplot(fit,palette = c("#E7B800", "#2E9FDF"))
an$plot + theme (legend.position = "right")
fit2 <- coxph(Surv(OS_month, OSS_) ~ class, data = d)
cbind(as.numeric(round(coef(summary(fit2))[,5],10)),as.numeric(round(exp(confint(fit2)), 2)[1]),as.numeric(round(exp(confint(fit2)), 2)[2]),as.numeric(round(exp(coef(fit2)), 2)))


##
da <- filter(data_clinical,data_clinical$cancer_type=='ucec_tcga')
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

i <- c(236,431,430,137)
da <- da[order(da[,236],decreasing = F),]
temp <- sum(is.na(da[,236])==FALSE)
aa <- round(temp/4)
da$num1 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,431],decreasing = F),]
temp <- sum(is.na(da[,431])==FALSE)
aa <- round(temp/4)
da$num2 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,430],decreasing = F),]
temp <- sum(is.na(da[,430])==FALSE)
aa <- round(temp/4)
da$num3 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
da <- da[order(da[,137],decreasing = T),]
temp <- sum(is.na(da[,137])==FALSE)
aa <- round(temp/4)
da$num4 <- c(rep(1,aa),rep(2,aa),rep(3,aa),rep(4,(temp-3*aa)),rep(NA,length(da[,1])-temp))
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
           palette = c( "#E7B800","#2E9FDF"))
an <- ggsurvplot(fit,palette = c("#E7B800", "#2E9FDF"))
an$plot + theme (legend.position = "right")
fit2 <- coxph(Surv(OS_month, OSS_) ~ class, data = d)
cbind(as.numeric(round(coef(summary(fit2))[,5],10)),as.numeric(round(exp(confint(fit2)), 2)[1]),as.numeric(round(exp(confint(fit2)), 2)[2]),as.numeric(round(exp(coef(fit2)), 2)))
