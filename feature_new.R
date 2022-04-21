setwd('D:\\wanglab\\smallexon\\new_60\\01features\\exon_class\\new')
a1 <- read.table('short_cancer.txt')
a2<- read.table('short_other.txt')
a3<- read.table('long_cancer.txt')
a4<- read.table('long_other.txt')
a1$V11 <-'short_cancer'
a2$V11 <-'short_other'
a3$V11 <-'long_cancer'
a4$V11 <-'long_other'
temp <- rbind(a1,a2,a3,a4)
write.table(temp,file='all_exon_class.txt',sep='\t')

exon <- read.table('all_exon_class.txt',head=TRUE,sep='\t')
up <- ggplot(exon,aes(y=exon$upintron_length,x=exon$type)) + geom_boxplot(width=0.6,fill='gray65')+theme_bw()+xlab("")+ylab("upstream intron length")+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+ scale_x_discrete(limits=c("short_cancer", "short_other", "long_cancer","long_other"))+ylim(0,10000)+theme(panel.grid=element_blank())
down <- ggplot(exon,aes(y=exon$downintron_length,x=exon$type)) + geom_boxplot(width=0.6,fill='gray65')+theme_bw()+xlab("")+ylab("downstram intron length")+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+ scale_x_discrete(limits=c("short_cancer", "short_other", "long_cancer","long_other"))+ylim(0,10000)+theme(panel.grid=element_blank())
library(ggpubr)
ggpubr::ggarrange(up,down, nrow = 1, ncol = 2)

##frame
exon$frame <- as.factor(exon$exon_length%%3)
exon$frame <- factor(exon$frame ,levels=c('2','1','0'))
ggplot(exon,aes(y=1,x=exon$type,fill=exon$frame)) +geom_bar(stat = "identity", position = "fill",width=0.6)+theme_bw()+ scale_x_discrete(limits=c("short_cancer", "short_other", "long_cancer","long_other"))+ theme_bw()+scale_fill_manual(values=cc,labels = c("frame2", "frame1", "frame0")) +xlab("")+ylab("")+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),legend.title=element_blank())+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),legend.key.size=unit(0.4,'cm'))+theme(panel.grid=element_blank())

cc <- c('dodgerblue4','dodgerblue1','cadetblue1')
cc <- c('dodgerblue4','dodgerblue4','cadetblue1')


###
setwd('D:\\wanglab\\smallexon\\new_60\\01features\\exon_class\\new')
short_c <- read.table('short_cancer_5ss_score.txt')
short_o <- read.table('short_other_5ss_score.txt')
long_c <- read.table('long_cancer_5ss_score.txt')
long_o <- read.table('long_other_5ss_score.txt')

z <- rbind(short_c,short_o,long_c,long_o)
type <- c(rep('short_cancer',443),rep('short_other',10998),rep('long_cancer',1355),rep('long_other',113959))
z <- as.data.frame(z)
colnames(z)<- c('sequence','5SS_score')
z$type <- type
colnames(z)

#ggplot(z,aes(y=z$`5SS_score`,x=z$type)) + geom_boxplot(width=0.8,fill='moccasin')+theme_bw()+xlab("")+ylab("3' splice site score")+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+ scale_x_discrete(limits=c("short_cancer", "short_other", "long_cancer","long_other"))+theme(panel.grid=element_blank())+ylim(-20,20)

short_c <- read.table('short_cancer_3ss_score.txt')
short_o <- read.table('short_other_3ss_score.txt')
long_c <- read.table('long_cancer_3ss_score.txt')
long_o <- read.table('long_other_3ss_score.txt')

short_o$ss <- str_sub(short_o$V1,start=4L,end=5L)
long_o$ss <- str_sub(long_o$V1,start=4L,end=5L)
long_c$ss <- str_sub(long_c$V1,start=4L,end=5L)
short_c$ss <- str_sub(short_c$V1,start=4L,end=5L)

zz <- rbind(short_c,short_o,long_c,long_o)
type <- c(rep('short_cancer',443),rep('short_other',10998),rep('long_cancer',1355),rep('long_other',113959))
zz <- as.data.frame(zz)
colnames(zz)<- c('sequence','3SS_score')
zz$type <- type
colnames(zz)


ss5 <- ggplot(z,aes(y=z$`3SS_score`,x=z$type)) + geom_boxplot(width=0.8,fill='moccasin')+theme_bw()+xlab("")+ylab("3' splice site score")+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+ scale_x_discrete(limits=c("short_cancer", "short_other", "long_cancer","long_other"))+theme(panel.grid=element_blank())+ylim(-20,30)+stat_compare_means(comparisons = my_comparisons,method='t.test',label.y = c(20,24,28),label = "p.signif",hide.ns = FALSE)
ss3 <- ggplot(zz,aes(y=zz$`5SS_score`,x=zz$type)) + geom_boxplot(width=0.8,fill='moccasin')+theme_bw()+xlab("")+ylab("5' splice site score")+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+ scale_x_discrete(limits=c("short_cancer", "short_other", "long_cancer","long_other"))+theme(panel.grid=element_blank())+ylim(-20,30)+stat_compare_means(comparisons = my_comparisons,method='t.test',label.y = c(20,24,28),label = "p.signif",hide.ns = FALSE)
ggpubr::ggarrange(ss5,ss3, nrow = 1, ncol = 2)

zz$SS3_score <- z$`3SS_score`
zz$SCORE_SUM <- zz$`5SS_score`+zz$SS3_score
ggplot(zz,aes(y=zz$SCORE_SUM,x=zz$type)) + geom_boxplot(width=0.8,fill='moccasin')+theme_bw()+xlab("")+ylab("3ss score + 5ss score")+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+ scale_x_discrete(limits=c("short_cancer", "short_other", "long_cancer","long_other"))+theme(panel.grid=element_blank())

ss5 <- ggplot(z,aes(y=z$`5SS_score`,x=z$type)) + geom_boxplot(width=0.8,fill='moccasin')+theme_bw()+xlab("")+ylab("5' splice site score")+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+ scale_x_discrete(limits=c("short_cancer", "short_other", "long_cancer","long_other"))+theme(panel.grid=element_blank())+ylim(-20,20)
ss3 <- ggplot(zz,aes(y=zz$`3SS_score`,x=zz$type)) + geom_boxplot(width=0.8,fill='moccasin')+theme_bw()+xlab("")+ylab("3' splice site score")+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+ scale_x_discrete(limits=c("short_cancer", "short_other", "long_cancer","long_other"))+theme(panel.grid=element_blank())+ylim(-20,20)
ggpubr::ggarrange(ss5,ss3, nrow = 1, ncol = 2)




###################




d <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\exon_class\\new\\up_everyscore.txt',head=TRUE,row.names = 6)
d$site <- rownames(d)
d<-d[,-1]
library(reshape2)
mydata <- melt(d,id.vars="site",variable.name="type",value.name="score")
mydata$type <- factor(mydata$type,levels=c('short_cancer','short_other','long_cancer','long_other'))
up <- ggplot(data = mydata, mapping = aes(x = as.numeric(site), y = as.numeric(mydata$score), group=type,colour = type)) + geom_line(size=1.1)+theme_bw()+theme(panel.grid=element_blank(),axis.text.x = element_text( angle = 90))+ expand_limits(y=c(0,1))+ scale_x_continuous(breaks=c(-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0))+scale_color_manual(values=c("#EC7357","gold3","forestgreen","#9DC3C1"))+ xlab("") + ylab("") 


d <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\exon_class\\new\\down_everyscore.txt',head=TRUE,row.names = 6)
d$site <- rownames(d)
d<-d[,-1]
mydata2 <- melt(d,id.vars="site",variable.name="type",value.name="score")
mydata2$type <- factor(mydata2$type,levels=c('short_cancer','short_other','long_cancer','long_other'))
down <- ggplot(data = mydata2, mapping = aes(x = as.numeric(site), y = as.numeric(mydata2$score), group=type,colour = type)) + geom_line(size=1.1)+theme_bw()+theme(panel.grid=element_blank(),axis.text.x = element_text( angle = 90))+ expand_limits(y=c(0,1))+ scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100))+scale_color_manual(values=c("#EC7357","gold3","forestgreen","#9DC3C1"))+ xlab("") + ylab("") 


ggpubr::ggarrange(up,down, nrow = 1, ncol = 2)


# 
temp <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\exon_class\\new\\all_exon_class.txt',head=TRUE)

