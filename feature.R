

#####motif
seq <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\meme_out\\motif1.txt',head=FALSE)

library("ggseqlogo")
ggplot()+geom_logo(as.character(seq$V6))+theme_classic()


seq <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\frequency.txt',head=FALSE)
library(dplyr)
seq <- filter(seq,seq$V5>=4)
ggplot()+geom_logo(as.character(seq$V1))+theme_classic()


seq <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\exon_1.txt',head=FALSE)
library(dplyr)
m1 <- ggplot()+geom_logo(as.character(seq$V1))+theme_classic()
seq <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\exon_2.txt',head=FALSE)
m2 <- ggplot()+geom_logo(as.character(seq$V1))+theme_classic()
ggpubr::ggarrange(m1,m2, nrow =2, ncol = 1)
seq <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\exon_motif.txt',head=FALSE)
seq1 <- seq[1:9,]
seq2 <- seq[10:14,]
m1 <- ggplot()+geom_logo(as.character(seq1$V2))+theme_classic()
m2 <- ggplot()+geom_logo(as.character(seq2$V2))+theme_classic()
ggpubr::ggarrange(m1,m2, nrow =2, ncol = 1)




seq <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\up_motif.txt',head=FALSE)
seq1 <- seq[1:18,]
seq2 <- seq[21:35,]
m1 <- ggplot()+geom_logo(as.character(seq1$V2))+theme_classic()+ scale_x_discrete(limits=c('','1','2','3','4','5','6','7',''))
m2 <- ggplot()+geom_logo(as.character(seq2$V2))+theme_classic()+ scale_x_discrete(limits=c('','1','2','3','4','5','6','7',''))
ggpubr::ggarrange(m1,m2, nrow =2, ncol = 1)

seq <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\down_motif.txt',head=FALSE)
seq1 <- seq[1:36,]
seq2 <- seq[37:46,]
m1 <- ggplot()+geom_logo(as.character(seq1$V2))+theme_classic()+ coord_cartesian(xlim = c(4,11))
m2 <- ggplot()+geom_logo(as.character(seq2$V2))+theme_classic()+ coord_cartesian(xlim = c(2,10))
ggpubr::ggarrange(m1,m2, nrow =2, ncol = 1)




seq <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\up300_motif.txt',head=FALSE)
seq1 <- seq[1:7,]
seq2 <- seq[8:43,]
m1 <- ggplot()+geom_logo(as.character(seq1$V2))+theme_classic()+ scale_x_discrete(limits=c('','','1','2','3','4','5','6','7','8'))
m2 <- ggplot()+geom_logo(as.character(seq2$V2))+theme_classic()+ scale_x_discrete(limits=c('','1','2','3','4','5','6','7','8','9'))
ggpubr::ggarrange(m1,m2, nrow =2, ncol = 1)

seq <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\down300_motif.txt',head=FALSE)
seq1 <- seq[1:6,]
seq2 <- seq[7:31,]
m1 <- ggplot()+geom_logo(as.character(seq1$V2))+theme_classic()+ scale_x_discrete(limits=c('','1','2','3','4','5','6','7',''))
m2 <- ggplot()+geom_logo(as.character(seq2$V2))+theme_classic()+ scale_x_discrete(limits=c('','','1','2','3','4','5','6','7','8',''))
ggpubr::ggarrange(m1,m2, nrow =2, ncol = 1)

seq1 <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\down100_sl.txt',head=FALSE)
a <- ggplot()+geom_logo(as.character(seq1$V2))+theme_classic()+ggtitle('down100_long')
seq2 <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\down_300_sl_motif.txt',head=FALSE)
b<- ggplot()+geom_logo(as.character(seq2$V2[2:7]))+theme_classic()+ggtitle('down300_long')
seq3 <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\down_500_sl_motif.txt',head=FALSE)
c <- ggplot()+geom_logo(as.character(seq3$V2[1:3]))+theme_classic()+ggtitle('down300_long_m1')
d <- ggplot()+geom_logo(as.character(seq3$V2[4:6]))+theme_classic()+ggtitle('down300_long_m2')
ggpubr::ggarrange(a,b,c,d, nrow =1, ncol = 4)


seq4 <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\up100_sl_motif.txt',head=FALSE)
e <- ggplot()+geom_logo(as.character(seq4$V2))+theme_classic()+ggtitle('up100_long')
seq5 <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\up300_motif.txt',head=FALSE)
f <- ggplot()+geom_logo(as.character(seq5$V2))+theme_classic()+ggtitle('up300_long')
seq6 <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\up500_sl_motif.txt',head=FALSE)
g <- ggplot()+geom_logo(as.character(seq6$V2))+theme_classic()+ggtitle('up500_long')



seq1 <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\down100_other_motif.txt',head=FALSE)
a <- ggplot()+geom_logo(as.character(seq1$V2[1:10]))+theme_classic()+ggtitle('down100_other_m1')
b<- ggplot()+geom_logo(as.character(seq1$V2[11:16]))+theme_classic()+ggtitle('down100_other_m2')
seq2 <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\down300_other_motif.txt',head=FALSE)
c <- ggplot()+geom_logo(as.character(seq2$V2[1:6]))+theme_classic()+ggtitle('down300_other_m1')
d <- ggplot()+geom_logo(as.character(seq2$V2[c(7,8,9,10,11,12,13,14,15,22,23,24,25,26,27,28,29)]))+theme_classic()+ggtitle('down300_other_m2')
e <- ggplot()+geom_logo(as.character(seq2$V2[16:21]))+theme_classic()+ggtitle('down300_other_m3')

seq3 <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\down500_other_motif.txt',head=FALSE)
g <- ggplot()+geom_logo(as.character(seq3$V2[1:30]))+theme_classic()+ggtitle('down500_other_m1')
h <- ggplot()+geom_logo(as.character(seq3$V2[31:40]))+theme_classic()+ggtitle('down500_other_m2')

ggpubr::ggarrange(a,b,c,d,e,g,h, nrow =2, ncol = 4)



seq1 <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\up100_other_motif.txt',head=FALSE)
a <- ggplot()+geom_logo(as.character(seq1$V2[31:58]))+theme_classic()+ggtitle('up100_other_m1')
b<- ggplot()+geom_logo(as.character(seq1$V2[9:30]))+theme_classic()+ggtitle('up100_other_m2')
c<- ggplot()+geom_logo(as.character(seq1$V2[1:8]))+theme_classic()+ggtitle('up100_other_m3')
seq2 <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\up300_other_motif.txt',head=FALSE)
d <- ggplot()+geom_logo(as.character(seq2$V2[29:54]))+theme_classic()+ggtitle('up300_other_m1')
e<- ggplot()+geom_logo(as.character(seq2$V2[14:28]))+theme_classic()+ggtitle('up300_other_m2')
f<- ggplot()+geom_logo(as.character(seq2$V2[1:13]))+theme_classic()+ggtitle('up300_other_m3')
seq3<- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\up500_other_motif.txt',head=FALSE)
g <- ggplot()+geom_logo(as.character(seq3$V2[1:19]))+theme_classic()+ggtitle('up500_other_m1')
h<- ggplot()+geom_logo(as.character(seq3$V2[20:40]))+theme_classic()+ggtitle('up500_other_m2')

ggpubr::ggarrange(a,b,c,d,e,f,g,h, nrow =2, ncol = 4)



###
seq1 <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\up300_',head=FALSE)
a <- ggplot()+geom_logo(as.character(seq1$V2[28:44]))+theme_classic()+ggtitle('up300_cs_vs_so_m1')
b<- ggplot()+geom_logo(as.character(seq1$V2[1:27]))+theme_classic()+ggtitle('up300_cs_vs_so_m2')

seq2<- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\sc_so_motif.txt',head=FALSE)
c <- ggplot()+geom_logo(as.character(seq2$V2))+theme_classic()+ggtitle('cs_vs_so_m1')

seq3<- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\down300_sc_so.txt',head=FALSE)
d<- ggplot()+geom_logo(as.character(seq3$V2[1:6]))+theme_classic()+ggtitle('down300_sc_so_m1')
e<- ggplot()+geom_logo(as.character(seq3$V2[7:19]))+theme_classic()+ggtitle('down300_sc_so_m2')
f<- ggplot()+geom_logo(as.character(seq3$V2[20:25]))+theme_classic()+ggtitle('down300_sc_so_m3')

ggpubr::ggarrange(a,b,c,d,e,f, nrow =2, ncol = 3)



seq1 <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\up300_lc_lo.txt',head=FALSE)
a <- ggplot()+geom_logo(as.character(seq1$V2[32:74]))+theme_classic()+ggtitle('up300_lc_vs_lo_m1')
b<- ggplot()+geom_logo(as.character(seq1$V2[1:31]))+theme_classic()+ggtitle('up300_lc_vs_lo_m2')

seq2<- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\lc_lo.txt',head=FALSE)
c <- ggplot()+geom_logo(as.character(seq2$V2[1:25]))+theme_classic()+ggtitle('lc_vs_lo_m1')
d <- ggplot()+geom_logo(as.character(seq2$V2[26:37]))+theme_classic()+ggtitle('lc_vs_lo_m2')

seq3<- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\down300_lc_lo.txt',head=FALSE)
e<- ggplot()+geom_logo(as.character(seq3$V2[1:5]))+theme_classic()+ggtitle('down300_lc_lo_m1')
f<- ggplot()+geom_logo(as.character(seq3$V2[11:17]))+theme_classic()+ggtitle('down300_lc_lo_m2')

ggpubr::ggarrange(a,b,c,d,e,f, nrow =2, ncol = 3)
####

seq1 <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\up300_sc_so.txt',head=FALSE)
a <- ggplot()+geom_logo(as.character(seq1$V2[1:7]))+theme_classic()+ggtitle('up300_cs_vs_slother_m1')
b<- ggplot()+geom_logo(as.character(seq1$V2[8:24]))+theme_classic()+ggtitle('up300_cs_vs_slother_m2')
c<- ggplot()+geom_logo(as.character(seq1$V2[25:54]))+theme_classic()+ggtitle('up300_cs_vs_slother_m3')






seq<- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\exon_sc_vs_lc.txt',head=FALSE)
exon<- ggplot()+geom_logo(as.character(seq$V2))+theme_classic()+ggtitle('exon_sc_vs_lc')


d <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\up300_other_motifsite.txt',head=FALSE,row.names = 1)
colnames(d)<-c('motif3','motif2','motif1')
d<- d[c("motif1", "motif2", "motif3")]
pheatmap(t(d),cluster_rows = F,cluster_cols = F,color = c(colorRampPalette(colors = c("white","firebrick3"))(500)),show_colnames = F)

d <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\down300_other_motifsite.txt',head=FALSE,row.names = 1)
colnames(d)<-c('motif1','motif2','motif3')
pheatmap(t(d),cluster_rows = F,cluster_cols = F,color = c(colorRampPalette(colors = c("white","firebrick3"))(500)),show_colnames = F)

d <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\exon_motifsite.txt',head=FALSE,row.names = 1)
colnames(d)<-c('motif1','motif2')
pheatmap(t(d),cluster_rows = F,cluster_cols = F,color = c(colorRampPalette(colors = c("white","skyblue"))(500)),show_colnames = F)





#######

seq1 <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\up300_so_lo.txt',head=FALSE)
a <- ggplot()+geom_logo(as.character(seq1$V2[1:67]))+theme_classic()+ggtitle('up300_so_lo_m1')
b <- ggplot()+geom_logo(as.character(seq1$V2[68:130]))+theme_classic()+ggtitle('up300_so_lo_m2')

seq2 <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\down300_so_lo.txt',head=FALSE)
ggplot()+geom_logo(as.character(seq2$V2[1:23]))+theme_classic()+ggtitle('up300_so_lo_m1')
c <- ggplot()+geom_logo(as.character(seq2$V2[24:56]))+theme_classic()+ggtitle('up300_so_lo_m1')
d <- ggplot()+geom_logo(as.character(seq2$V2[57:133]))+theme_classic()+ggtitle('up300_so_lo_m1')

seq3 <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\so_lo_exon.txt',head=FALSE)
e <- ggplot()+geom_logo(as.character(seq3$V2[1:21]))+theme_classic()+ggtitle('exonso_lo_m1')
ggplot()+geom_logo(as.character(seq3$V2[44:113]))+theme_classic()+ggtitle('up300_so_lo_m2')
ggplot()+geom_logo(as.character(seq3$V2[114:296]))+theme_classic()+ggtitle('up300_so_lo_m3')
f<-ggplot()+geom_logo(as.character(seq3$V2[297:361]))+theme_classic()+ggtitle('up300_so_lo_m4')


##############
exon <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\exon_class\\allexon_length_class.txt',head=TRUE)
up <- ggplot(exon,aes(y=exon$upintron_length,x=exon$type)) + geom_boxplot(width=0.6,fill='gray65')+theme_bw()+xlab("")+ylab("upstream intron length")+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+ scale_x_discrete(limits=c("short_cancer", "short_other", "long_cancer","long_other"))+ylim(0,10000)+theme(panel.grid=element_blank())
down <- ggplot(exon,aes(y=exon$downintron_length,x=exon$type)) + geom_boxplot(width=0.6,fill='gray65')+theme_bw()+xlab("")+ylab("downstram intron length")+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+ scale_x_discrete(limits=c("short_cancer", "short_other", "long_cancer","long_other"))+ylim(0,10000)+theme(panel.grid=element_blank())
#library(ggpubr)
ggpubr::ggarrange(up,down, nrow = 1, ncol = 2)
t1<- filter(exon,exon$type=='short_cancer')
t2<- filter(exon,exon$type=='short_other')
t3 <-filter(exon,exon$type=='long_cancer')
t4 <-filter(exon,exon$type=='long_other')

exon$frame <- as.factor(exon$exon_length%%3)
exon$frame <- factor(exon$frame ,levels=c('2','1','0'))
ggplot(exon,aes(y=1,x=exon$type,fill=exon$frame)) +geom_bar(stat = "identity", position = "fill",width=0.6)+theme_bw()+ scale_x_discrete(limits=c("short_cancer", "short_other", "long_cancer","long_other"))+ theme_bw()+scale_fill_manual(values=cc,labels = c("frame2", "frame1", "frame0")) +xlab("")+ylab("")+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),legend.title=element_blank())+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),legend.key.size=unit(0.4,'cm'))+theme(panel.grid=element_blank())

cc <- c('dodgerblue4','dodgerblue1','cadetblue1')
cc <- c('dodgerblue4','dodgerblue4','cadetblue1')

#############splice site
setwd('D:\\wanglab\\smallexon\\new_60\\01features\\014splicesite')
short_c <- read.table('short_cancer_3ss_score.txt')
short_o <- read.table('short_other_3ss_score.txt')
long_c <- read.table('long_cancer_3ss_score.txt')
long_o <- read.table('long_other_3ss_score.txt')

z <- rbind(short_c,short_o,long_c,long_o)
type <- c(rep('short_cancer',443),rep('short_other',16303),rep('long_cancer',1355),rep('long_other',171829))
z <- as.data.frame(z)
colnames(z)<- c('sequence','3SS_score')
z$type <- type
colnames(z)

#ggplot(z,aes(y=z$`3SS_score`,x=z$type)) + geom_boxplot(width=0.8,fill='moccasin')+theme_bw()+xlab("")+ylab("3' splice site score")+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+ scale_x_discrete(limits=c("short_cancer", "short_other", "long_cancer","long_other"))+theme(panel.grid=element_blank())+ylim(-20,20)

short_c <- read.table('short_cancer_5ss_score.txt')
short_o <- read.table('short_other_5ss_score.txt')
long_c <- read.table('long_cancer_5ss_score.txt')
long_o <- read.table('long_other_5ss_score.txt')

short_o$ss <- str_sub(short_o$V1,start=4L,end=5L)
long_o$ss <- str_sub(long_o$V1,start=4L,end=5L)
long_c$ss <- str_sub(long_c$V1,start=4L,end=5L)
short_c$ss <- str_sub(short_c$V1,start=4L,end=5L)

zz <- rbind(short_c,short_o,long_c,long_o)
type <- c(rep('short_cancer',443),rep('short_other',16303),rep('long_cancer',1355),rep('long_other',171829))
zz <- as.data.frame(zz)
colnames(zz)<- c('sequence','5SS_score')
zz$type <- type
colnames(zz)


ss3 <- ggplot(z,aes(y=z$`3SS_score`,x=z$type)) + geom_boxplot(width=0.8,fill='moccasin')+theme_bw()+xlab("")+ylab("3' splice site score")+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+ scale_x_discrete(limits=c("short_cancer", "short_other", "long_cancer","long_other"))+theme(panel.grid=element_blank())+ylim(-20,20)
ss5 <- ggplot(zz,aes(y=zz$`5SS_score`,x=zz$type)) + geom_boxplot(width=0.8,fill='moccasin')+theme_bw()+xlab("")+ylab("5' splice site score")+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+ scale_x_discrete(limits=c("short_cancer", "short_other", "long_cancer","long_other"))+theme(panel.grid=element_blank())+ylim(-20,20)
ggpubr::ggarrange(ss5,ss3, nrow = 1, ncol = 2)

zz$SS3_score <- z$`3SS_score`
zz$SCORE_SUM <- zz$`5SS_score`+zz$SS3_score
ggplot(zz,aes(y=zz$SCORE_SUM,x=zz$type)) + geom_boxplot(width=0.8,fill='moccasin')+theme_bw()+xlab("")+ylab("3ss score + 5ss score")+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+ scale_x_discrete(limits=c("short_cancer", "short_other", "long_cancer","long_other"))+theme(panel.grid=element_blank())







###########conscore

d <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\015conservation\\up_everyscore.txt',head=TRUE,row.names = 6)
d$site <- rownames(d)
d<-d[,-1]
library(reshape2)
mydata <- melt(d,id.vars="site",variable.name="type",value.name="score")
mydata$type <- factor(mydata$type,levels=c('short_cancer','short_other','long_cancer','long_other'))
up <- ggplot(data = mydata, mapping = aes(x = as.numeric(site), y = as.numeric(mydata$score), group=type,colour = type)) + geom_line(size=1.1)+theme_bw()+theme(panel.grid=element_blank(),axis.text.x = element_text( angle = 90))+ expand_limits(y=c(0,1))+ scale_x_continuous(breaks=c(-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0))+scale_color_manual(values=c("#EC7357","gold3","forestgreen","#9DC3C1"))+ xlab("") + ylab("") 


d <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\015conservation\\down_everyscore.txt',head=TRUE,row.names = 6)
d$site <- rownames(d)
d<-d[,-1]
mydata2 <- melt(d,id.vars="site",variable.name="type",value.name="score")
mydata2$type <- factor(mydata2$type,levels=c('short_cancer','short_other','long_cancer','long_other'))
down <- ggplot(data = mydata2, mapping = aes(x = as.numeric(site), y = as.numeric(mydata2$score), group=type,colour = type)) + geom_line(size=1.1)+theme_bw()+theme(panel.grid=element_blank(),axis.text.x = element_text( angle = 90))+ expand_limits(y=c(0,1))+ scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100))+scale_color_manual(values=c("#EC7357","gold3","forestgreen","#9DC3C1"))+ xlab("") + ylab("") 


ggpubr::ggarrange(up,down, nrow = 1, ncol = 2)

##########intron length
 data <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\exon_class\\allexon_length_class.txt',head=TRUE)
long_c <- filter(data,data$type=='long_cancer') 
long_o <- filter(data,data$type=='long_other') 
short_c <- filter(data,data$type=='short_cancer') 
short_o <- filter(data,data$type=='short_other') 


###########exon ATCG 1027
data <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\exon_class\\exon_atcg.txt',head=TRUE)
mydata1 <- melt(data,id.vars=c('name','class'),variable.name="necleutide",measure.vars=c('A','T','C','G'))
mydata1$class <- factor(mydata1$class,levels=c('short_cancer','short_other','long_cancer','long_other'))
exon <- ggplot(mydata1,aes(y=as.numeric(as.character(mydata1$value)),x=mydata1$necleutide,fill=mydata1$class)) + geom_boxplot()+ theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.title = element_blank(), legend.key = element_blank())+ labs(x = '', y = 'frequency')+ggtitle('exon')+     scale_fill_manual(values = c('darkorange', 'tan','darkgreen','lightgreen'))

data <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\atcg\\atcg_down100.txt',head=TRUE)
mydata2 <- melt(data,id.vars=c('name','class'),variable.name="necleutide",measure.vars=c('A','T','C','G'))
mydata2$class <- factor(mydata2$class,levels=c('short_cancer','short_other','long_cancer','long_other'))
down100 <- ggplot(mydata2,aes(y=as.numeric(as.character(mydata2$value)),x=necleutide,fill=class)) + geom_boxplot()+ theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.title = element_blank(), legend.key = element_blank())+ labs(x = '', y = 'log2(normalized_counts)')+ggtitle('down100')+     scale_fill_manual(values = c('darkorange', 'gray','darkgreen','gray'))

data <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\atcg\\atcg_down300.txt',head=TRUE)
mydata3 <- melt(data,id.vars=c('name','class'),variable.name="necleutide",measure.vars=c('A','T','C','G'))
mydata3$class <- factor(mydata3$class,levels=c('short_cancer','short_other','long_cancer','long_other'))
down300 <- ggplot(mydata3,aes(y=as.numeric(as.character(mydata3$value)),x=necleutide,fill=class)) + geom_boxplot()+ theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.title = element_blank(), legend.key = element_blank())+ labs(x = '', y = 'frequency')+ggtitle('down300')+     scale_fill_manual(values =c('darkorange', 'tan','darkgreen','lightgreen'))

data <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\atcg\\atcg_down500.txt',head=TRUE)
mydata4 <- melt(data,id.vars=c('name','class'),variable.name="necleutide",measure.vars=c('A','T','C','G'))
mydata4$class <- factor(mydata4$class,levels=c('short_cancer','short_other','long_cancer','long_other'))
down500 <- ggplot(mydata4,aes(y=as.numeric(as.character(mydata4$value)),x=necleutide,fill=class)) + geom_boxplot()+ theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.title = element_blank(), legend.key = element_blank())+ labs(x = '', y = 'log2(normalized_counts)')+ggtitle('down500')+     scale_fill_manual(values = c('darkorange', 'gray','darkgreen','gray'))


data <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\atcg\\atcg_up100.txt',head=TRUE)
mydata5 <- melt(data,id.vars=c('name','class'),variable.name="necleutide",measure.vars=c('A','T','C','G'))
mydata5$class <- factor(mydata5$class,levels=c('short_cancer','short_other','long_cancer','long_other'))
up100 <- ggplot(mydata5,aes(y=as.numeric(as.character(mydata5$value)),x=necleutide,fill=class)) + geom_boxplot()+ theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.title = element_blank(), legend.key = element_blank())+ labs(x = '', y = 'log2(normalized_counts)')+ggtitle('up100')+     scale_fill_manual(values = c('darkorange', 'gray','darkgreen','gray'))

data <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\atcg\\atcg_up300.txt',head=TRUE)
mydata6 <- melt(data,id.vars=c('name','class'),variable.name="necleutide",measure.vars=c('A','T','C','G'))
mydata6$class <- factor(mydata6$class,levels=c('short_cancer','short_other','long_cancer','long_other'))
up300 <- ggplot(mydata6,aes(y=as.numeric(as.character(mydata6$value)),x=necleutide,fill=class)) + geom_boxplot()+ theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.title = element_blank(), legend.key = element_blank())+ labs(x = '', y = 'frequency')+ggtitle('up300')+     scale_fill_manual(values =c('darkorange', 'tan','darkgreen','lightgreen'))

data <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\atcg\\atcg_up500.txt',head=TRUE)
mydata7 <- melt(data,id.vars=c('name','class'),variable.name="necleutide",measure.vars=c('A','T','C','G'))
mydata7$class <- factor(mydata7$class,levels=c('short_cancer','short_other','long_cancer','long_other'))
up500 <- ggplot(mydata7,aes(y=as.numeric(as.character(mydata7$value)),x=necleutide,fill=class)) + geom_boxplot()+ theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.title = element_blank(), legend.key = element_blank())+ labs(x = '', y = 'log2(normalized_counts)')+ggtitle('up500')+     scale_fill_manual(values = c('darkorange', 'gray','darkgreen','gray'))

######aa
a <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\017protein\\swiss_aa.txt',head=TRUE)
a$aa <- factor(a$aa,levels=c('L','I','F','W','V','M','C','Y','A','T','E','G','S','Q','D','R','K','N','H','P'))
a$aa <- factor(a$aa,levels=c('I','V','L','F','C','M','A','G','T','S','W','Y','P','H','E','Q','D','N','K','R'))

ggplot(data = a, mapping = aes(x = aa, y = as.numeric(as.character(percent)), fill = class)) + geom_bar(stat = 'identity', position = 'dodge',width=0.7) + scale_fill_brewer(palette = 'Accent')+theme_bw()+labs(x="",y="")

######iupred

a <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\017protein\\iupred\\a.txt',head=TRUE)
cc <- c('slateblue3','steelblue3','yellow3')
a$type = factor(a$type, levels=c('high','mid','low')) ## 设置柱条的顺序
a$class = factor(a$class, levels=c('short_cancer','short_other','long_cancer','long_other')) ## 设置柱条的顺序
p <- ggplot(a, aes(x =class, y = iupred_score, fill =as.factor(type)) )+    geom_bar(stat = "identity", position = "fill",width=0.6,col='black')
p + theme_bw()+scale_fill_manual(values=cc) +theme(axis.text.x  = element_text(angle=45, hjust=1,vjust=1 ),legend.title=element_blank())+ labs(x = "")+theme(panel.grid=element_blank())


p <- ggplot(a, aes(x =class, y = anchor_score, fill =as.factor(type)) )+    geom_bar(stat = "identity", position = "fill",width=0.6)
p + theme_bw()+scale_fill_manual(values=cc) +theme(axis.text.x  = element_text(angle=45, hjust=1,vjust=1 ),legend.title=element_blank())+ labs(x = "")+theme(panel.grid=element_blank())



a<-c(34.69,37.62,29.87,55.78)
b <- c('short_cancer','short_other','long_cancer','long_other')
c <- as.data.frame(cbind(a,b))
c$b = factor(c$b, levels=c('short_cancer','short_other','long_cancer','long_other'))
ggplot(c,aes(x=b,y=as.numeric(as.character(a))))+geom_bar(stat = "identity",width=0.7,fill="grey66",col='black')+theme_bw()+theme(axis.text.x  = element_text(angle=45, hjust=1,vjust=1 ))+labs(x="",y="")+ylim(0,60)+theme(panel.grid=element_blank())



#########
a <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\number.txt',head=TRUE)
a <- filter(a,a$type=='a')
a$class<- factor(a$class,levels=c('cancer','other_SE','other_all'))
p1 <- ggplot(a,aes(y=number,x=class,fill=length))+ geom_bar(stat = "identity", position = "fill",width=0.5)+theme_bw()+scale_fill_manual(name='',values = c("darkgreen", "#E69F00"))  +theme(panel.grid=element_blank()) +ylab('')+xlab("")+theme_classic()
fisher.test(matrix(c(494,1571,17467,179284),nrow=2))$p.value
#[1] 1.67052e-90

b <- filter(a,a$type=='b')
b$class <- factor(b$class,levels=c('all_exon','SE','cancer_exon'))
p2 <- ggplot(b,aes(y=number,x=class,fill=length))+ geom_bar(stat = "identity", position = "fill",width=0.5)+theme_bw()+scale_fill_manual(name='',values = c("darkgreen", "#E69F00"))  +theme(panel.grid=element_blank()) +ylab('')+xlab("")+theme_classic()
ggpubr::ggarrange(p1,p2, nrow = 1, ncol = 2)

###
gene <- read.table('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\some_file\\517_gene_name.txt')
gene_ <- bitr(gene$V1,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db",drop = T)$ENTREZID
ego <- enrichKEGG(gene = gene_ ,keyType ="kegg",organism ='hsa', pvalueCutoff = 0.05, pAdjustMethod = "BH",qvalueCutoff = 0.05)
write.table(summary(ego),file='D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\some_file\\517_gene_kegg.txt',sep='\t')
d5 <- read.table('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\some_file\\517_gene_kegg.txt',head=TRUE,sep='\t')
d5$logP <- -log10(d5$p.adjust)
d5 <- d5[order(d5$logP),]
d5$Description<- factor(d5$Description,levels=d5$Description)
d5 <- d5[,c(2,10)]
p5 <- ggplot() +  geom_bar(data = d5,  aes(x = Description, y = logP,fill=logP),colour='black',stat = "identity",width = 0.6,  position = position_dodge(width = 0.9))+ coord_flip()+theme_bw()+theme(panel.grid=element_blank()) +ylab("-logP")+xlab("")+ theme(axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 12))+ scale_fill_gradient(low = "pink", high = "darkred")+ guides(fill=FALSE)+ggtitle('sgR10_AS')


###
par(mfrow = c(1,3))
info = c(494,1571)
piepercent = paste(round(100*info/sum(info)), "%")
names = c("short", "long")
cols = c( "#FFC90E","#22B14C")
pie(info, labels=piepercent, col=cols)
legend("topright", names, fill=cols)


info3 = c(17031,139729)
piepercent3 = paste(round(100*info3/sum(info3)), "%")
pie(info3, labels=piepercent3, col=cols)
legend("topright", names, fill=cols)


info2 = c(17467,179284)
piepercent2 = paste(round(100*info2/sum(info2)), "%")
pie(info2, labels=piepercent2, col=cols)
legend("topright", names, fill=cols)



a <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\number_2.txt',head=TRUE)
a$class<- factor(a$class,levels=c('cancer','SE','all'))
a$length<- factor(a$length,levels=c('short','long'))
ggplot(a,aes(y=number,x=class,fill=length))+ geom_bar(stat = "identity", position = "fill",width=0.5)+theme_bw()+scale_fill_manual(name='',values = c( "#E69F00","darkgreen"))  +theme(panel.grid=element_blank()) +ylab('')+xlab("")+theme_classic()
a$length<- factor(a$length,levels=c('long','short'))
