setwd('D:\\wanglab\\smallexon\\new_60\\03rbp')
d <- read.table('all_shrbp_psi_60.txt',head=TRUE)
d <- read.table('all_shrbp_short_psi.txt',head=TRUE)


f_na<-function(x) sum(is.na(x))
na <- apply(d,1,f_na)
d$NA_number <- na
library(dplyr)
d <- filter(d,d$NA_number<=150)

name <- read.table('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\some_file\\517_shortexon.txt',head=TRUE)
n <-c()
library(stringr)
for (i in seq(from=1,to=517)){
  n <- c(n,str_c(str_split_fixed(name$name[i],':',9)[1,1],str_c(str_split_fixed(name$name[i],':',9)[1,5], str_split_fixed(name$name[i],':',9)[1,6],sep='-'),sep=':'))
}
name <- as.data.frame(name[,3])
name$name_2 <- n
name <- unique(name)

n <-c()
for (i in seq(from=1,to=7618)){
  n <- c(n,str_c(str_split_fixed(d$Node[i],':',2)[1,1],str_c(str_split_fixed(d$Node[i],'-',4)[1,2], str_split_fixed(d$Node[i],'-',4)[1,3],sep='-'),sep=':'))
}
d$name <-n

d_ <- d[which(d$name %in% name$name_2),]

f_na<-function(x) sum(is.na(x))
na <- apply(d_,1,f_na)
d_$NA_number <- na
library(dplyr)
d_[is.na(d_)]=0
d_ <- filter(d_,d_$NA_number<=220)
pheatmap(t(d_[,4:230]),color = c(colorRampPalette(colors = c("royalblue","white"))(length(bk)/2),colorRampPalette(colors = c("white","maroon1"))(length(bk)/2)),na_col="snow2",show_colnames = F,fontsize_row=6)

write.table(d_,file='rbp_CASE.txt',sep='\t',row.names = F)
dt <- read.table('rbp_CASE.txt',head=TRUE,row.names = 1)






f_<-function(x) (sum(x> (0.1),na.rm=TRUE)+sum(x< (-0.1),na.rm=TRUE))
num <- apply(d_[,4:230],1,f_)
d_$num <- num

TEMP <-c()
a <- apply(d_[,4:230],2,f_)
for (i in seq(from=1,to=227)){
  if (a[i]>=111){
    TEMP<- c(TEMP,i+3)
  }
}
###517/3=172
d_2 <- d_[,TEMP]

pheatmap(t(d_2[,4:45]),color = c(colorRampPalette(colors = c("royalblue","white"))(length(bk)/2),colorRampPalette(colors = c("white","maroon1"))(length(bk)/2)),na_col="snow2",show_colnames = F,fontsize_row=8)

d_2$num <- d_$num
d_2$name <-d_$name
d_2$gene_name <-d_$Gene
d_2 <- filter(d_2,d_2$num>=10)

colnames(d_2) <- str_split_fixed(colnames(d_2),'_',2)[,1]


###
D <- read.table('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\some_file\\shortexon_merge_t.txt',sep='\t',head=TRUE)

myMat.cor <-cor(as.matrix(D[,2:518]), method=c("pearson"),use="pairwise.complete.obs")

bk <- c(seq(-1,0,by=0.01),seq(0,1,by=0.01))

pheatmap(myMat.cor,show_rownames = F,show_colnames = F, color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),na_col="snow2")

out <- pheatmap(myMat.cor,show_rownames =T,show_colnames =T, color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),na_col="snow2")
out2 <- myMat.cor[out$tree_row[["order"]],]
write.table(myMat.cor[out$tree_row[["order"]],],file='D:\\wanglab\\smallexon\\new_60\\03rbp\\shortexon_cluster_cor.txt',sep='\t')

pdf("D:\\wanglab\\smallexon\\new_60\\03rbp\\exon_cluster_red.pdf", width=6, height=10)
print(p)
graphics.off()


myMat.cor <-abs(myMat.cor)
myMat.cor <- as.data.frame(myMat.cor)
p2 <- pheatmap(myMat.cor,show_rownames = F,show_colnames = F, color = c(colorRampPalette(colors = c("white","firebrick3"))(500)),na_col="snow2")

pdf("D:\\wanglab\\smallexon\\new_60\\03rbp\\exon_cluster_red_abs.pdf", width=6, height=5.5)
print(p2)
graphics.off()

out2 <- myMat.cor[p2$tree_row$order,]
out2 <- out2[p2$tree_row$order]
write.table(out2,file='D:\\wanglab\\smallexon\\new_60\\03rbp\\shortexon_cluster_abs_cor.txt',sep='\t')

name <- read.table('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\some_file\\517_shortexon.txt',head=TRUE)
n <-c()
library(stringr)
for (i in seq(from=1,to=517)){
  n <- c(n,str_c(str_split_fixed(name$name[i],':',9)[1,1],str_c(str_split_fixed(name$name[i],':',9)[1,5], str_split_fixed(name$name[i],':',9)[1,6],sep='-'),sep=':'))
}
name <- as.data.frame(name[,3])
name$name_2 <- n
name <- unique(name)

exon_cluster <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\shortexon_cluster_abs_cor.txt',head=TRUE)
exon_order <- as.character(exon_cluster$a)
exon_order <- as.data.frame(exon_order)
colnames(exon_order)<- c('name_2')
t <- merge(exon_order,name,by='name_2')
t <- t[p2$tree_row$order]
rownames(t) <- t$name_2
write.table(name,file='D:\\wanglab\\smallexon\\new_60\\03rbp\\gene_name.txt',sep='\t')


setwd('D:\\wanglab\\smallexon\\new_60\\03rbp\\rbp_expression_tcga')
BLCA <- read.table('BLCA_rbp_expre.txt',head=TRUE,sep='\t')
BRCA <- read.table('BRCA_rbp_expre.txt',head=TRUE,sep='\t')
COAD <- read.tabl('COAD_rbp_expre.txt',head=TRUE,sep='\t')
ESCA <- read.table('ESCA_rbp_expre.txt',head=TRUE,sep='\t')
GBM <- read.table('GBM_rbp_expre.txt',head=TRUE,sep='\t')
HNSC <- read.table('HNSC_rbp_expre.txt',head=TRUE,sep='\t')
KICH <- read.table('KICH_rbp_expre.txt',head=TRUE,sep='\t')
KIRC <- read.table('KIRC_rbp_expre.txt',head=TRUE,sep='\t')
KIRP <- read.table('KIRP_rbp_expre.txt',head=TRUE,sep='\t')
LIHC <- read.table('LIHC_rbp_expre.txt',head=TRUE,sep='\t')
LUAD <- read.table('LUAD_rbp_expre.txt',head=TRUE,sep='\t')
LUSC <- read.table('LUSC_rbp_expre.txt',head=TRUE,sep='\t')
PAAD <- read.table('PAAD_rbp_expre.txt',head=TRUE,sep='\t')
PRAD <- read.table('PRAD_rbp_expre.txt',head=TRUE,sep='\t')
READ <- read.table('READ_rbp_expre.txt',head=TRUE,sep='\t')
STAD <- read.table('STAD_rbp_expre.txt',head=TRUE,sep='\t')
THCA <- read.table('THCA_rbp_expre.txt',head=TRUE,sep='\t')
UCEC <- read.table('UCEC_rbp_expre.txt',head=TRUE,sep='\t')
t <- merge(BLCA,BRCA,by='Hybridization.REF')
t <- merge(t,COAD,by='Hybridization.REF')
t <- merge(t,ESCA,by='Hybridization.REF')
t <- merge(t,GBM,by='Hybridization.REF')
t <- merge(t,HNSC,by='Hybridization.REF')
t <- merge(t,KICH,by='Hybridization.REF')
t <- merge(t,KIRC,by='Hybridization.REF')
t <- merge(t,KIRP,by='Hybridization.REF')
t <- merge(t,LIHC,by='Hybridization.REF')
t <- merge(t,LUAD,by='Hybridization.REF')
t <- merge(t,LUSC,by='Hybridization.REF')
t <- merge(t,PAAD,by='Hybridization.REF')
t <- merge(t,PRAD,by='Hybridization.REF')
t <- merge(t,READ,by='Hybridization.REF')
t <- merge(t,STAD,by='Hybridization.REF')
t <- merge(t,THCA,by='Hybridization.REF')
t <- merge(t,UCEC,by='Hybridization.REF')

write.table(t,file='D:\\wanglab\\smallexon\\new_60\\03rbp\\rbp_geneexpression_tcga.txt',sep='\t')


D <- read.table('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\some_file\\shortexon_merge_t.txt',sep='\t',head=TRUE,row.names = 1)
rbp <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\rbp_geneexpression_tcga_t.txt',head=TRUE,row.names=1)


myMat.cor <-cor(rbp, method=c("pearson"),use="pairwise.complete.obs")
bk <- c(seq(-1,0,by=0.01),seq(0.01,1,by=0.01))
pheatmap(myMat.cor,show_rownames = F,show_colnames = F, color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),na_col="snow2",breaks=bk)


D <- read.table('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\some_file\\shortexon_merge_t.txt',sep='\t',head=TRUE,row.names = 1)
rbp <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\rbp_geneexpression_tcga_t.txt',head=TRUE,row.names=1)
D$sample_name <- row.names(D)
rbp$sample_name <- row.names(rbp)
exon_rbp <- merge(D,rbp,by='sample_name')
a <- c()
b <- c()
exon_name <-c()
rbp_name <- c()
for (i in seq(from=2,to=518)){
  for (j in seq(from=519,to=1866)){
    temp <-c(colnames(exon_rbp)[i],colnames(exon_rbp)[j],cor(exon_rbp[,i],exon_rbp[,j],use="complete.obs"))
   # print(temp)
    exon_name <-c(exon_name,colnames(exon_rbp)[i])
    a <- c(a,cor(exon_rbp[,i],exon_rbp[,j],use="complete.obs"))
  }
  rbp_name<-c(rbp_name,colnames(exon_rbp)[j])
  b <- rbind(b,a)
  a <-c()
}
b<-as.data.frame(b)
colnames(b)<-colnames(exon_rbp)[519:1866]
rownames(b)<-colnames(exon_rbp)[2:518]


exon_cluster <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\shortexon_cluster_abs_cor.txt',head=TRUE)
exon_order <- as.character(exon_cluster$a)
exon_order <- as.data.frame(exon_order)

cluster1_name <- colnames(exon_cluster)[423:517]
cluster1 <- as.data.frame(t(b[rownames(b) %in% cluster1_name,]))
f<-function(x) sum(abs(x)>0.5)
num <- apply(cluster1,1,f)
cluster1$num05 <- num
cluster1$name <- rownames(cluster1)
cluster1 <- filter(cluster1,cluster1$num05>=48)
rownames(cluster1)<- cluster1$name
bk <- c(seq(-1,0,by=0.01),seq(0.01,1,by=0.01))
pheatmap(cluster1[1:31,1:95],show_colnames = F,border_color = NA,color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),breaks=bk)



cluster2_name <- colnames(exon_cluster)[373:410]
cluster2 <- as.data.frame(t(b[rownames(b) %in% cluster2_name,]))
f<-function(x) sum(abs(x)>0.5)
num <- apply(cluster2,1,f)
cluster2$num05 <- num
cluster2$name <- rownames(cluster2)
cluster2 <- filter(cluster2,cluster2$num05>=13)
rownames(cluster2)<- cluster2$name
bk <- c(seq(-1,0,by=0.01),seq(0.01,1,by=0.01))
pheatmap(cluster2[,1:38],show_colnames = F,border_color = NA,color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),breaks=bk)


cluster3_name <- colnames(exon_cluster)[344:372]
cluster3 <- as.data.frame(t(b[rownames(b) %in% cluster3_name,]))
f<-function(x) sum(abs(x)>0.5)
num <- apply(cluster3,1,f)
cluster3$num05 <- num
cluster3$name <- rownames(cluster3)
cluster3 <- filter(cluster3,cluster3$num05>=15)
rownames(cluster3)<- cluster3$name
bk <- c(seq(-1,0,by=0.01),seq(0.01,1,by=0.01))
pheatmap(cluster3[,1:29],show_colnames = F,border_color = NA,color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),breaks=bk,cellwidth=1.3, cellheight=8,fontsize_row=8)
pheatmap(cluster3[,1:29],show_colnames = F,border_color = NA,color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),breaks=bk,cellwidth=4.5, cellheight=6.6,fontsize_row=7)


cluster4_name <- colnames(exon_cluster)[193:218]
cluster4 <- as.data.frame(t(b[rownames(b) %in% cluster4_name,]))
f<-function(x) sum(abs(x)>0.5)
num <- apply(cluster4,1,f)
cluster4$num05 <- num
cluster4$name <- rownames(cluster4)
cluster4 <- filter(cluster4,cluster4$num05>=13)
rownames(cluster4)<- cluster4$name
bk <- c(seq(-1,0,by=0.01),seq(0.01,1,by=0.01))
pheatmap(cluster4[,1:29],show_colnames = F,border_color = NA,color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),breaks=bk)



cluster5_name <- colnames(exon_cluster)[181:192]
cluster5 <- as.data.frame(t(b[rownames(b) %in% cluster5_name,]))
f<-function(x) sum(abs(x)>0.5)
num <- apply(cluster5,1,f)
cluster5$num05 <- num
cluster5$name <- rownames(cluster5)
cluster5 <- filter(cluster5,cluster5$num05>=13)
rownames(cluster5)<- cluster5$name
bk <- c(seq(-1,0,by=0.01),seq(0.01,1,by=0.01))
pheatmap(cluster4[,1:29],show_colnames = F,border_color = NA,color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),breaks=bk)


a <-c(c(193:218),c(249:274),c(423:517))
cluster14_name <- colnames(exon_cluster)[a]
cluster14 <- as.data.frame(t(b[rownames(b) %in% cluster14_name,]))
f<-function(x) sum(abs(x)>0.5)
num <- apply(cluster14,1,f)
cluster14$num05 <- num
cluster14$name <- rownames(cluster14)
cluster14 <- filter(cluster14,cluster14$num05>=70)
rownames(cluster14)<- cluster14$name
bk <- c(seq(-1,0,by=0.01),seq(0.01,1,by=0.01))
pheatmap(cluster14[,1:147],show_colnames = F,border_color = NA,color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),breaks=bk)
pheatmap(cluster14[,1:147],show_colnames = F,border_color = NA,color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),breaks=bk,cellwidth=1.3, cellheight=8,fontsize_row=8)



go <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\cluster1_go.txt',head=TRUE,sep='\t')
barplot(go$P,names=go$name,horiz=TRUE,border=FALSE)
go$name = factor(go$name, levels=go$name)
ggplot(data =go, mapping = aes(y = go$P, x = go$name,width = 0.5)) + geom_bar(stat = 'identity')+coord_flip()+theme_classic()+ylab("")+xlab("")+theme(panel.grid=element_blank())




exon_cluster <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\shortexon_cluster_abs_cor.txt',head=TRUE)
exon_order <- as.character(exon_cluster$a)
exon_order <- as.data.frame(exon_order)

cluster_14 <- c(c(193:218),c(249:274),c(423:517))
cluster_14_ <- exon_cluster[cluster_14,]
cluster_14_ <- cluster_14_[,cluster_14+1]
rownames(cluster_14_) <-colnames(cluster_14_)

cluster_14_ <- cluster_14_[,-c(8,29,34,40,41,45,118)]
cluster_14_ <- cluster_14_[-c(8,29,34,40,41,45,118),]
cluster_14_ <- cluster_14_[,-37]
cluster_14_ <- cluster_14_[-37,]
pheatmap(cluster_14_,cluster_cols = F,cluster_rows=F,show_rownames = F,show_colnames = F,color = c(colorRampPalette(colors = c("white","firebrick3"))(500)), border_color=NA)
write.table(cluster_14_,file='D:\\wanglab\\smallexon\\new_60\\03rbp\\cluster\\cluster_14_new.txt',sep='\t')


c2 <- c(373:410)
cluster_2_ <- exon_cluster[c2,]
cluster_2_ <- cluster_2_[,c2+1]
rownames(cluster_2_) <-colnames(cluster_2_)
pheatmap(cluster_2_,cluster_cols = F,cluster_rows=F,show_rownames = F,show_colnames = F,color = c(colorRampPalette(colors = c("white","firebrick3"))(500)), border_color=NA)
write.table(cluster_2_,file='D:\\wanglab\\smallexon\\new_60\\03rbp\\cluster\\cluster_2.txt',sep='\t')


cluster_2_ <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\cluster\\cluster_2.txt',head=TRUE)
cluster_2_net <- as.data.frame(rbp_clip[rbp_clip$V13 %in% cluster_2_$name,])
da2 <- cbind(as.character(cluster_2_net$V1),as.character(cluster_2_net$V12))
g <- graph_from_data_frame(da2, directed = FALSE, vertices = NULL)
V(g)
V(g)$color=c(rep("darkolivegreen3",37),rep("tan2",19))
plot(g,vertex.size=map(degree(g),c(1,20)),vertex.label.cex=map(degree(g),c(0.6,0.8)),
     layout=layout_on_sphere,vertex.frame.color=V(g)$color,vertex.label.color='grey2')
write.table(cluster_2_net,file='D:\\wanglab\\smallexon\\new_60\\03rbp\\cluster\\cluster_2_net.txt',sep='\t')






c3 <- c(344:372)
cluster_3_ <- exon_cluster[c3,]
cluster_3_ <- cluster_3_[,c3+1]
rownames(cluster_3_) <-colnames(cluster_3_)
cluster_3_ <- cluster_3_[,-25]
cluster_3_ <- cluster_3_[-25,]
pheatmap(cluster_3_,cluster_cols = F,cluster_rows=F,show_rownames = F,show_colnames = F,color = c(colorRampPalette(colors = c("white","firebrick3"))(500)), border_color=NA)
write.table(cluster_3_,file='D:\\wanglab\\smallexon\\new_60\\03rbp\\cluster\\cluster_3.txt',sep='\t')


cluster_3_ <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\cluster\\cluster_3.txt',head=TRUE)
cluster_3_net <- as.data.frame(rbp_clip[rbp_clip$V13 %in% cluster_3_$name,])
da2 <- cbind(as.character(cluster_3_net$V1),as.character(cluster_3_net$V12))
g <- graph_from_data_frame(da2, directed = FALSE, vertices = NULL)
V(g)
V(g)$color=c(rep("darkolivegreen3",37),rep("tan2",17))
plot(g,vertex.size=map(degree(g),c(1,20)),vertex.label.cex=map(degree(g),c(0.6,0.8)),
     layout=layout_on_sphere,vertex.frame.color=V(g)$color,vertex.label.color='grey2')
write.table(cluster_3_net,file='D:\\wanglab\\smallexon\\new_60\\03rbp\\cluster\\cluster_3_net.txt',sep='\t')



#########
AQR <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\exon_length\\AQR_ENCSR624OUI_SE.MATS.JCEC.txt',head=TRUE)
AQR$length <- AQR$exonEnd-AQR$exonStart_0base
AQR_ <- filter(AQR,abs(AQR$IncLevelDifference)>=0.1 & AQR$PValue<0.05)
boxplot(AQR$length,AQR_$length,ylim=c(0,300),names=c('all','AQR'),col='lightskyblue')

AQR <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\exon_length\\AQR.txt',head=TRUE)
ggplot(AQR,aes(y=as.numeric(AQR$number),x=AQR$length,fill=AQR$type))+geom_bar(stat = "identity", position = "fill",width=0.7)+ geom_text(mapping = aes(label = number), size = 3.5,vjust = 3.5, hjust = .5, colour = 'black', position = position_fill())+theme_bw()+scale_fill_manual(name='',values = c("violet", "skyblue2"))+ylab("")+xlab("exon_length")+theme(panel.grid=element_blank()) 



RBFOX2 <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\exon_length\\RBFOX2_ENCSR336DFS_SE.MATS.JCEC.txt',head=TRUE)
RBFOX2$length <- RBFOX2$exonEnd-RBFOX2$exonStart_0base
RBFOX2_ <- filter(RBFOX2,abs(RBFOX2$IncLevelDifference)>=0.1 & RBFOX2$PValue<0.05)
boxplot(RBFOX2$length,RBFOX2_$length,ylim=c(0,300),names=c('all','RBFOX2'),col='lightskyblue')

RBFOX2 <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\exon_length\\RBFOX2_ENCSR336DFS_SE.txt',head=TRUE)
ggplot(RBFOX2,aes(y=as.numeric(RBFOX2$number),x=RBFOX2$length,fill=RBFOX2$type))+geom_bar(stat = "identity", position = "fill",width=0.7)+ geom_text(mapping = aes(label = number), size = 3.5,vjust = 3.5, hjust = .5, colour = 'black', position = position_fill())+theme_bw()+scale_fill_manual(name='',values = c("violet", "skyblue2"))+ylab("")+xlab("exon_length")+theme(panel.grid=element_blank()) 


celf1 <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\exon_length\\CELF1.txt',head=TRUE)
ggplot(celf1,aes(y=as.numeric(celf1$number),x=celf1$length,fill=type))+geom_bar(stat = "identity", position = "fill",width=0.7)+ geom_text(mapping = aes(label = number), size = 3.5,vjust = 3.5, hjust = .5, colour = 'black', position = position_fill())+theme_bw()+scale_fill_manual(name='',values = c("violet", "skyblue2"))+ylab("")+xlab("exon_length")+theme(panel.grid=element_blank()) 

CELF1 <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\exon_length\\CELF1_ENCSR605MFS_SE.MATS.JCEC.txt',head=TRUE)
CELF1$length <- CELF1$exonEnd-CELF1$exonStart_0base
CELF1_ <- filter(CELF1,abs(CELF1$IncLevelDifference)>=0.1 & CELF1$PValue<0.05)
boxplot(CELF1$length,CELF1_$length,ylim=c(0,300),names=c('all','CELF1'),col='lightskyblue')


celf1 <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\exon_length\\U2AF2.txt',head=TRUE)
ggplot(celf1,aes(y=as.numeric(celf1$number),x=celf1$length,fill=type))+geom_bar(stat = "identity", position = "fill",width=0.7)+ geom_text(mapping = aes(label = number), size = 3.5,vjust = 3.5, hjust = .5, colour = 'black', position = position_fill())+theme_bw()+scale_fill_manual(name='',values = c("violet", "skyblue2"))+ylab("")+xlab("exon_length")+theme(panel.grid=element_blank()) 

CELF1 <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\exon_length\\U2AF2_ENCSR904CJQ_SE.MATS.JCEC.txt',head=TRUE)
CELF1$length <- CELF1$exonEnd-CELF1$exonStart_0base
CELF1_ <- filter(CELF1,abs(CELF1$IncLevelDifference)>=0.1 & CELF1$PValue<0.05)
boxplot(CELF1$length,CELF1_$length,ylim=c(0,300),names=c('all','CELF1'),col='lightskyblue')



PTBP1 <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\exon_length\\PTBP1_ENCSR527IVX_SE.MATS.JCEC.txt',head=TRUE)
PTBP1$length <- PTBP1$exonEnd-PTBP1$exonStart_0base
PTBP1_ <- filter(PTBP1,abs(PTBP1$IncLevelDifference)>=0.1 & PTBP1$PValue<0.05)
boxplot(PTBP1$length,PTBP1_$length,ylim=c(0,300),names=c('all','PTBP1'),col='lightskyblue')

PTBP1 <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\exon_length\\PTBP1_ENCSR527IVX_SE.txt',head=TRUE)
ggplot(PTBP1,aes(y=as.numeric(PTBP1$number),x=PTBP1$length,fill=type))+geom_bar(stat = "identity", position = "fill",width=0.7)+ geom_text(mapping = aes(label = number), size = 3.5,vjust = 3.5, hjust = .5, colour = 'black', position = position_fill())+theme_bw()+scale_fill_manual(name='',values = c("violet", "skyblue2"))+ylab("")+xlab("exon_length")+theme(panel.grid=element_blank()) 



SRSF1 <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\exon_length\\SRSF1_ENCSR066VOO_SE.txt',head=TRUE)
ggplot(SRSF1,aes(y=as.numeric(SRSF1$number),x=SRSF1$length,fill=type))+geom_bar(stat = "identity", position = "fill",width=0.7)+ geom_text(mapping = aes(label = number), size = 3.5,vjust = 3.5, hjust = .5, colour = 'black', position = position_fill())+theme_bw()+scale_fill_manual(name='',values = c("violet", "skyblue2"))+ylab("")+xlab("exon_length")+theme(panel.grid=element_blank()) 



SRSF1 <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\exon_length\\SRSF1_ENCSR066VOO_SE.MATS.JCEC.txt',head=TRUE)
SRSF1$length <- SRSF1$exonEnd-SRSF1$exonStart_0base
SRSF1_ <- filter(SRSF1,abs(SRSF1$IncLevelDifference)>=0.1 & SRSF1$PValue<0.05)
boxplot(SRSF1$length,SRSF1_$length,ylim=c(0,300),names=c('all','SRSF1'),col='lightskyblue')


d <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\exon_length\\rbp_length_score.txt',head=FALSE)
####
d_fast <- read.table('D:\\wanglab\\smallexon\\new_60\\02pol2\\fast_data_short.txt',head=TRUE)
d_slow <-read.table('D:\\wanglab\\smallexon\\new_60\\02pol2\\slow_data_short.txt',head=TRUE)
d_case <- read.table('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\some_file\\exon_name.txt',head=F)
d_aqr <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\exon_length\\AQR_ENCSR624OUI_SE.MATS.JCEC.txt',head=TRUE)
d_rbfox2 <-read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\exon_length\\RBFOX2_ENCSR336DFS_SE.MATS.JCEC.txt',head=TRUE)
d_long <- read.table('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\some_file\\longexon_name.txt',head=F)

d_aqr <-filter(d_aqr,abs(d_aqr$IncLevelDifference)>=0.1) 
d_aqr <-filter(d_aqr,abs(d_aqr$PValue)<0.05) 
d_aqr$exon_length <- d_aqr$exonEnd - d_aqr$exonStart_0base
#d_aqr <- filter(d_aqr,d_aqr$exon_length<=60)
d_aqr$exon_name <-  paste(d_aqr$chr,paste(d_aqr$exonStart_0base,d_aqr$exonEnd,sep='-'),sep=':')
d_rbfox2 <-filter(d_rbfox2,abs(d_rbfox2$IncLevelDifference)>=0.1) 
d_rbfox2<-filter(d_rbfox2,abs(d_rbfox2$PValue)<0.05) 
d_rbfox2$exon_length <- d_rbfox2$exonEnd - d_rbfox2$exonStart_0base
#d_rbfox2 <- filter(d_rbfox2,d_rbfox2$exon_length<=60)
d_rbfox2$exon_name <-  paste(d_rbfox2$chr,paste(d_rbfox2$exonStart_0base,d_rbfox2$exonEnd,sep='-'),sep=':')


library(VennDiagram)
library(RColorBrewer)
venn.diagram(list(A=d_case$V1,B=d_aqr$exon_name,C=d_rbfox2$exon_name), "D:\\wanglab\\smallexon\\new_60\\Venn_2set_simple.jpeg",category.names = c( "CASEs", "AQR", "RBFOX2"),fill = brewer.pal(3, "Set2"))
venn.diagram(list(A = d_fast$exon_name, B = d_slow$exon_name,C=d_case$V1,D=d_aqr$exon_name,E=d_rbfox2$exon_name), "D:\\wanglab\\smallexon\\new_60\\Venn.jpeg",category.names = c("slow" , "fast " , "CASEs", "AQR", "RBFOX2"),fill = brewer.pal(5, "Set2"))


temp1 <- d_aqr[which(d_aqr$exon_name %in% d_case$V1),]
temp2 <- d_aqr[which(d_aqr$exon_name %in% d_long$V1),]
length(temp1[temp1$IncLevelDifference<0,][,1])
length(temp1[temp1$IncLevelDifference>0,][,1])
length(temp2[temp1$IncLevelDifference<0,][,1])
length(temp2[temp1$IncLevelDifference>0,][,1])

a <- c(length(temp1[temp1$IncLevelDifference<0,][,1]),length(temp1[temp1$IncLevelDifference>0,][,1]),length(temp2[temp1$IncLevelDifference<0,][,1]),length(temp2[temp1$IncLevelDifference>0,][,1]))
b <- c('3~60','3~60','61+','61+')   
c <- c('PSI_decrease','PSI_increase','PSI_decrease','PSI_increase')
dd <- as.data.frame(cbind(a,b,c),stringsAsFactors = F)
colnames(dd) <- c('number','length','type')
dd$type <- factor(dd$type,levels=c('PSI_increase','PSI_decrease'))

temp1 <- d_rbfox2[which(d_rbfox2$exon_name %in% d_case$V1),]
temp2 <- d_rbfox2[which(d_rbfox2$exon_name %in% d_long$V1),]

a <- c(length(temp1[temp1$IncLevelDifference<0,][,1]),length(temp1[temp1$IncLevelDifference>0,][,1]),length(temp2[temp1$IncLevelDifference<0,][,1]),length(temp2[temp1$IncLevelDifference>0,][,1]))
b <- c('3~60','3~60','61+','61+')   
c <- c('PSI_decrease','PSI_increase','PSI_decrease','PSI_increase')
dd <- as.data.frame(cbind(a,b,c),stringsAsFactors = F)
colnames(dd) <- c('number','length','type')
dd$type <- factor(dd$type,levels=c('PSI_increase','PSI_decrease'))


ggplot(dd,aes(y=as.numeric(dd$number),x=dd$length,fill=dd$type))+geom_bar(stat = "identity", position = "fill",width=0.7)+ geom_text(mapping = aes(label = number), size = 3.5,vjust = 3.5, hjust = .5, colour = 'black', position = position_fill())+theme_bw()+scale_fill_manual(name='',values = c("violet", "skyblue2"))+ylab("")+xlab("exon_length")+theme(panel.grid=element_blank()) 


###

a <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\RBFOX2ko_PCR\\PSI_pcr.txt',head=TRUE)
a.mean<-aggregate(as.numeric(a$PSI),by=list(a$class,a$gene),FUN=mean)
a.sd<-aggregate(as.numeric(a$PSI),by=list(a$class,a$gene),FUN=sd)
aa<-data.frame(a.mean,sd=a.sd$x)
#aa$Group.1<- factor(aa$Group.1,levels=c('WT','MUT'))
#UPF3B <- ggplot(data=aa,aes(x=Group.1,y=x,fill=Group.1))+geom_bar(stat="identity",position="dodge",width=0.6)+geom_errorbar(aes(ymax=x+sd,ymin=x-sd),position=position_dodge(0.9),width=0.15)+scale_fill_brewer(palette="Set1")+theme_classic()+xlab('UPF3B')+ylab('PSI')+scale_fill_manual(values=c('gray40','hotpink3'))+ guides(fill=FALSE)


aa1 <- filter(aa,aa$Group.2=='ARHGAP17')
ARHGAP17 <- ggplot(data=aa1,aes(x=Group.2,y=x,fill=Group.1))+geom_bar(stat="identity",position=position_dodge(0.8),width=0.7)+geom_errorbar(aes(ymax=x+sd,ymin=x-sd),position=position_dodge(0.8),width=0.15)+theme_classic()+xlab('')+ylab('PSI')+scale_fill_manual(values=c('gray37','coral3','coral2'))+ guides(fill=guide_legend(title=NULL))

aa2 <- filter(aa,aa$Group.2=='PTBP2')
PTBP2 <- ggplot(data=aa2,aes(x=Group.2,y=x,fill=Group.1))+geom_bar(stat="identity",position=position_dodge(0.8),width=0.7)+geom_errorbar(aes(ymax=x+sd,ymin=x-sd),position=position_dodge(0.8),width=0.15)+theme_classic()+xlab('')+ylab('PSI')+scale_fill_manual(values=c('gray37','coral3','coral2'))+ guides(fill=guide_legend(title=NULL))

aa1 <- filter(aa,aa$Group.2=='ARHGAP17')
ARHGAP17 <- ggplot(data=aa1,aes(x=Group.2,y=x,fill=Group.1))+geom_bar(stat="identity",position=position_dodge(0.8),width=0.7)+geom_errorbar(aes(ymax=x+sd,ymin=x-sd),position=position_dodge(0.8),width=0.15)+theme_classic()+xlab('')+ylab('PSI')+scale_fill_manual(values=c('gray37','coral3','coral2'))+ guides(fill=guide_legend(title=NULL))

aa3 <- filter(aa,aa$Group.2=='GNAS')
GNAS <- ggplot(data=aa3,aes(x=Group.2,y=x,fill=Group.1))+geom_bar(stat="identity",position=position_dodge(0.8),width=0.7)+geom_errorbar(aes(ymax=x+sd,ymin=x-sd),position=position_dodge(0.8),width=0.15)+theme_classic()+xlab('')+ylab('PSI')+scale_fill_manual(values=c('gray37','coral3','coral2'))+ guides(fill=guide_legend(title=NULL))


ggpubr::ggarrange(ARHGAP17,PTBP2,GNAS, nrow =1, ncol = 3)


############
data <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\Ï¸°û\\RBP.txt',head=TRUE)
ggplot(data=data, mapping=aes(x = CASE, y = PSI,fill=RBP))+ geom_boxplot(aes(fill = RBP),position=position_dodge(0.5),width=0.6)

AQR <- filter(data,data$TYPE=='AQR')
ggplot(data=AQR, mapping=aes(x = CASE, y = PSI,fill=RBP))+ geom_boxplot(aes(fill = RBP),position=position_dodge(0.5),width=0.6)
aqr <- ggplot(data=AQR, mapping=aes(x = CASE, y = PSI,fill=RBP))+ geom_boxplot(aes(fill = RBP),position=position_dodge(0.5),width=0.6)+scale_fill_manual(values = c('grey80', 'tomato3'))+theme_bw()+theme(panel.grid=element_blank())+stat_compare_means(aes(group=RBP), label = "p.signif",label.y = 0.5,method='t.test') 


U2AF2 <- filter(data,data$TYPE=='U2AF2')
u2af2 <- ggplot(data=U2AF2, mapping=aes(x = CASE, y = PSI,fill=RBP))+ geom_boxplot(aes(fill = RBP),position=position_dodge(0.5),width=0.6)+scale_fill_manual(values = c('grey80', 'tomato3'))+theme_bw()+theme(panel.grid=element_blank()) +stat_compare_means(aes(group=RBP), label = "p.signif",label.y = 0.55,method='t.test') 

PTBP1 <- filter(data,data$TYPE=='PTBP1')
ptbp1 <- ggplot(data=PTBP1, mapping=aes(x = CASE, y = PSI,fill=RBP))+ geom_boxplot(aes(fill = RBP),position=position_dodge(0.5),width=0.6)+scale_fill_manual(values = c('grey80', 'tomato3'))+theme_bw()+theme(panel.grid=element_blank()) +stat_compare_means(aes(group=RBP), label = "p.signif",label.y = 1,method='t.test') 


RBFOX2<- filter(data,data$TYPE=='RBFOX2')
rbfox2 <- ggplot(data=RBFOX2, mapping=aes(x = CASE, y = PSI,fill=RBP))+ geom_boxplot(aes(fill = RBP),position=position_dodge(0.5),width=0.6)+scale_fill_manual(values = c('grey80', 'tomato3'))+theme_bw()+theme(panel.grid=element_blank()) +stat_compare_means(aes(group=RBP), label = "p.signif",label.y = 0.5,method='t.test') 


ggpubr::ggarrange(rbfox2,aqr,u2af2,ptbp1, nrow =2, ncol = 2)





PTBP2<- filter(data,data$CASE=='PTBP2')
ptbp2 <- ggplot(data=PTBP2, mapping=aes(x =TYPE, y = PSI,fill=TYPE_2))+ ggtitle('PTBP2')+geom_boxplot(aes(fill = TYPE_2),position=position_dodge(0.7),width=0.7)+geom_jitter(position=position_dodge(0.7),aes(fill= TYPE_2),shape = 21,size=2)+scale_fill_manual(values = c('grey80', 'tomato3'))+theme_bw()+theme(panel.grid=element_blank())+theme_classic()+ guides(fill=guide_legend(title=NULL))+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+stat_compare_means(aes(group=TYPE_2), label = "p.signif",label.y = 1.1,method='t.test')+xlab(" ")

ARHGAP17 <- filter(data,data$CASE=='ARHGAP17')
arhgap17 <- ggplot(data=ARHGAP17, mapping=aes(x =TYPE, y = PSI,fill=TYPE_2))+ ggtitle('ARHGAP17')+ geom_boxplot(aes(fill = TYPE_2),position=position_dodge(0.7),width=0.7)+geom_jitter(position=position_dodge(0.7),aes(fill= TYPE_2),shape = 21,size=2)+scale_fill_manual(values = c('grey80', 'tomato3'))+theme_bw()+theme(panel.grid=element_blank())+theme_classic()+ guides(fill=guide_legend(title=NULL))+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+stat_compare_means(aes(group=TYPE_2), label = "p.signif",label.y = 0.35,method='t.test')+xlab(" ")


MYL6<- filter(data,data$CASE=='MYL6')
myl6 <- ggplot(data=MYL6, mapping=aes(x =TYPE, y = PSI,fill=TYPE_2))+ ggtitle('MYL6')+ geom_boxplot(aes(fill = TYPE_2),position=position_dodge(0.7),width=0.7)+geom_jitter(position=position_dodge(0.7),aes(fill= TYPE_2),shape = 21,size=2)+scale_fill_manual(values = c('grey80', 'tomato3'))+theme_bw()+theme(panel.grid=element_blank())+theme_classic()+ guides(fill=guide_legend(title=NULL))+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+stat_compare_means(aes(group=TYPE_2), label = "p.signif",label.y = 0.5,method='t.test')+xlab(" ")


GNAS<- filter(data,data$CASE=='GNAS')
gnas <- ggplot(data=GNAS, mapping=aes(x =TYPE, y = PSI,fill=TYPE_2))+ ggtitle('GNAS')+ geom_boxplot(aes(fill = TYPE_2),position=position_dodge(0.7),width=0.7)+geom_jitter(position=position_dodge(0.7),aes(fill= TYPE_2),shape = 21,size=2)+scale_fill_manual(values = c('grey80', 'tomato3'))+theme_bw()+theme(panel.grid=element_blank())+theme_classic()+ guides(fill=guide_legend(title=NULL))+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+stat_compare_means(aes(group=TYPE_2), label = "p.signif",label.y = 0.6,method='t.test')+xlab(" ")


RPS24<- filter(data,data$CASE=='RPS24')
rps24 <- ggplot(data=RPS24, mapping=aes(x =TYPE, y = PSI,fill=TYPE_2))+ ggtitle('RPS24')+ geom_boxplot(aes(fill = TYPE_2),position=position_dodge(0.7),width=0.7)+geom_jitter(position=position_dodge(0.7),aes(fill= TYPE_2),shape = 21,size=2)+scale_fill_manual(values = c('grey80', 'tomato3'))+theme_bw()+theme(panel.grid=element_blank())+theme_classic()+ guides(fill=guide_legend(title=NULL))+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+stat_compare_means(aes(group=TYPE_2), label = "p.signif",label.y = 0.65,method='t.test')+xlab(" ")


ggpubr::ggarrange(ptbp2,rps24,gnas,myl6,arhgap17, nrow =2, ncol = 3)

####
data <- read.table('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\some_file\\all_merge.txt',head=TRUE)
f_na<-function(x) sum(is.na(x))
na <- apply(data,1,f_na)
data$NA_number <- na
library(dplyr)
data <- filter(data,data$NA_number<=5000)
write.table(data,file='D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\some_file\\all_merge_na.txt',row.names=F,sep='\t')

data_2 <- read.table('D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\some_file\\all_merge_na_t.txt',head=TRUE)
data<-merge(data_2,SF1,by='name')
x <- c()
y <- c()
for (i in 2:4470){
  x <- rbind(x,(colnames(data)[i]))
  y <- rbind(y,(cor(data$`RBM9|23543`,data[,i],use = 'complete.obs')))
}
rownames(y)<- x
y <- as.data.frame(y)
write.table(y,file='D:\\wanglab\\smallexon\\new_60\\03rbp\\cor_rbfox2.txt',row.names=T,sep='\t')

data <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\motif_rbfox2_cor.txt',head=TRUE)
data2 <- read.table('D:\\wanglab\\smallexon\\new_60\\01features\\exon_class\\allexon_class.txt',head=TRUE)
library(ggplot2)
data <- filter(data,data$class != 'NA')
c <- ggplot(data, aes(x = cor))+ geom_density(aes(fill = class), alpha=0.4)+theme_bw()+theme(panel.grid=element_blank()) +xlab('correlation')


data3 <- merge(data,data2,by='exon')

data3$type <- factor(data3$type,levels=c('long_other','long_cancer','short_other','short_cancer'))
a <-ggplot(data3, aes(x = cor))+ geom_density(aes(fill =type), alpha=0.6)+theme_bw()+theme(panel.grid=element_blank()) +xlab('correlation')+scale_fill_manual(values=c("gold3","#EC7357","#9DC3C1","forestgreen"))
data4 <- filter(data3,data3$type=='short_cancer')

b<- ggplot(data4, aes(x = cor))+ geom_density(aes(fill = class), alpha=0.4)+theme_bw()+theme(panel.grid=element_blank()) +xlab('correlation')



####
setwd('D:\\wanglab\\smallexon\\new_60\\01features\\014splicesite')
short_c <- read.table('short_cancer_3ss_score.txt')
short_o <- read.table('short_other_3ss_score.txt')
long_c <- read.table('long_cancer_3ss_score.txt')
long_o <- read.table('long_other_3ss_score.txt')

z_3 <- rbind(short_c,short_o,long_c,long_o)
u2af_3 <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\exon_length\\splicesite\\U2AF2_3ss_score.txt')
ptbp_3 <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\exon_length\\splicesite\\PTBP1_3ss_score.txt')
rbfox_3 <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\exon_length\\splicesite\\RBFOX2_3ss_score.txt')
aqr_3 <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\exon_length\\splicesite\\AQR_3ss_score.txt')

type <- c(rep('all',189930),rep('U2AF2',4931),rep('PTBP1',3793),rep('AQR',8412),rep('RBFOX2',1713))
z <- as.data.frame(rbind(z_3,u2af_3,ptbp_3,aqr_3,rbfox_3))
colnames(z)<- c('sequence','3SS_score')
z$type <- type
colnames(z)
ss3 <- ggplot(z,aes(y=z$`3SS_score`,x=z$type)) + geom_boxplot(width=0.8,fill='lightpink')+theme_bw()+xlab("")+ylab("3' splice site score")+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+theme(panel.grid=element_blank())+ylim(-25,20)


short_c <- read.table('short_cancer_5ss_score.txt')
short_o <- read.table('short_other_5ss_score.txt')
long_c <- read.table('long_cancer_5ss_score.txt')
long_o <- read.table('long_other_5ss_score.txt')

z_5 <- rbind(short_c,short_o,long_c,long_o)
u2af_5 <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\exon_length\\splicesite\\U2AF2_5ss_score.txt')
ptbp_5 <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\exon_length\\splicesite\\PTBP1_5ss_score.txt')
rbfox_5 <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\exon_length\\splicesite\\RBFOX2_5ss_score.txt')
aqr_5 <- read.table('D:\\wanglab\\smallexon\\new_60\\03rbp\\exon_length\\splicesite\\AQR_5ss_score.txt')

type <- c(rep('all',189930),rep('U2AF2',4931),rep('PTBP1',3793),rep('AQR',8412),rep('RBFOX2',1713))
z5 <- as.data.frame(rbind(z_5,u2af_5,ptbp_5,aqr_5,rbfox_5))
colnames(z5)<- c('sequence','3SS_score')
z5$type <- type
colnames(z5)
ss5 <- ggplot(z5,aes(y=z5$`3SS_score`,x=z5$type)) + geom_boxplot(width=0.8,fill='lightpink')+theme_bw()+xlab("")+ylab("5' splice site score")+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+theme(panel.grid=element_blank())+ylim(-25,20)
