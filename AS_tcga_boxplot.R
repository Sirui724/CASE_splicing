####TCGA AS PSI boxplot

a <- 'chr3:152163071:152163328:+@chr3:152164493:152164546:+@chr3:152165409:152165562:+'


BLCA <-read.table('D:\\wanglab\\smallexon\\TCGA_smallexon\\calculated_data_summary\\calculated_BLCA.txt',sep="\t",row.names = 2,head=TRUE)
cancer_type <- rep('BLCA',length(BLCA)-5)
d_temp <- as.data.frame(cbind(t(BLCA[a,])[-c(1,2,3,4,5),],cancer_type))


BRCA <-read.table('D:\\wanglab\\smallexon\\TCGA_smallexon\\calculated_data_summary\\calculated_BRCA.txt',sep="\t",row.names = 2,head=TRUE)
cancer_type <- rep('BRCA',length(BRCA)-5)
TEMP<- as.data.frame(cbind(t(BRCA[a,])[-c(1,2,3,4,5),],cancer_type))
d_temp<- rbind(d_temp,TEMP)

temp <-read.table('D:\\wanglab\\smallexon\\TCGA_smallexon\\calculated_data_summary\\calculated_COAD.txt',sep="\t",row.names = 2,head=TRUE)
cancer_type <- rep('COAD',length(temp)-5)
TEMP<- as.data.frame(cbind(t(temp[a,])[-c(1,2,3,4,5),],cancer_type))
d_temp<- rbind(d_temp,TEMP)

temp <-read.table('D:\\wanglab\\smallexon\\TCGA_smallexon\\calculated_data_summary\\calculated_ESCA.txt',sep="\t",row.names = 2,head=TRUE)
cancer_type <- rep('ESCA',length(temp)-5)
TEMP<- as.data.frame(cbind(t(temp[a,])[-c(1,2,3,4,5),],cancer_type))
d_temp<- rbind(d_temp,TEMP)


temp <-read.table('D:\\wanglab\\smallexon\\TCGA_smallexon\\calculated_data_summary\\calculated_GBM.txt',sep="\t",row.names = 2,head=TRUE)
cancer_type <- rep('GBM',length(temp)-5)
TEMP<- as.data.frame(cbind(t(temp[a,])[-c(1,2,3,4,5),],cancer_type))
d_temp<- rbind(d_temp,TEMP)


temp <-read.table('D:\\wanglab\\smallexon\\TCGA_smallexon\\calculated_data_summary\\calculated_HNSC.txt',sep="\t",row.names = 2,head=TRUE)
cancer_type <- rep('HNSC',length(temp)-5)
TEMP<- as.data.frame(cbind(t(temp[a,])[-c(1,2,3,4,5),],cancer_type))
d_temp<- rbind(d_temp,TEMP)

temp <-read.table('D:\\wanglab\\smallexon\\TCGA_smallexon\\calculated_data_summary\\calculated_KICH.txt',sep="\t",row.names = 2,head=TRUE)
cancer_type <- rep('KICH',length(temp)-5)
TEMP<- as.data.frame(cbind(t(temp[a,])[-c(1,2,3,4,5),],cancer_type))
d_temp<- rbind(d_temp,TEMP)

temp <-read.table('D:\\wanglab\\smallexon\\TCGA_smallexon\\calculated_data_summary\\calculated_KIRC.txt',sep="\t",row.names = 2,head=TRUE)
cancer_type <- rep('KIRC',length(temp)-5)
TEMP<- as.data.frame(cbind(t(temp[a,])[-c(1,2,3,4,5),],cancer_type))
d_temp<- rbind(d_temp,TEMP)

temp <-read.table('D:\\wanglab\\smallexon\\TCGA_smallexon\\calculated_data_summary\\calculated_KIRP.txt',sep="\t",row.names = 2,head=TRUE)
cancer_type <- rep('KIRP',length(temp)-5)
TEMP<- as.data.frame(cbind(t(temp[a,])[-c(1,2,3,4,5),],cancer_type))
d_temp<- rbind(d_temp,TEMP)

temp <-read.table('D:\\wanglab\\smallexon\\TCGA_smallexon\\calculated_data_summary\\calculated_LIHC.txt',sep="\t",row.names = 2,head=TRUE)
cancer_type <- rep('LIHC',length(temp)-5)
TEMP<- as.data.frame(cbind(t(temp[a,])[-c(1,2,3,4,5),],cancer_type))
d_temp<- rbind(d_temp,TEMP)

temp <-read.table('D:\\wanglab\\smallexon\\TCGA_smallexon\\calculated_data_summary\\calculated_LUAD.txt',sep="\t",row.names = 2,head=TRUE)
cancer_type <- rep('LUAD',length(temp)-5)
TEMP<- as.data.frame(cbind(t(temp[a,])[-c(1,2,3,4,5),],cancer_type))
d_temp<- rbind(d_temp,TEMP)

temp <-read.table('D:\\wanglab\\smallexon\\TCGA_smallexon\\calculated_data_summary\\calculated_LUSC.txt',sep="\t",row.names = 2,head=TRUE)
cancer_type <- rep('LUSC',length(temp)-5)
TEMP<- as.data.frame(cbind(t(temp[a,])[-c(1,2,3,4,5),],cancer_type))
d_temp<- rbind(d_temp,TEMP)

temp <-read.table('D:\\wanglab\\smallexon\\TCGA_smallexon\\calculated_data_summary\\calculated_PAAD.txt',sep="\t",row.names = 2,head=TRUE)
cancer_type <- rep('PAAD',length(temp)-5)
TEMP<- as.data.frame(cbind(t(temp[a,])[-c(1,2,3,4,5),],cancer_type))
d_temp<- rbind(d_temp,TEMP)

temp <-read.table('D:\\wanglab\\smallexon\\TCGA_smallexon\\calculated_data_summary\\calculated_PRAD.txt',sep="\t",row.names = 2,head=TRUE)
cancer_type <- rep('PRAD',length(temp)-5)
TEMP<- as.data.frame(cbind(t(temp[a,])[-c(1,2,3,4,5),],cancer_type))
d_temp<- rbind(d_temp,TEMP)

temp <-read.table('D:\\wanglab\\smallexon\\TCGA_smallexon\\calculated_data_summary\\calculated_READ.txt',sep="\t",row.names = 2,head=TRUE)
cancer_type <- rep('READ',length(temp)-5)
TEMP<- as.data.frame(cbind(t(temp[a,])[-c(1,2,3,4,5),],cancer_type))
d_temp<- rbind(d_temp,TEMP)

temp <-read.table('D:\\wanglab\\smallexon\\TCGA_smallexon\\calculated_data_summary\\calculated_STAD.txt',sep="\t",row.names = 2,head=TRUE)
cancer_type <- rep('STAD',length(temp)-5)
TEMP<- as.data.frame(cbind(t(temp[a,])[-c(1,2,3,4,5),],cancer_type))
d_temp<- rbind(d_temp,TEMP)

temp <-read.table('D:\\wanglab\\smallexon\\TCGA_smallexon\\calculated_data_summary\\calculated_THCA.txt',sep="\t",row.names = 2,head=TRUE)
cancer_type <- rep('THCA',length(temp)-5)
TEMP<- as.data.frame(cbind(t(temp[a,])[-c(1,2,3,4,5),],cancer_type))
d_temp<- rbind(d_temp,TEMP)

temp <-read.table('D:\\wanglab\\smallexon\\TCGA_smallexon\\calculated_data_summary\\calculated_UCEC.txt',sep="\t",row.names = 2,head=TRUE)
cancer_type <- rep('UCEC',length(temp)-5)
TEMP<- as.data.frame(cbind(t(temp[a,])[-c(1,2,3,4,5),],cancer_type))
d_temp<- rbind(d_temp,TEMP)


library('stringr')
library('dplyr')
tumor<-c()
normal<-c()
type <- c()
for (i in seq(from=1,to=7528)){
  if (str_sub(rownames(d_temp)[i],start=14L,end=15L)=='01' | str_sub(rownames(d_temp)[i],start=14L,end=15L)=='02'){
    tumor<- c(tumor,i)
    type <-c(type,'tumor')
  }
  else {
    normal<-c(normal,i)
    type <-c(type,'normal')
  }
}
#if (str_sub(rownames(SF1)[i],start=14L,end=15L)=='11')
d_temp<- cbind(d_temp,type)

ggplot(d_temp,aes(y=as.numeric(as.character(d_temp$V1)),x=d_temp$cancer_type,fill=type)) + geom_boxplot()+
  scale_fill_manual(values = c('grey80', 'tomato3'))+
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.title = element_blank(), legend.key = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  labs(x = '', y = 'MBNL1_PSI')

ggsave("D:\\wanglab\\smallexon\\new_60\\00tcga_shortexon\\some_plots\\ARHGAP17.pdf", width = 13, height =7.5, units = "cm")

