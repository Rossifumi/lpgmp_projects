######### read in fpkm
# t72<-read.table('~/WORK/projects/LR_GN/diceseq_exp/72h_hisat_sort.count.txt',header = T,sep = "\t")
# t72<-t72[,-(2:3)]
# t72$gene_length=72
# colnames(t72)<-c("gene_id","sampling_time","count","FPKM")
# 
# Bna_fpkm<-rbind(t0,t0_5,t1,t2,t3,t4,t5,t6,t8,t10,t12,t14,t16,t20,t24,t28,t32,t36,t40,t44,t48,t56,t64,t72)
#> dim(Bna_fpkm)
#[1] 2424960       4

# savelist<-c("savelist","Bna_fpkm")
# save(list = savelist,file = "~/WORK/projects/LR_GN/diceseq_exp/Bna_exp_profile.RData")

load(file = "~/WORK/projects/LR_GN/diceseq_exp/Bna_exp_profile.RData")

######### plot expression pattern for certain gene
library(ggplot2)

Bna_fpkm[Bna_fpkm$gene_id=="BnaA02g05070D",]

######### read in qPCR delta delta cq value

ddcq<-read.table(file="/Users/Jiajia/WORK/projects/LR_GN/qPCR/DeltaDeltaCq_value.txt",head=T,sep="\t")

# ddcq[ddcq$Gene=="ARF7-A1-2",]->tmp

library(plotrix)
# rescale(tmp[tmp$Repeat=="A",]$Value_mod,range(Bna_fpkm[Bna_fpkm$gene_id=="BnaA02g05070D",]$FPKM))
# plot(Bna_fpkm[Bna_fpkm$gene_id=="BnaA02g05070D",]$sampling_time,Bna_fpkm[Bna_fpkm$gene_id=="BnaA02g05070D",]$FPKM,type="l")
# lines(Bna_fpkm[Bna_fpkm$gene_id=="BnaA02g05070D",]$sampling_time,rescale(tmp[tmp$Repeat=="A",]$Value_mod,range(Bna_fpkm[Bna_fpkm$gene_id=="BnaA02g05070D",]$FPKM)))

######### compare fpkm with cq, cor test, plot for 16 genes

gene_names<-levels(ddcq$Gene)
gene_ids<-c("BnaC05g14880D","BnaA02g05070D","BnaA10g14760D","BnaC09g37130D","BnaC03g20410D","BnaAnng12210D","BnaC08g50150D","BnaA04g13630D","BnaC04g35880D","BnaA09g37320D","BnaC08g29120D","BnaA06g14040D","BnaC05g15390D","BnaC02g44180D")

pdf(file="/Users/Jiajia/WORK/projects/LR_GN/qPCR/qPCR_cor_fpkm_cq.pdf")
for (i in 1:length(gene_names))
{
  g_fpkm = Bna_fpkm[Bna_fpkm$gene_id==gene_ids[i],]$FPKM
  g_cq_a = rescale(subset(ddcq,ddcq$Repeat=="A" & ddcq$Gene==gene_names[i])$Value_mod,range(g_fpkm))
  g_cq_b = rescale(subset(ddcq,ddcq$Repeat=="B" & ddcq$Gene==gene_names[i])$Value_mod,range(g_fpkm))
  g_cq_c = rescale(subset(ddcq,ddcq$Repeat=="C" & ddcq$Gene==gene_names[i])$Value_mod,range(g_fpkm))
  g_cq<-(g_cq_a+g_cq_b+g_cq_c)/3
  g_p_value = format(cor.test(g_fpkm,g_cq)$p.value,scientific = TRUE,digits = 3)
  g_cor = round(cor(g_fpkm,g_cq),digits = 3)
  
  plot(x=Bna_fpkm[Bna_fpkm$gene_id==gene_ids[i],]$sampling_time,y=g_fpkm,type="l",lwd=3,main=paste(gene_names[i]," cor = ",g_cor," p-value = ",g_p_value,sep=""),xlab="Time(h)",ylab="Expression Level")
  lines(x=Bna_fpkm[Bna_fpkm$gene_id==gene_ids[i],]$sampling_time,y=g_cq,type="l",lwd=3,col="grey")
}
dev.off()

### fetch data by SQL

library(sqldf)

# sqldf("select Time,Value_mod from ddcq where Gene='AT2G36930-C1-1' and Repeat='A'")


######### LR gene list

LR_list<-read.table(file="/Users/Jiajia/WORK/projects/LR_GN/exp_profile/LR_GRN_Bna.list.sort.append.tab",head=F,sep="\t")

# savelist<-c("savelist","Bna_fpkm","LR_list")
# save(list = savelist,file = "~/WORK/projects/LR_GN/diceseq_exp/Bna_exp_profile.RData")

Bna_fpkm[Bna_fpkm$gene_id %in% LR_list$V1,]->Bna_fpkm_LR


######### filter expression profile with certian criteria and then interpolation 

# python LR_GRN_interpolation_fpkm.py


######## k-means clustering 

source('/Users/Jiajia/WORK/projects/LR_GN/diceseq_exp/k-means_single_species_func.R')

Mat_Ori<-read.csv('/Users/Jiajia/WORK/projects/LR_GN/diceseq_exp/Bna_exp_fpkm_filter_interp.csv',header = TRUE,sep = ",")

a<-normalize.rows.f(Mat_Ori[,-1])

#30 is cluster number; 20 is iteration times; Rice_30line_c40 is desired file name
setwd('/Users/Jiajia/WORK/projects/LR_GN/diceseq_exp/Bna_fi_c30/')
k_means(a,30,50,"Bna_fi_c30",age=seq(0,72,by=0.5)) 


######## TDCor check time delay

library(TDCor)
library(varhandle)
Bna_mat<-read.table(file = "~/WORK/projects/LR_GN/diceseq_exp/Bna_exp_fpkm.tab",head=F,sep="\t")
Bna_LR_mat<-Bna_mat[Bna_mat[,1]%in%LR_list[,1],]
rownames(Bna_LR_mat)<-Bna_LR_mat[,1]
Bna_LR_mat<-Bna_LR_mat[,-1]
colnames(Bna_LR_mat)<-c(1:24)

times<-c(0,0.5,1,2,3,4,5,6,8,10,12,14,16,20,24,28,32,36,40,44,48,56,64,72)

as.numeric(unfactor(Bna_LR_mat[1,]))


