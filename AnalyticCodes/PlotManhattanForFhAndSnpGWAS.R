# plot fh / snp GWAS
setwd("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Seprate/")

library(data.table)
pheList = list.files("E:/QiGA/document/CottonLab/XLD/Phe/EMMAX/Seprated/Phe/","EMMAX.*.txt")
#PHE = sub("EMMAX_","",sapply(strsplit(pheList,".",fixed=T),'[[',1))
phe_field1 = sapply(strsplit(pheList,".",fixed=T),'[[',2)
phe_field2 = sapply(strsplit(pheList,".",fixed=T),'[[',3)
PHE = paste0(phe_field1,"_",phe_field2)
PHE1 = paste0(phe_field1,".",phe_field2)

gff = fread("E:/QiGA/document/CottonLab/Cotton_GO/TM-1_V2.1.gene.gff",header=F,data.table = F)
gff = gff[gff$V3=='mRNA',]
gff$V9 = sub(";","",sub("ID=","",gff$V9,fixed=T),fixed=T)
gff = gff[,c(1,4,5,7,9)]
colnames(gff) = c("CHR","START","END","STRAND","Gene")

#source("https://raw.githubusercontent.com/YinLiLin/CMplot/master/R/CMplot.r")
source("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/plot/cmplot.R")
#for(i in 1:length(PHE)){
for(i in 1:4){
  pheName = PHE[i]
  fh.glmm.File = paste0("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Seprate/Phe_",pheName,".GLMM.csv")
  snp.glmm.file = paste0("E:/QiGA/document/CottonLab/XLD/Phe/EMMAX/Seprated/Regression/", PHE1[i], ".ps.gz")
  if(!file.exists(fh.glmm.File)){
    cat(paste0("GLMM results for trait ",pheName," not found.\n"))
    stop()
  }
  
  fh.glmm = fread(file = fh.glmm.File,data.table = F,sep=',',header = T,stringsAsFactors = F)
  snp.glmm = fread(file = snp.glmm.file,data.table = F,header = F,sep='\t',stringsAsFactors = F,encoding = 'UTF-8')
  
  fh.glmm = merge(fh.glmm,gff,'Gene')
  fh.glmm$P = as.numeric(fh.glmm$P)
  fh.glmm$Chisq = as.numeric(fh.glmm$Chisq)
  fh.glmm$Padj = pchisq(fh.glmm$Chisq/mean(fh.glmm$Chisq,na.rm=T),1,0,lower.tail = F)
  colnames(fh.glmm)[1] = 'Variants'
  # if(i < 5){
  #   fh.glmm.plot = fh.glmm[,c(1,15,16,19)] # p.adj
  # } else {
  #   fh.glmm.plot = fh.glmm[,c(1,15,16,14)] # p
  # }
  fh.glmm.plot = fh.glmm[,c(1,15,16,14)] # p
  colnames(fh.glmm.plot)[4] = 'P'
  fh.glmm.plot$P[which(fh.glmm.plot$P=='HighVaried')] = 1
  fh.glmm.plot$P = as.numeric(fh.glmm.plot$P)
  
  snp.glmm = na.omit(snp.glmm)
  colnames(snp.glmm) = c("Variants","Beta","SE","P-SNP")
  snp.glmm$CHR = sapply(strsplit(snp.glmm$Variants,'_',fixed=T),'[[',1)
  snp.glmm$POS = sapply(strsplit(snp.glmm$Variants,'_',fixed=T),'[[',2)
  snp.glmm$Ref_Allele = sapply(strsplit(snp.glmm$Variants,'_',fixed=T),'[[',3)
  snp.glmm$Alt_Allele = sapply(strsplit(snp.glmm$Variants,'_',fixed=T),'[[',4)
  snp.glmm.plot = snp.glmm[,c(1,5,6,4)]
  
  colnames(snp.glmm.plot) = colnames(fh.glmm.plot) = c("SNP","Chromosome","Position",'P')
  
  # plot2
  #fh.glmm.plot$P = p.adjust(fh.glmm.plot$P)
  plotDT = rbind(snp.glmm.plot,fh.glmm.plot)
  colnames(plotDT)[4] = 'SNP-P'
  plotDT$`FH-P`= NA
  plotDT$`FH-P`[grep("GH",plotDT$SNP)] = plotDT$`SNP-P`[grep("GH",plotDT$SNP)]
  plotDT$`SNP-P`[grep("GH",plotDT$SNP)] = NA
  plotDT = plotDT[,c(1,2,3,5,4)]
  cmplot(plotDT, plot.type="m", multraits=TRUE,cex = c(1.1,1), col = c("#FFC20A","#5DB1DD"),mar = c(0.5,1,1,1),
         LOG10=TRUE, threshold=c(1e-5),axis.cex = 1.3,legend.pos = 'none',
         threshold.lty=2, threshold.lwd=1, threshold.col="darkblue", amplify=F, pch = c(16,16),
         chr.den.col=NULL, points.alpha = 120,chr.labels.angle = 45,
         file="jpg",file.name=paste0(pheName,".croped_2"),dpi=600,file.output=TRUE,verbose=TRUE,width=16,height=3)
  cmplot(plotDT, plot.type="m", multraits=TRUE,cex = c(1.1,1), col = c("#FFC20A","#5DB1DD"),
         LOG10=TRUE, threshold=c(1e-5),axis.cex = 1.3,legend.cex = 1.3,
         threshold.lty=2, threshold.lwd=1, threshold.col="darkblue", amplify=F,pch = c(16,16),
         chr.den.col=NULL, points.alpha = 120,chr.labels.angle = 45,
         file="jpg",file.name=paste0(pheName,"_2"),dpi=600,file.output=TRUE,verbose=TRUE,width=16,height=3)
  
  # valid.snp = nrow(snp.glmm)
  # threshold = 0.05/valid.snp
  # 
  # sig.snp.glmm = snp.glmm[snp.glmm$V4<=threshold,]
  # colnames(sig.snp.glmm) = c("Variants","Beta","SE","P-SNP")
}
rm(list = ls())
