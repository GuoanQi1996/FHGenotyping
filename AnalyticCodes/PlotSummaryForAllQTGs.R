# FH hap analysis
library(data.table)
library(ggplot2)
library(ggplotify)
library(ggsci)
library(grid)
library(ggpubr)
library(cowplot)
library(scales)
library(agricolae)
library(Biostrings)
library(gridExtra)
source("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/plot/LocusZooms-2.1/functions/locus_zoom.R")

setwd("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM")

GT = fread("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/FHAP_whole.type.2.txt", sep = '\t', data.table = F, stringsAsFactors = F)
colnames(GT) = sub("0_","",colnames(GT))
GT1 = t(GT[,-1])
colnames(GT1) = GT$Gene
GT1 = data.frame(ID = rownames(GT1),GT1)

plotDT = read.csv("SignificantGene.20240909.csv")
unitDT = read.csv("Units.csv")

gff = fread("E:/QiGA/document/CottonLab/Cotton_GO/TM-1_V2.1.gene.gff",header=F,data.table = F)
gff = gff[gff$V3=='mRNA',]
gff$V9 = sub(";","",sub("ID=","",gff$V9,fixed=T),fixed=T)
gff = gff[,c(1,4,5,9)]
gff$Coding = 'proteincoding'
colnames(gff) = c("Chrom","Start","End","Gene","Coding")

Bims = fread("E:/QiGA/document/CottonLab/FH/data/core_Gh_accession.union.MAF0001GENO20.bim",header=F,stringsAsFactors = F,data.table = F,nThread = 14)

name2paper = data.frame(Name = c("DXM1202", "HZG486", "MZY1081", "MZY419"),
                        Paper = c("ICRCAAS_NG_2021.N1202","ZJU_ICP_2022.N486","HEBAU_NG_2021.N1081","HEBAU_NG_2018.N419"))
empty_df = data.frame(x = 1, y = 1, label = "Except for the reference sequence,\nno FH appears in more than 2 samples/lines in the population\n(Find more details in the FH summary of this gene)")
#CDS = readDNAStringSet("E:/QiGA/document/CottonLab/Cotton_GO/TM-1_V2.1.gene.cds.fa")
for(i in 1:nrow(plotDT)){
#for(i in 1:10){
  gene = plotDT$Gene[i]
  phe_short = strsplit(plotDT$Phe[i],"_",fixed=T)[[1]][2]
  phe_labs = unitDT$Units[which(unitDT$Phe1==phe_short)]
  
  dataset = strsplit(plotDT$Phe[i],"_",fixed=T)[[1]][1]
  outputName = name2paper$Paper[which(name2paper$Name==dataset)]
  field1 = strsplit(plotDT$Phe[i],paste0(dataset,"_"),fixed=T)[[1]][2]
  outputName = paste0(outputName,".",field1)
  dataset = paste0(dataset,".",field1)
  manhattan_names = sub(".","_",dataset,fixed=T)
  
  print(paste0(i,"/27991."))
  
  if(file.exists(paste0("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/CombinedFigures/FH.", gene,".", outputName, ".png"))){
    next
  }

  p_manhat = jpeg::readJPEG(paste0("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Seprate/Multi-traits_Manhtn.", manhattan_names, "_2.jpg"))
  p_manhat = ggplot()+background_image(p_manhat)+theme_void()

  phe = fread(paste0("E:/QiGA/document/CottonLab/XLD/Phe/Seprate_phenotype/Phe/EMMAX.",dataset,".txt"),header=F,sep = '\t',data.table = F)
  #phe = phe[,c("Sample","PH_2016_SHZ")]
  colnames(phe) = c("FID","Sample","Phenotype")
  #colnames(phe) = c("Sample","Phenotype")
  #phe$Sample = paste0("HZG486_",phe$Sample)
  
  targetFH = GT1[,gene]
  SampleNames = rownames(GT1)
  DT = data.frame(Sample = SampleNames, FH = targetFH)
  DT = merge(DT,phe,"Sample")
  DT = na.omit(DT)
  
  DTStat = as.data.frame(table(DT$FH))
  DTStat$Var1 = as.character(DTStat$Var1)
  lowDiversity = as.character(DTStat$Var1[which(DTStat$Freq<2)])
  rmFH = which(!is.na(match(as.character(DT$FH),lowDiversity)))
  if(length(rmFH)>0){
    DT = DT[-rmFH,]
    DTStat = DTStat[-which(!is.na(match(as.character(DTStat$Var1),lowDiversity))),]
  }
  
  DTAve = aggregate(DT$Phenotype,list(DT$FH),mean,na.rm=T)
  DTMedi = aggregate(DT$Phenotype,list(DT$FH),median,na.rm=T)
  DTQuantile = aggregate(DT$Phenotype,list(DT$FH),function(x){
    iqr = IQR(x)
    Q3 = quantile(x,probs=0.75)
    min(max(x), Q3 + 1.5 * iqr)
  })
  colnames(DTAve)[2] = "average"
  colnames(DTMedi)[2] = "median"
  colnames(DTQuantile)[2] = "quantile"
  DTStat = merge(DTStat,DTMedi,by.x = 'Var1',by.y = 'Group.1')
  DTStat = merge(DTStat,DTAve,by.x = 'Var1',by.y = 'Group.1')
  DTStat = merge(DTStat,DTQuantile,by.x = 'Var1',by.y = 'Group.1')
  DTStat = DTStat[order(DTStat$average,decreasing = F),]
  xlabel = paste0("FH-",DTStat$Var1,"\n(n=",DTStat$Freq,")")
  DT$FH = factor(DT$FH,levels = as.character(DTStat$Var1),labels = xlabel)
  DTStat$Var1 = factor(DTStat$Var1,levels = as.character(DTStat$Var1),labels = xlabel)
  
  if(nrow(DTStat)==1){
    p <- ggplot(empty_df, aes(x = x, y = y)) +
      geom_rect(aes(xmin = 0, xmax = 2, ymin = 0, ymax = 2), fill = "transparent") +
      geom_text(aes(label = label), color = "black",size = 3) +
      theme_void()
  } else {
    variance = aov(Phenotype~FH,data = DT)
    MC = LSD.test(variance,"FH", p.adj="BH")
    MC = MC[["groups"]]
    MC$FH = rownames(MC)
    #MC
    MC = MC[order(MC$FH,decreasing = F),]
    MC = merge(MC,DTStat,by.x = 'FH',by.y = 'Var1')
    if(nrow(MC)<=20){
      p = ggplot(DT,aes(x=FH, y = Phenotype, fill = FH))+
        geom_boxplot(width = 0.2,outlier.colour = NA)+
        scale_x_discrete(breaks = MC$FH[order(MC$Phenotype,decreasing = T)])+
        xlab("")+
        ylab(phe_labs)+
        theme_bw()+
        coord_flip()+
        geom_text(data=MC,size = 6,
                  aes(x=factor(FH),y=quantile*1.03,label=groups))+
        scale_fill_d3(palette = 'category20c',alpha = 0.8)+
        theme(legend.position = 'none',
              plot.margin = unit(c(0.5,1,0,1),'line'),
              text = element_text(family = 'serif',color = 'black'),
              axis.title.y = element_text(family = 'serif',size = 12,color = 'black'),
              axis.title.x = element_text(family = 'serif',size = 12,color = 'black'),
              axis.text.y = element_text(size=12,color = 'black'),
              axis.text.x = element_text(size=10,color = 'black'))
    } else if(nrow(MC)>20 & nrow(MC)<=60){
      p = ggplot(DT,aes(x=FH, y = Phenotype, fill = FH))+
        geom_boxplot(width = 0.2,outlier.colour = NA)+
        scale_x_discrete(breaks = MC$FH[order(MC$Phenotype,decreasing = T)],label = sub('\n',' ', MC$FH[order(MC$Phenotype,decreasing = T)]))+
        xlab("")+
        ylab(phe_labs)+
        theme_bw()+
        coord_flip()+
        geom_text(data=MC,size = 5,
                  aes(x=factor(FH),y=quantile*1.03,label=groups))+
        scale_fill_igv(alpha = 0.8)+
        theme(legend.position = 'none',
              plot.margin = unit(c(0.5,1,0,1),'line'),
              text = element_text(family = 'serif',color = 'black'),
              axis.title.y = element_text(family = 'serif',size = 9,color = 'black'),
              axis.title.x = element_text(family = 'serif',size = 9,color = 'black'),
              axis.text.y = element_text(size=9,color = 'black'),
              axis.text.x = element_text(size=6,color = 'black'))
    } else if (nrow(MC)>60 & nrow(MC)<=110){
      p = ggplot(DT,aes(x=FH, y = Phenotype, fill = FH))+
        geom_boxplot(width = 0.1,outlier.colour = NA)+
        scale_x_discrete(breaks = MC$FH[order(MC$Phenotype,decreasing = T)],label = sub('\n',' ', MC$FH[order(MC$Phenotype,decreasing = T)]))+
        xlab("")+
        ylab(phe_labs)+
        theme_bw()+
        coord_flip()+
        geom_text(data=MC,size = 2,
                  aes(x=factor(FH),y=quantile*1.03,label=groups))+
        scale_fill_igv(alpha = 0.8)+
        theme(legend.position = 'none',
              plot.margin = unit(c(0.5,1,0,1),'line'),
              text = element_text(family = 'serif',color = 'black'),
              axis.title.y = element_text(family = 'serif',size = 6,color = 'black'),
              axis.title.x = element_text(family = 'serif',size = 6,color = 'black'),
              axis.text.y = element_text(size=6,color = 'black'),
              axis.text.x = element_text(size=4,color = 'black'))
    } else {
      p = ggplot(DT,aes(x=FH, y = Phenotype, fill = FH))+
        geom_boxplot(width = 0.1,outlier.colour = NA)+
        scale_x_discrete(breaks = MC$FH[order(MC$Phenotype,decreasing = T)],label = sub('\n',' ', MC$FH[order(MC$Phenotype,decreasing = T)]))+
        xlab("")+
        ylab(phe_labs)+
        theme_bw()+
        coord_flip()+
        geom_text(data=MC,size = 1,
                  aes(x=factor(FH),y=quantile*1.03,label=groups))+
        scale_fill_igv(alpha = 0.8)+
        theme(legend.position = 'none',
              plot.margin = unit(c(0.5,1,0,1),'line'),
              text = element_text(family = 'serif',color = 'black'),
              axis.title.y = element_text(family = 'serif',size = 3,color = 'black'),
              axis.title.x = element_text(family = 'serif',size = 3,color = 'black'),
              axis.text.y = element_text(size=3,color = 'black'),
              axis.text.x = element_text(size=2,color = 'black'))
    }
  }
  p_hap = p
  
  
  snp.gwas.file = paste0("E:/QiGA/document/CottonLab/XLD/Phe/EMMAX/Seprated/Regression/", dataset, ".ps.gz")
  fh.gwas.file = paste0("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Seprate/Phe_", manhattan_names, ".GLMM.csv")
  chr = gff$Chrom[which(gff$Gene==gene)]
  flank = 500000
  gene_start = gff$Start[which(gff$Gene==gene)]
  gene_end = gff$End[which(gff$Gene==gene)]
  gene_middle = floor(gene_start + (gene_end-gene_start)/2)
  
  start = gene_middle-flank
  end = gene_middle+flank
  if(start<0){
    start=0
  }
  
  snp.gwas = fread(snp.gwas.file, data.table = F,header = F,sep='\t',stringsAsFactors = F,encoding = 'UTF-8',nThread = 12)
  fh.gwas = fread(fh.gwas.file, data.table = F, sep=',', header = T, stringsAsFactors = F)
  snp.gwas = na.omit(snp.gwas)
  rawRows = nrow(snp.gwas)
  colnames(snp.gwas) = c("Variants","Beta","SE","P-SNP")
  variantSplit = strsplit(snp.gwas$Variants,'_',fixed=T)
  snp.gwas$CHR = sapply(variantSplit,'[[',1)
  snp.gwas$POS = sapply(variantSplit,'[[',2)
  snp.gwas$POS = as.numeric(snp.gwas$POS)
  snp.gwas = snp.gwas[which(snp.gwas$CHR==chr & snp.gwas$POS >= start & snp.gwas$POS <= end),]
  snp.gwas = snp.gwas[which(!is.na(match(snp.gwas$Variants,Bims$V2))),]
  top = as.numeric(snp.gwas$POS[which(snp.gwas$`P-SNP`==min(snp.gwas$`P-SNP`))])[1]
  # snp.gwas$Ref_Allele = sapply(strsplit(snp.gwas$Variants,'_',fixed=T),'[[',3)
  # snp.gwas$Alt_Allele = sapply(strsplit(snp.gwas$Variants,'_',fixed=T),'[[',4)
  snp.gwas = snp.gwas[,c(5,1,6,4)]
  colnames(snp.gwas) = c("CHR","SNP","BP","P")
  
  target.gwas = snp.gwas
  ld.snps = paste0("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/plot/LD/",chr,"_",round(start/1000000,4),"_",round(end/1000000,4),".snps.txt")
  write.table(target.gwas$SNP,ld.snps,quote=F,col.names = F,row.names = F)
  
  top.snp = target.gwas$SNP[which(target.gwas$P==min(target.gwas$P))][1]
  cmd = paste0("plink2 --r2-phased --silent --aec --chr ", chr, 
               " --from-bp ", start, " --to-bp ", end, 
               " --ld-snp ", top.snp, 
               " --ld-window-kb 9999999 --ld-window-r2 0 --bfile E:/QiGA/document/CottonLab/FH/data/core_Gh_accession.union.MAF0001GENO20 --out E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/plot/LD/"
               ,chr,"_",round(start/1000000,4),"_",round(end/1000000,4))
  system(cmd)
  ld.plink.file = paste0("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/plot/LD/",chr,"_",round(start/1000000,4),"_",round(end/1000000,4),".vcor")
  ld.plink = fread(ld.plink.file,header=T,data.table = F)
  #ld.plink = ld.plink[(grepl(top.snp,ld.plink$ID_A) | grepl(top.snp,ld.plink$ID_B)),c(3,6,7)]
  ld.plink = ld.plink[,c(3,6,7)]
  colnames(ld.plink) = c("SNP_A","SNP_B","R2")
  
  fh.gwas = merge(fh.gwas,gff,'Gene')
  fh.gwas$P = as.numeric(fh.gwas$P)
  fh.gwas$Chisq = as.numeric(fh.gwas$Chisq)
  fh.gwas$Padj = pchisq(fh.gwas$Chisq/mean(fh.gwas$Chisq,na.rm=T),1,0,lower.tail = F)
  gene.target = fh.gwas[fh.gwas$Chrom == chr & fh.gwas$Start<=end & fh.gwas$End>=start,]
  gene.target1 = gene.target[,c(1,15,16,17,18)]
  gene.P.dt = gene.target[,c(1,14)]
  
  target.gwas$BP = target.gwas$BP/1000000
  gene.target1$Start = gene.target1$Start/1000000
  gene.target1$End = gene.target1$End/1000000
  p3 = locus.zoom(data = target.gwas,
                  region = c(chr,start/1000000,end/1000000),
                  offset_bp = 0,nominal = 5,significant = 5,
                  ld.file = ld.plink,
                  genes.data = gene.target1,
                  human=F,
                  genes.pvalue = gene.P.dt,colour.genes = T,rsid.check = F)
  #p3.grob = base2grob(p3)
  
  P = plot_grid(p_manhat,plot_grid(p_hap,p3,ncol = 2,rel_widths = c(1,2)),
                align = 'hv',ncol = 1,rel_heights = c(1,3))

  ggsave(filename = paste0("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/CombinedFigures/FH.", gene,".", outputName, ".png"),
         plot = P,dpi = 300,
         width = 12,height = 9)
  ggsave(filename = paste0("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/CombinedFigures/FH.", gene,".", outputName, ".pdf"),
         plot = P,
         width = 12,height = 9)
}
