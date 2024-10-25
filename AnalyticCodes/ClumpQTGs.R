# clump QTGs
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

for(i in 1:length(PHE)){
  #for(i in 1:4){
  pheName = PHE[i]
  fh.glmm.File = paste0("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Seprate/Phe_",pheName,".GLMM.csv")

  fh.glmm = fread(file = fh.glmm.File,data.table = F,sep=',',header = T,stringsAsFactors = F)

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
  fh.glmm.plot = fh.glmm[,c(1,15,16,14)]
  
  colnames(fh.glmm.plot)[4] = 'P'
  fh.glmm.plot$P[which(fh.glmm.plot$P=='HighVaried')] = 1
  fh.glmm.plot$P = as.numeric(fh.glmm.plot$P)
  
  #sigFH = fh.glmm.plot[which(fh.glmm.plot$P<=0.05/nrow(fh.glmm.plot)),]
  sigFH = fh.glmm.plot[which(fh.glmm.plot$P<=1e-5),]
  sigFH = sigFH[order(sigFH$Variants,decreasing = F),]
  clumpDT = data.frame()
  if(nrow(sigFH)>0){
    chrList = unique(sigFH$CHR)
    for(j in 1:length(chrList)){
      chr = chrList[j]
      tempDT = sigFH[which(sigFH$CHR == chr),]
      if(nrow(tempDT)==1){
        clumpDT = rbind(clumpDT,tempDT)
      } else {
        rmIndex = c()
        basePOS = tempDT$START[1]
        for(l in 2:nrow(tempDT)){
          if(tempDT$START[l]-basePOS<500000){
            rmIndex = c(rmIndex,l)
          } else {
            basePOS = tempDT$START[l]
          }
        }
        if(length(rmIndex)>0){
          tempDT = tempDT[-rmIndex,]
        }
        clumpDT = rbind(clumpDT,tempDT)
      }
    }
  }
  write.table(clumpDT,paste0("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Seprate/Clumps/",PHE1[i],".FH.clump.txt"),quote=F,col.names = T,row.names = F,sep = '\t')
}

# summary significant FH and significant SNP
setwd("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Seprate/")

library(data.table)
pheList = list.files("E:/QiGA/document/CottonLab/XLD/Phe/EMMAX/Seprated/Phe/","EMMAX.*.txt")
#PHE = sub("EMMAX_","",sapply(strsplit(pheList,".",fixed=T),'[[',1))
phe_field1 = sapply(strsplit(pheList,".",fixed=T),'[[',2)
phe_field2 = sapply(strsplit(pheList,".",fixed=T),'[[',3)
PHE = paste0(phe_field1,"_",phe_field2)
PHE1 = paste0(phe_field1,".",phe_field2)

FH.STATUS = data.frame()
SNP.STATUS = data.frame()
PropSum = data.frame()
for(flank in seq(500000,5000000,500000)){
  for(i in 1:length(PHE)){
    #for(i in 1:4){
    pheName = PHE[i] 
    fh.glmm.File = paste0("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Seprate/Clumps/", PHE1[i], ".FH.clump.txt")
    snp.glmm.file = paste0("E:/QiGA/document/CottonLab/XLD/Phe/EMMAX/Seprated/Clump/", PHE1[i], ".cl.clumps")
    
    if(file.size(fh.glmm.File)==2 | file.size(snp.glmm.file)==2){
      next
    }
    sigFH = fread(file = fh.glmm.File,data.table = F,sep='\t',header = T,stringsAsFactors = F)
    snp.glmm = fread(file = snp.glmm.file,data.table = F,header = T,sep='\t',stringsAsFactors = F)
    
    sigFH$P = as.numeric(sigFH$P)
    snp.glmm = na.omit(snp.glmm)
    
    sigSNP = snp.glmm[,c(1,2,3)]
    colnames(sigSNP) = c("CHR","POS","Variants")
    
    FH.status = data.frame()
    SNP.status = data.frame()
    for(j in 1:nrow(sigFH)){
      gene = sigFH$Variants[j]
      gene_chr = sigFH$CHR[j]
      gene_pos = sigFH$START[j]
      
      interval.left = gene_pos - flank
      interval.right = gene_pos + flank
      
      targetSNP = sigSNP[which(sigSNP$CHR==gene_chr & sigSNP$POS>=interval.left & sigSNP$POS<=interval.right),]
      cat(nrow(targetSNP))
      if(nrow(targetSNP)==0){
        record.FH = c(flank,gene,pheName,'Unique-Sig-FH','NA')
      } else {
        record.FH = c(flank,gene,pheName,'Shared-Sig-FH',mean(abs(gene_pos-targetSNP$POS)))
        record.SNP = data.frame(Flank = flank,SNP = targetSNP$Variants,Source = pheName,Status = 'Shared-Sig-SNP',Distance = targetSNP$POS-gene_pos)
        SNP.status = rbind(SNP.status,record.SNP)
      }
      FH.status = rbind(FH.status,record.FH)
    }
    colnames(FH.status) = c('Flank','Gene','Source','Status','Distance')
    if(length(which(duplicated(FH.status$Gene)))>0){
      FH.status = FH.status[-which(duplicated(FH.status$Gene)),]
    }
    if(nrow(SNP.status)>0){
      colnames(SNP.status) = c('Flank','SNP','Source','Status','Distance')
    }
    if(length(which(duplicated(SNP.status$SNP)))>0){
      SNP.status = SNP.status[-which(duplicated(SNP.status$SNP)),]
    }
    SNP.unique = data.frame(Flank = flank, SNP = setdiff(sigSNP$Variants,SNP.status$SNP),Source = pheName, Status = "Unique-Sig-SNP", Distance = 'NA')
    SNP.status = rbind(SNP.status,SNP.unique)
    
    FH.STATUS = rbind(FH.STATUS,FH.status)
    SNP.STATUS = rbind(SNP.STATUS,SNP.status)

    uniqueFH = length(which(FH.status$Status=='Unique-Sig-FH'))
    sharedFH = length(which(FH.status$Status!='Unique-Sig-FH'))
    uniqueSNP = length(which(SNP.status$Status!='Shared-Sig-SNP'))
    sharedSNP = length(which(SNP.status$Status=='Shared-Sig-SNP'))
    
    shared_in_SigSNP_Propor = sharedSNP/(sharedSNP+uniqueSNP)
    shared_in_SigFH_Propor = sharedFH/(sharedFH+uniqueFH)
    record = c(flank,pheName,uniqueSNP,sharedSNP,uniqueFH,sharedFH,shared_in_SigSNP_Propor,shared_in_SigFH_Propor)
    PropSum = rbind(PropSum,record)
  }
  cat('\n')
}
SNP.STATUS$Distance = as.numeric(SNP.STATUS$Distance)
FH.STATUS$Distance = as.numeric(FH.STATUS$Distance)
SNP.STATUS$Flank = as.numeric(SNP.STATUS$Flank)
FH.STATUS$Flank = as.numeric(FH.STATUS$Flank)
colnames(PropSum) = c('Flank','Source','uniqueSNP','sharedSNP','uniqueFH','sharedFH','OverlappedSNPProp','OverlappedFHProp')
PropSum$OverlappedSNPProp = as.numeric(PropSum$OverlappedSNPProp)
PropSum$OverlappedFHProp = as.numeric(PropSum$OverlappedFHProp)
PropSum$Flank = as.numeric(PropSum$Flank)

plot(aggregate(PropSum$OverlappedSNPProp,list(PropSum$Flank),mean)[,'x'])
plot(aggregate(PropSum$OverlappedFHProp,list(PropSum$Flank),mean)[,'x'])
PropSum.SNPMean = aggregate(PropSum$OverlappedSNPProp,list(PropSum$Flank),mean)
PropSum.FHMean = aggregate(PropSum$OverlappedFHProp,list(PropSum$Flank),mean)
PropSum.SNPMean$IncreaseRate = "-"
PropSum.FHMean$IncreaseRate = "-"
for(i in 2:nrow(PropSum.FHMean)){
  x1 = (PropSum.SNPMean$x[i]-PropSum.SNPMean$x[i-1])/PropSum.SNPMean$x[i-1]
  x2 = (PropSum.FHMean$x[i]-PropSum.FHMean$x[i-1])/PropSum.FHMean$x[i-1]
  PropSum.SNPMean$IncreaseRate[i] = x1
  PropSum.FHMean$IncreaseRate[i] = x2
}
colnames(PropSum.FHMean) = colnames(PropSum.SNPMean) = c("Flank","Representative","IncreaseRate")

write.table(SNP.STATUS,"E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Overlap_with_SNPs.txt",quote=F,col.names = T,row.names = F,sep = '\t')
write.table(FH.STATUS,"E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Overlap_with_FHs.txt",quote=F,col.names = T,row.names = F,sep = '\t')
write.table(PropSum,"E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Overlap_proportion.txt",quote=F,col.names = T,row.names = F,sep = '\t')
write.table(PropSum.SNPMean,"E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Overlap_proportion_SNP_increase.txt",quote=F,col.names = T,row.names = F,sep = '\t')
write.table(PropSum.FHMean,"E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Overlap_proportion_FH_increase.txt",quote=F,col.names = T,row.names = F,sep = '\t')

# plot dot
library(ggplot2)

PropSum.SNPMean = read.table("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Overlap_proportion_SNP_increase.txt",header=T)
PropSum.FHMean = read.table("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Overlap_proportion_FH_increase.txt",header=T)

PropSum.FHMean$Type = 'QTG'
PropSum.SNPMean$Type = 'QTL'

DT = rbind(PropSum.SNPMean,PropSum.FHMean)
DT$Type = factor(DT$Type,levels = c('QTL','QTG'))
DT = DT[which(!is.na(match(DT$Flank,c(1000000,2000000,3000000,4000000,5000000)))),]
DT$Flank2 = paste0(DT$Flank/1000000)
DT$Representative = DT$Representative*100
cols = c("#0073C2","#EFC000")
names(cols) = c('QTL','QTG')

p = 
ggplot(data = DT, aes(x = Flank2, y = Representative, color = Type, group = Type))+
  geom_point(aes(shape = Type),size = 4)+
  geom_line()+
  scale_color_manual(values = cols)+
  labs(x = 'Flanking region (Mb)', y = 'Percentage of overlapping (%)')+
  theme_bw()+
  theme(text = element_text(family = 'sans',colour = 'black'),
        legend.background = element_rect(fill = 'transparent',colour = 'transparent'),
        legend.position = c(0.9,0.1),
        legend.title = element_blank(),
        plot.background = element_rect(fill = 'transparent',colour = 'transparent'),
        panel.background = element_rect(fill = 'transparent',colour = 'transparent'),
        axis.title = element_text(size = 18),axis.text = element_text(size = 16,color='black'),legend.text = element_text(size = 16))
ggsave("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/OverlappPercentage.pdf",p,width = 7,height = 4)

# Venn plot
SNP.Status = read.table("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Overlap_with_SNPs.txt",header=T,sep='\t')
FH.Status = read.table("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Overlap_with_FHs.txt",header=T,sep='\t')

colnames(SNP.Status)[2] = colnames(FH.Status)[2] = 'Variants'

flanks = unique(SNP.Status$Flank)
cols = c("#0073C2","#EFC000")
for(i in 1:length(flanks)){
  flank = flanks[i]
  
  SNPSet = SNP.Status[which(SNP.Status$Flank==flank),]
  FHSet = FH.Status[which(FH.Status$Flank==flank),]

  snp.prop = as.data.frame(table(SNPSet$Status)/nrow(SNPSet))
  fh.prop = as.data.frame(table(FHSet$Status)/nrow(FHSet))
  
  p1=
  ggplot(snp.prop, aes(x = 0, y = Freq, fill = Var1)) + 
    geom_bar(stat = "identity", width = 1) +    
    coord_polar(theta = "y",start = 0) + 
    theme(axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          plot.margin=unit(rep(0,4),'lines'),
          legend.title = element_blank(),
          legend.position = "none",
          panel.grid  = element_blank(),
          panel.border=element_blank(),
          strip.background = element_rect(fill='transparent'),
          panel.background = element_rect(fill='transparent'),
          plot.background = element_rect(fill='transparent',colour = NA),
          legend.background = element_rect(fill='transparent'))+
    scale_fill_manual(values = rev(cols),name="")
  ggsave(plot = p1,filename = paste0("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Clump/SNP.flank_",flank,".pdf"),width = 1.5,height = 1.5,bg = "transparent")
  
  p2=
    ggplot(fh.prop, aes(x = 0, y = Freq, fill = Var1)) + 
    geom_bar(stat = "identity", width = 1) +    
    coord_polar(theta = "y",start = 0) + 
    theme(axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          plot.margin=unit(rep(0,4),'lines'),
          legend.title = element_blank(),
          legend.position = "none",
          panel.grid  = element_blank(),
          panel.border=element_blank(),
          strip.background = element_rect(fill='transparent'),
          panel.background = element_rect(fill='transparent'),
          plot.background = element_rect(fill='transparent',colour = NA),
          legend.background = element_rect(fill='transparent'))+
    scale_fill_manual(values = cols,name="")
  ggsave(plot = p2,filename = paste0("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Clump/FH.flank_",flank,".pdf"),width = 1.5,height = 1.5,bg = "transparent")
  
}

# plot bar
# 1. snp VS clumped snp
# snps
typeDT = data.frame(pheShort = c("FE","FL","FM","FS","BN","BW","FBN","FFBP","FT","FU","HFFBF","LI","LP","PH","SI","VW","FWPB","MAT","SCI","Wilt","Salt"),
                    Type = c(rep("Fiber quality",4),rep("Lint yield",3),"Plant architecture","Maturity","Fiber quality","Plant architecture","Lint yield","Lint yield","Plant architecture","Lint yield","Resistence","Lint yield","Maturity","Fiber quality","Resistence","Resistence"))

snp.glmm.files = list.files("E:/QiGA/document/CottonLab/XLD/Phe/EMMAX/Seprated/Regression/","*.*ps$",full.names = T)
raw.snps = data.frame()
raw.snps.nums = c()
for(i in 1:length(snp.glmm.files)){
  file = snp.glmm.files[i]
  cat(paste0("Processing ",i,"/",length(snp.glmm.files)," ...\n"))
  source = strsplit(strsplit(file,"/",fixed=T)[[1]][10],".",fixed=T)[[1]][1]
  pheName = strsplit(strsplit(file,"/",fixed=T)[[1]][10],".",fixed=T)[[1]][2]
  if(file.size(file)==2){
    next
  }
  dt = fread(file,header=T,sep='\t',data.table=F)
  dt$P = as.numeric(dt$P)
  dt = na.omit(dt)
  
  sigDT = dt[which(dt$P<=1e-5),]
  raw.snps.nums = c(raw.snps.nums,nrow(sigDT))
  if(nrow(sigDT)==0){
    next
  }
  sigDT$phe = pheName
  sigDT$Source = source
  raw.snps = rbind(raw.snps,sigDT)
}
raw.snps$pheShort = sapply(strsplit(raw.snps$phe,'_',fixed=T),'[[',1)
raw.snps2 = merge(raw.snps,typeDT,by = 'pheShort',sort=F)
write.table(raw.snps2,"E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Summary/Raw.Sig.SNP.txt",quote=F,col.names = T,row.names = F,sep = '\t')

unique_index = which(!duplicated(raw.snps2$ID))
raw.snps2$CHR = sapply(strsplit(raw.snps2$ID,'_',fixed=T),'[[',1)
table(raw.snps2$CHR[unique_index])

raw.snps2$Field1 = paste0(raw.snps2$ID,"_",raw.snps2$Type)
if(length(duplicated(raw.snps2$Field1))>0){
  raw.snps3 = raw.snps2[-which(duplicated(raw.snps2$Field1)),]
}
table(raw.snps3$CHR,raw.snps3$Type)

# clumped snps
typeDT = data.frame(pheShort = c("FE","FL","FM","FS","BN","BW","FBN","FFBP","FT","FU","HFFBF","LI","LP","PH","SI","VW","FWPB","MAT","SCI","Wilt","Salt"),
                    Type = c(rep("Fiber quality",4),rep("Lint yield",3),"Plant architecture","Maturity","Fiber quality","Plant architecture","Lint yield","Lint yield","Plant architecture","Lint yield","Resistence","Lint yield","Maturity","Fiber quality","Resistence","Resistence"))

snp.glmm.files = list.files("E:/QiGA/document/CottonLab/XLD/Phe/EMMAX/Seprated/Clump/","*.*clumps$",full.names = T)
clump.snps = data.frame()
clump.snps.nums = c()
for(i in 1:length(snp.glmm.files)){
  file = snp.glmm.files[i]
  cat(paste0("Processing ",i,"/",length(snp.glmm.files)," ...\n"))
  source = strsplit(strsplit(file,"/",fixed=T)[[1]][10],".",fixed=T)[[1]][1]
  pheName = strsplit(strsplit(file,"/",fixed=T)[[1]][10],".",fixed=T)[[1]][2]
  if(file.size(file)==2){
    next
  }
  dt = fread(file,header=T,sep='\t',data.table=F)
  dt = na.omit(dt)
  
  sigDT = dt
  clump.snps.nums = c(clump.snps.nums,nrow(sigDT))
  if(nrow(sigDT)==0){
    next
  }
  sigDT$phe = pheName
  sigDT$Source = source
  clump.snps = rbind(clump.snps,sigDT)
}
clump.snps$pheShort = sapply(strsplit(clump.snps$phe,'_',fixed=T),'[[',1)
clump.snps2 = merge(clump.snps,typeDT,by = 'pheShort',sort=F)
length(unique(clump.snps2$ID))
write.table(clump.snps2,"E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Summary/Clump.Sig.SNP.txt",quote=F,col.names = T,row.names = F,sep = '\t')

unique_index = which(!duplicated(clump.snps2$ID))
clump.snps2$CHR = sapply(strsplit(clump.snps2$ID,'_',fixed=T),'[[',1)
table(clump.snps2$CHR[unique_index])

clump.snps2$Field1 = paste0(clump.snps2$ID,"_",clump.snps2$Type)
if(length(duplicated(clump.snps2$Field1))>0){
  clump.snps3 = clump.snps2[-which(duplicated(clump.snps2$Field1)),]
}
table(clump.snps3$CHR,clump.snps3$Type)

# 2. FH VS clumped FH
# FH
typeDT = data.frame(pheShort = c("FE","FL","FM","FS","BN","BW","FBN","FFBP","FT","FU","HFFBF","LI","LP","PH","SI","VW","FWPB","MAT","SCI","Wilt","Salt"),
                    Type = c(rep("Fiber quality",4),rep("Lint yield",3),"Plant architecture","Maturity","Fiber quality","Plant architecture","Lint yield","Lint yield","Plant architecture","Lint yield","Resistence","Lint yield","Maturity","Fiber quality","Resistence","Resistence"))
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

qtgs = data.frame()
for(i in 1:length(PHE)){
  #for(i in 1:4){
  cat(paste0("Processing ",i,"/",length(PHE)," ...\n"))
  pheName = PHE[i]
  fh.glmm.File = paste0("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Seprate/Phe_",pheName,".GLMM.csv")
  
  fh.glmm = fread(file = fh.glmm.File,data.table = F,sep=',',header = T,stringsAsFactors = F)
  
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
  
  sigDT = fh.glmm.plot[which(fh.glmm.plot$P<=1e-5),]
  if(nrow(sigDT)==0){
    next
  }
  sigDT$phe = phe_field2[i]
  sigDT$Source = phe_field1[i]
  qtgs = rbind(qtgs,sigDT)
  
}

qtgs$pheShort = sapply(strsplit(qtgs$phe,'_',fixed=T),'[[',1)
qtgs2 = merge(qtgs,typeDT,by = 'pheShort',sort=F)
write.table(qtgs2,"E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Summary/Sig.FH.txt",quote=F,col.names = T,row.names = F,sep = '\t')

unique_index = which(!duplicated(qtgs2$Variants))
qtgs2$CHR = substr(qtgs2$Variants,4,6)
table(qtgs2$CHR[unique_index])

qtgs2$Field1 = paste0(qtgs2$Variants,"_",qtgs2$Type)
if(length(duplicated(qtgs2$Field1))>0){
  qtgs3 = qtgs2[-which(duplicated(qtgs2$Field1)),]
}
table(qtgs3$CHR,qtgs3$Type)

pleioQTGs = as.data.frame(table(qtgs3$Variants)[which(table(qtgs3$Variants)>=2)])

# Clumped FH
typeDT = data.frame(pheShort = c("FE","FL","FM","FS","BN","BW","FBN","FFBP","FT","FU","HFFBF","LI","LP","PH","SI","VW","FWPB","MAT","SCI","Wilt","Salt"),
                    Type = c(rep("Fiber quality",4),rep("Lint yield",3),"Plant architecture","Maturity","Fiber quality","Plant architecture","Lint yield","Lint yield","Plant architecture","Lint yield","Resistence","Lint yield","Maturity","Fiber quality","Resistence","Resistence"))

qtgs.clumped = data.frame()
clumpFiles = list.files("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Seprate/Clumps/","*.txt",full.names = T)
for(i in 1:length(clumpFiles)){
  file = clumpFiles[i]
  cat(paste0("Processing ",i,"/",length(clumpFiles)," ...\n"))
  source = strsplit(strsplit(file,"/",fixed=T)[[1]][10],".",fixed=T)[[1]][1]
  pheName = strsplit(strsplit(file,"/",fixed=T)[[1]][10],".",fixed=T)[[1]][2]
  if(file.size(file)==2){
    next
  }
  dt = fread(file,header=T,sep='\t',data.table=F)
  dt = na.omit(dt)
  
  sigDT = dt
  if(nrow(sigDT)==0){
    next
  }
  sigDT$phe = pheName
  sigDT$Source = source
  qtgs.clumped = rbind(qtgs.clumped,sigDT)
}
qtgs.clumped$pheShort = sapply(strsplit(qtgs.clumped$phe,'_',fixed=T),'[[',1)
qtgs.clumped2 = merge(qtgs.clumped,typeDT,by = 'pheShort',sort=F)
write.table(qtgs.clumped2,"E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Summary/Clump.Sig.FH.txt",quote=F,col.names = T,row.names = F,sep = '\t')

unique_index = which(!duplicated(qtgs.clumped2$Variants))
qtgs.clumped2$CHR = substr(qtgs.clumped2$Variants,4,6)
table(qtgs.clumped2$CHR[unique_index])

qtgs.clumped2$Field1 = paste0(qtgs.clumped2$Variants,"_",qtgs.clumped2$Type)
if(length(duplicated(qtgs.clumped2$Field1))>0){
  qtgs.clumped3 = qtgs.clumped2[-which(duplicated(qtgs.clumped2$Field1)),]
}
table(qtgs.clumped3$CHR,qtgs.clumped3$Type)

# plot clumped FH and snp
library(ggplot2)
library(ggsci)
library(cowplot)

dt = read.table("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Summary/Clumped_FH_SNP_Summary.txt",header=T,sep = '\t')
PHEs = unique(dt$Phenotype)
cols = c("#0073C2","#EFC000")
names(cols) = c("SNP/Indel","FH")
for(i in 1:length(PHEs)){
  phe = PHEs[i]
  plotDT = dt[which(dt$Phenotype==phe),]
  
  cat(cor(plotDT$Significant[which(plotDT$Type=='SNP/Indel')],plotDT$Significant[which(plotDT$Type=='FH')]))
  cat('\n')
  
  # assign(paste0("p",i),
  #        ggplot(plotDT,aes(x = Chromosome, y = Significant, fill = Type))+
  #          geom_bar(stat = 'identity', position = position_dodge())+
  #          scale_fill_manual(values = cols)+
  #          labs(x = NULL,y = 'Number of\nsignificant variants')+
  #          theme(legend.position = 'none',
  #                axis.line = element_line(colour = 'black',linewidth = 0.5),
  #                text = element_text(family = 'sans',color = 'black'),
  #                panel.grid = element_blank(),
  #                panel.background = element_rect(fill = 'transparent',colour = 'transparent')))
  # ggsave(filename = paste0("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Summary/Compare_FH_SNP.",phe,".pdf"),
  #        plot = get(paste0("p",i)),
  #        width = 8,height = 2)
  # 
  # 
  # assign(paste0("pe",i),
  #        ggplot(plotDT,aes(x = Chromosome, y = Significant, fill = Type))+
  #          geom_bar(stat = 'identity', position = position_dodge())+
  #          scale_fill_manual(values = cols)+
  #          labs(x = NULL,y = 'Number of\nsignificant variants')+
  #          theme(legend.position = 'none',
  #                axis.line = element_line(colour = 'black',linewidth = 0.5),
  #                text = element_blank(),
  #                panel.grid = element_blank(),
  #                panel.background = element_rect(fill = 'transparent',colour = 'transparent')))
  # ggsave(filename = paste0("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Summary/Compare_FH_SNP.",phe,".emptyAxis.pdf"),
  #        plot = get(paste0("pe",i)),
  #        width = 10,height = 2)
  
}
plot.list = list(p1,p2,p3,p4,p5)
P = plot_grid(plotlist = plot.list,align = 'hv',nrow = 5)
ggsave(filename = paste0("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Summary/Compare_FH_SNP.Combined.pdf"),
       plot = P,
       width = 8,height = 12)

# plot Subgenome FH-SNP
library(ggplot2)
library(ggsci)
library(cowplot)

dt = read.table("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Summary/Clump_FH_SNP_Summary_SubGenome.txt",header=T,sep = '\t')
PHEs = unique(dt$Phenotype)
TYPEs = unique(dt$Type)
cols = c("#EE0000CC","#631879CC")
for(i in 1:length(PHEs)){
  phe = PHEs[i]
  for(j in 1:length(TYPEs)){
    type = TYPEs[j]
    plotDT = dt[which(dt$Phenotype==phe & dt$Type==type),]
    
    plotDT$Prop = plotDT$Significant/sum(plotDT$Significant)
    print(plotDT)
    p1 = 
      ggplot(plotDT, aes(x = 0, y = Prop, fill = SubGenome)) + 
      geom_bar(stat = "identity", width = 1) +    
      coord_polar(theta = "y",start = 0) + 
      theme(axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.line = element_blank(),
            plot.margin=unit(rep(0,4),'lines'),
            legend.title = element_blank(),
            legend.position = "none",
            panel.grid  = element_blank(),
            panel.border=element_blank(),
            strip.background = element_rect(fill='transparent'),
            panel.background = element_rect(fill='transparent'),
            plot.background = element_rect(fill='transparent',colour = NA),
            legend.background = element_rect(fill='transparent'))+
      scale_fill_manual(values = cols,name="")
    if(type == 'SNP/Indel'){
      type = 'SNP'
    }
    phe2 = sub(" ","_",phe)
    ggsave(plot = p1,filename = paste0("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Summary/SubGenome.",phe2,".",type,".pie.pdf"),width = 1.5,height = 1.5,bg = "transparent")
    
  }
}
