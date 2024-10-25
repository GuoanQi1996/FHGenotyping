# calculate effect
# FH hap analysis
library(data.table)
library(scales)
library(stringr)
library(agricolae)
library(Biostrings)
has_common_letter <- function(str1, str2) {
  chars1 <- str_split(str1, "")[[1]]
  if(length(str2)==1){
    chars2 <- str_split(str2, "")[[1]]
    
    common_letters <- intersect(chars1, chars2)
    return(length(common_letters) > 0)
  }
  if(length(str2)>1){
    out = c()
    for(i in 1:length(str2)){
      e = str2[i]
      chars2 <- str_split(e, "")[[1]]
      
      common_letters <- intersect(chars1, chars2)
      out = c(out,length(common_letters) > 0)
    }
    return(out)
  }
}

setwd("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM")

GT = fread("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/FHAP_whole.type.2.txt", sep = '\t', data.table = F, stringsAsFactors = F)
colnames(GT) = sub("0_","",colnames(GT))
GT1 = t(GT[,-1])
colnames(GT1) = GT$Gene
GT1 = data.frame(ID = rownames(GT1),GT1)

plotDT = read.csv("SignificantGene.20240909.csv")
unitDT = read.csv("Units.csv")

name2paper = data.frame(Name = c("DXM1202", "HZG486", "MZY1081", "MZY419"),
                        Paper = c("ICRCAAS_NG_2021.N1202","ZJU_ICP_2022.N486","HEBAU_NG_2021.N1081","HEBAU_NG_2018.N419"))
#CDS = readDNAStringSet("E:/QiGA/document/CottonLab/Cotton_GO/TM-1_V2.1.gene.cds.fa")
higher_better = c("FL","FS","FE","BN","BW","FBN","FFBP","FU","HFFBF","LI","LP","FWPB","MAT","SCI")
lower_better = c("FM","FT","VW","SI","PH")
EffectDT = data.frame()
for(i in 1:nrow(plotDT)){
  #for(i in 1:10){
  cat(paste0("Processing ",i,"/27991.\n"))
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
  
  if(nrow(DTStat)>1){
    variance = aov(Phenotype~FH,data = DT)
    MC = LSD.test(variance,"FH", p.adj="BH")#结果显示：标记字母法out$group
    MC = MC[["groups"]]
    MC$FH = rownames(MC)
    #MC
    MC = MC[order(MC$FH,decreasing = F),]
    MC = merge(MC,DTStat,by.x = 'FH',by.y = 'Var1')
    
    info = MC$FH
    MC$FH = sapply(strsplit(info,"\n(",fixed=T),'[[',1)
    MC$Gene = gene
    MC$PhenotypeName = phe_short
    MC$Source = outputName
    MC = MC[,c("Source","PhenotypeName","Gene","FH","Freq","Phenotype","groups","median","average","quantile")]
    
    MC$EffectSize = 0
    baseIndex = which(MC$FH == 'FH-1')
    if(length(baseIndex)==0){
      baseIndex = which(MC$Freq==max(MC$Freq))
    }
    baseFH = MC$FH[baseIndex]
    basePhe = MC$Phenotype[baseIndex] 
    MC$EffectSize = MC$Phenotype-basePhe
    
    if(phe_short %in% higher_better){
      MC$EffectDirec = ifelse(MC$EffectSize>0,ifelse(!has_common_letter(MC$groups[baseIndex],MC$groups),'Superior','Superior-NS'),ifelse(!has_common_letter(MC$groups[baseIndex],MC$groups),'Inferior','Inferior-NS'))
    } else if (phe_short %in% lower_better){
      MC$EffectDirec = ifelse(MC$EffectSize<0,ifelse(!has_common_letter(MC$groups[baseIndex],MC$groups),'Superior','Superior-NS'),ifelse(!has_common_letter(MC$groups[baseIndex],MC$groups),'Inferior','Inferior-NS'))
    }
    if(baseFH == "FH-1"){
      MC$EffectDirec[baseIndex] = 'Ref'
    } else {
      MC$EffectDirec[baseIndex] = 'Ref-MaxFreq'
    }
    
    EffectDT = rbind(EffectDT,MC)
  }
}
write.table(EffectDT,"E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Effects/Effects.Summary.txt",quote=F,col.names = T,row.names = F,sep = '\t')

# FH effect and phenotype performance
EffectDT = read.table("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Effects/Effects.Summary.txt", header = T, sep = '\t')

length(unique(EffectDT$Gene))
length(unique(EffectDT$Gene[grep("A",EffectDT$Gene)]))
length(unique(EffectDT$Gene[grep("D",EffectDT$Gene)]))
# remove ambiguous
# the gene with opposite effects are the ambiguous one
geneList = unique(EffectDT$Gene)
coreEffectDT = data.frame()
badEffectDT = data.frame()
for(i in 1:length(geneList)){
  cat(paste0("Processing ",i,"/",length(geneList),".\n"))
  gene = geneList[i]
  tempDT = EffectDT[which(EffectDT$Gene==gene),]
  pheList = unique(tempDT$PhenotypeName)
  for(j in 1:length(pheList)){
    phe = pheList[j]
    tempDT2 = tempDT[which(tempDT$PhenotypeName==phe),]
    if(length(unique(tempDT2$Source))<2){
      coreEffectDT = rbind(coreEffectDT,tempDT2)
    } else {
      FHList = unique(tempDT2$FH)
      rmVec = c()
      for(l in 1:length(FHList)){
        FH = FHList[l]
        tempDT3 = tempDT2[which(tempDT2$FH==FH),]
        if("Superior" %in% tempDT3$EffectDirec & "Inferior" %in% tempDT3$EffectDirec){
          rmVec = c(rmVec, which(tempDT2$FH==FH))
          badEffectDT = rbind(badEffectDT,tempDT3)
        }
      }
      if(length(rmVec)>0){
        tempDT2 = tempDT2[-rmVec,]
        coreEffectDT = rbind(coreEffectDT,tempDT2)
      } else {
        coreEffectDT = rbind(coreEffectDT,tempDT2)
      }
    }
  }
}
write.table(coreEffectDT,"./Effects/Effects.Summary.remove_Ambiguous.txt",quote=F,col.names = T,row.names = F,sep = '\t')
write.table(badEffectDT,"./Effects/Effects.Summary.Ambiguous.txt",quote=F,col.names = T,row.names = F,sep = '\t')

# Statistics
library(data.table)
coreEffectDT = read.table("./Effects/Effects.Summary.remove_Ambiguous.txt",header = T,sep = '\t')

length(unique(paste0(coreEffectDT$Gene)[which(coreEffectDT$EffectDirec=='Inferior')]))
length(unique(paste0(coreEffectDT$Gene)[which(coreEffectDT$EffectDirec=='Superior')]))
length(unique(paste0(coreEffectDT$Gene,"_",coreEffectDT$FH)[which(coreEffectDT$EffectDirec=='Inferior')]))
length(unique(paste0(coreEffectDT$Gene,"_",coreEffectDT$FH)[which(coreEffectDT$EffectDirec=='Superior')]))

sigIndex = which(!grepl('NS',coreEffectDT$EffectDirec))
RefIndex = which(grepl('Ref',coreEffectDT$EffectDirec))

sigIndex = setdiff(sigIndex,RefIndex)
sigEffectDT = coreEffectDT[sigIndex,]

SigEffectGene = as.data.frame(table(sigEffectDT$Gene))
SigEffectGene = SigEffectGene[order(SigEffectGene$Freq,decreasing = T),]
colnames(SigEffectGene) = c("Gene","Freq")
write.table(SigEffectGene,"./Effects/SigHapEffect.GeneName.txt",quote=F,col.names = T,row.names = F,sep = '\t')

SigEffectFH = paste0(sigEffectDT$PhenotypeName,".",sigEffectDT$Gene,".",sigEffectDT$FH,".",sigEffectDT$EffectDirec)
SigEffectFHCount = as.data.frame(table(SigEffectFH))
# if(length(duplicated(SigEffectFH))>0){
#   SigEffectFH = SigEffectFH[-which(duplicated(SigEffectFH))]
# }

GT = fread("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/FHAP_whole.type.2.txt", sep = '\t', data.table = F, stringsAsFactors = F)
colnames(GT) = sub("0_","",colnames(GT))
GT1 = t(GT[,-1])
colnames(GT1) = GT$Gene
GT1 = data.frame(ID = rownames(GT1),GT1)

SigEffectFH.MarkVector = paste0(sigEffectDT$Gene,".",sigEffectDT$FH)
summaryDT1 = data.frame()
summaryDT2 = data.frame()
for(i in 1:nrow(GT1)){
  id = rownames(GT1)[i]
  genotype = GT1[i,-1]
  genotype[1,] = paste0(colnames(GT1)[-1],".FH-",genotype[1,])
  matchedSigFH.index = which(!is.na(match(SigEffectFH.MarkVector,genotype[1,])))
  
  indi.effect.dt = sigEffectDT[matchedSigFH.index,]
  write.table(indi.effect.dt,paste0("./Effects/Individual_effect/Sample.",id,".FH.txt"),quote=F,col.names = T,row.names = F,sep = '\t')
  
  matchedSigFH = unique(SigEffectFH[matchedSigFH.index])
  
  matchedSigFH.phe = sapply(strsplit(matchedSigFH,".",fixed=T),'[[',1)
  matchedSigFH.dir = sapply(strsplit(matchedSigFH,".",fixed=T),'[[',4)
  
  all.sig.FH = length(matchedSigFH)
  all.sig.superior = length(which(matchedSigFH.dir=='Superior'))
  all.sig.inferior = length(which(matchedSigFH.dir=='Inferior'))
  record = c(id,all.sig.FH,all.sig.superior,all.sig.inferior)
  summaryDT1 = rbind(summaryDT1,record)
  
  phes = unique(matchedSigFH.phe)
  for(j in 1:length(phes)){
    phe = phes[j]
    temp1 = length(which(matchedSigFH.phe == phe & matchedSigFH.dir == 'Superior'))
    temp2 = length(which(matchedSigFH.phe == phe & matchedSigFH.dir == 'Inferior'))
    record = c(id,phe,temp1,temp2)
    summaryDT2 = rbind(summaryDT2,record)
  }
}
colnames(summaryDT1) = c("ID","ALL_SIG_FH","SUPERIOR_FH","INFERIOR_FH")
colnames(summaryDT2) = c("ID","PHENOTYPE","SUPERIOR_FH","INFERIOR_FH")
summaryDT1$ALL_SIG_FH = as.numeric(summaryDT1$ALL_SIG_FH)
summaryDT1$SUPERIOR_FH = as.numeric(summaryDT1$SUPERIOR_FH)
summaryDT1$INFERIOR_FH = as.numeric(summaryDT1$INFERIOR_FH)

summaryDT2$SUPERIOR_FH = as.numeric(summaryDT2$SUPERIOR_FH)
summaryDT2$INFERIOR_FH = as.numeric(summaryDT2$INFERIOR_FH)
write.table(summaryDT1,"./Effects/Individual.FH.Summary.txt",quote=F,col.names = T,row.names = F, sep = '\t')
write.table(summaryDT2,"./Effects/Individual.FH.Phenotype.txt",quote=F,col.names = T,row.names = F, sep = '\t')

# add annotation
dt1 = read.table("./Effects/Individual.FH.Phenotype.txt",header=T,sep='\t',comment.char = "")
dt2 = read.table("./Effects/Individual.FH.Summary.txt",header=T,sep='\t',comment.char = "")

pop_info = fread("E:/QiGA/document/CottonLab/XLD/population_info/Detailed_Information.txt",header=T,data.table=F)[,c(1,10,19)]
colnames(pop_info)[1] = 'ID'

dt1 = merge(dt1,pop_info,'ID')
dt2 = merge(dt2,pop_info,'ID')

write.table(dt2,"./Effects/Individual.FH.Summary.withInfo.txt",quote=F,col.names = T,row.names = F, sep = '\t')
write.table(dt1,"./Effects/Individual.FH.Phenotype.withInfo.txt",quote=F,col.names = T,row.names = F, sep = '\t')


# correlations with phenotype
library(ggplot2)
library(ggsci)
units = read.csv("Units.csv",header=T)
summaryDT1 = read.table("./Effects/Individual.FH.Summary.txt",header=T,sep = '\t',comment.char = '')
summaryDT2 = read.table("./Effects/Individual.FH.Phenotype.txt",header=T,sep = '\t',comment.char = '')

sampleList = list.files("./Effects/Individual_effect/","*")
pop_info = read.table("E:/QiGA/document/CottonLab/XLD/population_info/Detailed_Information.txt",header=T,sep='\t',quote = "",comment.char = "@",encoding = 'utf_8')
pop_info = pop_info[,c(1,10,11,19)]
colnames(pop_info) = c('ID','Date','Data_Interval','Source')

higher_better = c("FL","FS","FE","BN","BW","FBN","FFBP","FU","HFFBF","LI","LP","FWPB","MAT","SCI")
lower_better = c("FM","FT","VW","SI","PH")

colors = c("#E64B35CC","#3C5488CC",'#C0C0C0CC')
names(colors) = c('positive','negative','Predicted')

pheList = unique(summaryDT2$PHENOTYPE)
for(i in 1:length(pheList)){
  phe = pheList[i]
  pheDT = read.table(paste0("E:/QiGA/document/CottonLab/XLD/Phe/EMMAX_", phe, ".txt"),header=F,sep='\t',comment.char = "")
  colnames(pheDT) = c("FID","ID","VALUE")
  pheDT = pheDT[,c(2,3)]
  
  summary.dt = summaryDT2[which(summaryDT2$PHENOTYPE == phe),]
  summary.dt = merge(summary.dt, pheDT, "ID", all.x = T)
  
  # get averaged effects
  sample_average_effects = data.frame()
  for(j in 1:length(sampleList)){
    sampleFile = paste0("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Effects/Individual_effect/",sampleList[j])
    sample = strsplit(sampleList[j],".",fixed=T)[[1]][2]
    
    temp = read.table(sampleFile,header=T,sep = '\t')
    temp = temp[which(temp$PhenotypeName == phe),]
    if(nrow(temp)==0){
      record = c(sample, phe, 0)
    } else {
      record = c(sample,phe,mean(temp$EffectSize))
    }
    sample_average_effects = rbind(sample_average_effects,record)
  }
  colnames(sample_average_effects) = c("ID","Phenotype","AverageEffects")
  sample_average_effects = sample_average_effects[,c(1,3)]
  
  summary.dt = merge(summary.dt,sample_average_effects,"ID",all.x = T)
  if(phe %in% higher_better){
    summary.dt$Direction = ifelse(summary.dt$AverageEffects>=0,'positive','negative')
  } else {
    summary.dt$Direction = ifelse(summary.dt$AverageEffects>=0,'negative','positive')
  }
  
  summary.dt$Direction[which(is.na(summary.dt$VALUE))] = 'Predicted'
  summary.dt$AverageEffects = as.numeric(summary.dt$AverageEffects)
  summary.dt$AbsAverageEffects = abs(summary.dt$AverageEffects)
  #summary.dt = na.omit(summary.dt)
  lm.Effect.Value = lm(VALUE ~ AverageEffects, data = summary.dt)
  summary.dt$PredictValueBasedOnEffect = predict(lm.Effect.Value,summary.dt)
  summary.dt$DEV = summary.dt$SUPERIOR_FH-summary.dt$INFERIOR_FH
  
  annotations <- data.frame(
    xpos = c(-Inf),
    ypos =  c(Inf),
    annotateText = paste0("r = ",round(cor(summary.dt$DEV,summary.dt$VALUE,use = 'complete.obs'),2)),
    hjustvar = c(-0.2) ,
    vjustvar = c(1.3)) #<- adjust
  annotations2 <- data.frame(
    xpos = c(-Inf),
    ypos =  c(Inf),
    annotateText = paste0("r = ",round(cor(summary.dt$AverageEffects,summary.dt$VALUE,use = 'complete.obs'),2)),
    hjustvar = c(-0.2) ,
    vjustvar = c(1.3)) #<- adjust
  
  summary.dt$VALUE = ifelse(is.na(summary.dt$VALUE),summary.dt$PredictValueBasedOnEffect,summary.dt$VALUE)
  

  p1 = 
    ggplot(summary.dt)+
    geom_point(aes(x = DEV, y = VALUE, color = Direction, size = AbsAverageEffects,fill = Direction))+
    geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText))+
    theme_bw()+
    labs(x = 'Relative numbers of superior FH', y = units$Units[which(units$Phe1 == phe)])+
    scale_color_manual(values = colors)+
    scale_fill_manual(values = colors)+
    theme(text = element_text(colour = 'black',family = 'sans'),
          panel.background = element_rect(fill = 'transparent',colour = 'transparent'),
          plot.background = element_rect(fill = 'transparent',colour = 'transparent'),
          axis.text = element_text(colour = 'black'))
  ggsave(paste0("./Effects/Phe_",phe,".value_vs_DEV.pdf"),p1,width = 6,height = 4.2)
  
  p2 = 
    ggplot(summary.dt)+
    geom_point(aes(x = AverageEffects, y = VALUE, color = Direction, size = AbsAverageEffects,fill = Direction))+
    geom_text(data=annotations2,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText))+
    theme_bw()+
    labs(x = 'Averaged effects per sample', y = units$Units[which(units$Phe1 == phe)])+
    scale_color_manual(values = colors)+
    theme(text = element_text(colour = 'black',family = 'sans'),
          panel.background = element_rect(fill = 'transparent',colour = 'transparent'),
          plot.background = element_rect(fill = 'transparent',colour = 'transparent'),
          axis.text = element_text(colour = 'black'))
  ggsave(paste0("./Effects/Phe_",phe,".value_vs_effect.pdf"),p2,width = 6,height = 4.2)
  
  
  summary.dt1 = merge(summary.dt,pop_info,'ID',all.x = T)
  summary.dt1$Date = as.numeric(summary.dt1$Date)
  
  write.table(summary.dt1,paste0("./Effects/Phenotype_effect/Phe_",phe,".effectPerSample.txt"),quote=F,col.names = T,row.names = F,sep = '\t')
  
  summary.dt1 = summary.dt1[-which(summary.dt1$Data_Interval==''),]
  summary.dt1$Data_Interval = factor(summary.dt1$Data_Interval,levels = c("Before 1900","1900-1949","1950-1959","1960-1969","1970-1979","1980-1989","1990-1999","2000-2009","2010-2019"))
  
  p3 = 
    ggplot(data = summary.dt1)+
    geom_boxplot(aes(x=Data_Interval,y=DEV,color = Data_Interval))+
    theme_bw()+
    labs(y = 'Relative numbers of superior FH')+
    scale_color_npg()+
    theme(legend.position = 'none',
          panel.background = element_rect(fill = 'transparent',colour = 'transparent'),
          plot.background = element_rect(fill = 'transparent',colour = 'transparent'),
          text = element_text(colour = 'black',family = 'sans'),
          axis.title.x = element_blank(),
          axis.text.x = element_text(color = 'black',angle = 45,hjust = 1),
          axis.text.y = element_text(color = 'black'))
  ggsave(paste0("./Effects/Phe_",phe,".Date_vs_DEV.pdf"),p3,width = 8,height = 3)
  
  p4 = 
    ggplot(data = summary.dt1)+
    geom_boxplot(aes(x=Data_Interval,y=AverageEffects,color = Data_Interval))+
    theme_bw()+
    labs(y = 'Averaged effects per sample')+
    scale_color_npg()+
    theme(legend.position = 'none',
          panel.background = element_rect(fill = 'transparent',colour = 'transparent'),
          plot.background = element_rect(fill = 'transparent',colour = 'transparent'),
          text = element_text(colour = 'black',family = 'sans'),
          axis.title.x = element_blank(),
          axis.text.x = element_text(color = 'black',angle = 45,hjust = 1),
          axis.text.y = element_text(color = 'black'))
  ggsave(paste0("./Effects/Phe_",phe,".Date_vs_effect.pdf"),p4,width = 8,height = 3)
}

# plot Date
library(ggplot2)
library(ggsci)
setwd("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM")
summaryDT2 = read.table("./Effects/Individual.FH.Phenotype.txt",header=T,sep = '\t',comment.char = '')

colors = c("#E64B35CC","#3C5488CC",'#C0C0C0CC')
names(colors) = c('positive','negative','Predicted')

pheList = unique(summaryDT2$PHENOTYPE)
for(i in 1:length(pheList)){
  phe = pheList[i]
  
  summary.dt1 = read.table(paste0("./Effects/Phenotype_effect/Phe_",phe,".effectPerSample.txt"),header=T,comment.char = "",sep = '\t')
  summary.dt1 = summary.dt1[-which(summary.dt1$Data_Interval==''),]
  summary.dt1$Data_Interval = factor(summary.dt1$Data_Interval,levels = c("Before 1900","1900-1949","1950-1959","1960-1969","1970-1979","1980-1989","1990-1999","2000-2009","2010-2019"))
  
  p3 = 
    ggplot(data = summary.dt1)+
    geom_boxplot(aes(x=Data_Interval,y=DEV,color = Data_Interval))+
    theme_bw()+
    labs(y = 'Relative numbers of superior FH')+
    scale_color_npg()+
    theme(legend.position = 'none',
          panel.background = element_rect(fill = 'transparent',colour = 'transparent'),
          plot.background = element_rect(fill = 'transparent',colour = 'transparent'),
          text = element_text(colour = 'black',family = 'sans'),
          axis.title.x = element_blank(),
          axis.text.x = element_text(color = 'black',angle = 45,hjust = 1),
          axis.text.y = element_text(color = 'black'))
  ggsave(paste0("./Effects/Phe_",phe,".Date_vs_DEV.pdf"),p3,width = 9,height = 3)
  
  p4 = 
    ggplot(data = summary.dt1)+
    geom_boxplot(aes(x=Data_Interval,y=AverageEffects,color = Data_Interval))+
    theme_bw()+
    labs(y = 'Averaged effects per sample')+
    scale_color_npg()+
    theme(legend.position = 'none',
          panel.background = element_rect(fill = 'transparent',colour = 'transparent'),
          plot.background = element_rect(fill = 'transparent',colour = 'transparent'),
          text = element_text(colour = 'black',family = 'sans'),
          axis.title.x = element_blank(),
          axis.text.x = element_text(color = 'black',angle = 45,hjust = 1),
          axis.text.y = element_text(color = 'black'))
  ggsave(paste0("./Effects/Phe_",phe,".Date_vs_effect.pdf"),p4,width = 9,height = 3)
}

# why FL and LP is negatively correlated
setwd("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM")

dt = read.table("./Effects/Effects.Summary.remove_Ambiguous.txt",header=T,sep = '\t',comment.char = "")

compareLevel = list(c('FL','LP'),c('FS','LP'),c('FL','LI'),c('FS','LI'),c('FL','FS'))

for(l in 1:length(compareLevel)){
  clevels = compareLevel[[l]]
  
  index = which(paste0(dt$PhenotypeName,"-",dt$EffectDirec) %in% c(paste0(clevels[1],"-",c("Superior","Inferior")),paste0(clevels[2],"-",c("Superior","Inferior"))))
  
  tempDT = dt[index,]
  tempDT$Name1 = paste0(tempDT$Gene,"-",tempDT$FH)
  
  geneCount = table(tempDT$Name1)
  targetGene = names(geneCount)[which(geneCount>1)]
  tempDT1 = tempDT[which(!is.na(match(tempDT$Name1,targetGene))),]
  
  summaryDT = data.frame()
  detailedDT = data.frame()
  for(j in 1:length(unique(tempDT1$Name1))){
    FH = unique(tempDT1$Name1)[j]
    gene = strsplit(FH,"-",fixed=T)[[1]][1]
    tempDT2 = tempDT1[which(tempDT1$Name1==FH),]
    
    tempDT2 = tempDT2[order(tempDT2$PhenotypeName,decreasing = F),]
    if(length(unique(tempDT2$PhenotypeName))==1){
      next
    } else {
      record = c(gene,FH,paste0(unique(paste0(tempDT2$PhenotypeName,"-",tempDT2$EffectDirec)),collapse = ','))
      tempDT2$Status = paste0(unique(paste0(tempDT2$PhenotypeName,"-",tempDT2$EffectDirec)),collapse = ',')
      tempDT2$Status2 = ifelse(tempDT2$Status==paste0(clevels[1],'-Inferior,', clevels[2], '-Superior') | tempDT2$Status == paste0(clevels[1],'-Superior,', clevels[2], '-Inferior'),'Opposite','Same')
      detailedDT = rbind(detailedDT,tempDT2)
    }
    summaryDT = rbind(summaryDT,record)
  }
  colnames(summaryDT) = c("Gene","FH","Status")
  x = table(summaryDT$Status)
  print(x)
  pdf(paste0("./Effects/PhenotypeCorrelation/",clevels[1],"_",clevels[2],".pie.pdf"),width = 8,height = 8)
  pie(x,col = rev(pal_npg()(4)),border = NA)
  dev.off()

  colors = c("#008B4599","#EE000099")
  names(colors) = c("Same",'Opposite')
  detailedDT1 = detailedDT[which(detailedDT$PhenotypeName==clevels[1]),]
  p1=
  ggplot(detailedDT1,aes(x = PhenotypeName, y = abs(EffectSize), fill = Status2)) +
    geom_boxplot(width = 0.5)+
    theme_bw()+
    scale_fill_manual(values = colors)+
    labs(y = paste0('Absolute effect size - ',clevels[1]))+
    theme(axis.title.x = element_blank(),
          text = element_text(family = 'sans',color = 'black'),
          axis.text = element_text(family = 'sans',color = 'black'),
          legend.position = 'none')
  ggsave(paste0("./Effects/PhenotypeCorrelation/",clevels[1],"_",clevels[2],".", clevels[1], ".averaged_effect.pdf"),p1,width = 5,height = 5)

  detailedDT2 = detailedDT[which(detailedDT$PhenotypeName==clevels[2]),]
  p2=
    ggplot(detailedDT2,aes(x = PhenotypeName, y = abs(EffectSize), fill = Status2)) +
    geom_boxplot(width = 0.5)+
    theme_bw()+
    scale_fill_manual(values = colors)+
    labs(y = paste0('Absolute effect size - ',clevels[2]))+
    theme(axis.title.x = element_blank(),
          text = element_text(family = 'sans',color = 'black'),
          axis.text = element_text(family = 'sans',color = 'black'),
          legend.position = 'none')
  ggsave(paste0("./Effects/PhenotypeCorrelation/",clevels[1],"_",clevels[2],".", clevels[2], ".averaged_effect.pdf"),p2,width = 5,height = 5)

}

# improvement potential
library(data.table)

dt1 = read.table("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Effects/Individual.FH.Phenotype.txt",header=T,sep='\t',comment.char = '')
dt2 = fread("E:/QiGA/document/CottonLab/XLD/population_info/Detailed_Information.txt",header=T,sep='\t',data.table = F)

dt3 = merge(dt1,dt2,by.x = 'ID',by.y = 'VCF_ID')
dt3$Ratio = dt3$INFERIOR_FH/dt3$SUPERIOR_FH
if(length(which(dt3$Ratio==Inf | dt3$Ratio==0))>0){
  dt3 = dt3[-which(dt3$Ratio==Inf | dt3$Ratio==0),]
}
aggregate(dt3$Ratio,list(dt3$PHENOTYPE,dt3$`Data interval`),mean) -> sumDT

write.table(sumDT,"E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Effects/Improment_potential.superior_slash_inferior.txt",quote=F,col.names = T,row.names = F,sep='\t')
