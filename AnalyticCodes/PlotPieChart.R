library(data.table)
library(ggsci)
sample.info = fread("E:/QiGA/document/CottonLab/XLD/population_info/Detailed_Information.txt",header = T,sep = '\t')[,c(1,11,14)]
colnames(sample.info) = c("Sample","Date","Zone")

sample.info$Date[which(sample.info$Date=="1950-1959" | sample.info$Date=="1960-1969")] = "1950-1969"
sample.info$Date[which(sample.info$Date=="1970-1979" | sample.info$Date=="1980-1989")] = "1970-1989"
sample.info$Date[which(sample.info$Date=="1990-1999" | sample.info$Date=="2000-2009")] = "1990-2009"
cols = c(`SWC`="#990080CC",`YZR`="#0099CCCC",`YER`="#BA6338CC",
         `NWC`="#339900CC",`SCL`="#8A4198CC",`NC`="#003399CC",
         USA="#020249CC",
         `AL`="#8C564BCC",
         GB="#EFC000CC",
         `EC`="#5DB1DDCC",
         FSU="#DC0000CC",Africa="#FF7F0ECC",Uncertain="#7F7F7FCC",
         Europe="#8F7700CC",
         Oceania="#17BECFCC",
         Asia="#008B45CC",
         CHN="#E64B35CC",
         Others="#7F7F7FCC")
pal_continuous = c("#7F7F7FCC",pal_gsea(palette = c("default"), n = length(unique(sample.info$Date)), alpha = 0.6, reverse = FALSE)(length(unique(sample.info$Date))))
#pal_continuous = pal_material("blue",n = length(unique(date_info$Date)),alpha = 0.8)(length(unique(date_info$Date)))
names(pal_continuous) = c("Before 1900","1900-1949","1950-1969","1970-1989","1990-2009","2010-2019")
for(i in 1:6){
  groupData = read.table(paste0("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/Phy/Group/G", i, ".FH.sample.list.txt"),header = F)
  
  dt = data.frame(Sample = groupData$V1)
  dt = merge(sample.info,dt,sort=F,all.y=T)
  dt$Date[which(is.na(dt$Date))] = "Uncertain"
  dt$Zone[which(is.na(dt$Zone))] = "Uncertain"
  
  sourceConsistution = table(dt$Zone)/sum(table(dt$Zone))
  dt$Zone[which(!is.na(match(dt$Zone,names(sourceConsistution)[which(sourceConsistution<0.05)])))]="Others"
  sourceConsistution = table(dt$Zone)/sum(table(dt$Zone))
  sourceConsistution = data.frame(Source = names(sourceConsistution),Proportion = as.numeric(sourceConsistution))
  sourceConsistution$Source[which(sourceConsistution$Source=='Uncertain')] = "Others"
  sourceConsistution = aggregate(sourceConsistution$Proportion,by = list(sourceConsistution$Source),sum)
  if(length(unique(sourceConsistution$Source))==1){
    sourceConsistution = data.frame(Source="Others",Proportion=1)
  }
  print(sourceConsistution)
  p1 = 
    ggplot(sourceConsistution, aes(x = 0, y = x, fill = Group.1)) + 
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
  ggsave(plot = p1,filename = paste0("./PIE/G",i,".Source_pie.pdf"),width = 1.5,height = 1.5,bg = "transparent")
  
  dateConsistution = table(dt$Date)/sum(table(dt$Date))
  dt$Date[which(!is.na(match(dt$Date,names(dateConsistution)[which(dateConsistution<0.01)])))]="Others"
  dateConsistution = table(dt$Date)/sum(table(dt$Date))
  dateConsistution = data.frame(Date = names(dateConsistution),Proportion = as.numeric(dateConsistution))
  dateConsistution$Date[which(dateConsistution$Date=='Uncertain')] = "Others"
  dateConsistution = aggregate(dateConsistution$Proportion,by = list(dateConsistution$Date),sum)
  if(length(unique(dateConsistution$Date))==1){
    dateConsistution = data.frame(Date="Others",Proportion=1)
  }
  
  p1 = 
    ggplot(dateConsistution, aes(x = 0, y = x, fill = Group.1)) + 
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
    scale_fill_manual(values = pal_continuous,name="")
  ggsave(plot = p1,filename = paste0("./PIE/G",i,".Date_pie.pdf"),width = 1.5,height = 1.5, bg = "transparent")
}

# SNP pie
# plot pie
library(ggplot2)
setwd("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/Phy/PIE/")
cols = c(`SWC`="#990080CC",`YZR`="#0099CCCC",`YER`="#BA6338CC",
         `NWC`="#339900CC",`SC`="#8A4198CC",`NC`="#003399CC",
         USA="#020249CC",
         `AL`="#8C564BCC",
         GB="#EFC000CC",
         `EC`="#5DB1DDCC",
         FSU="#DC0000CC",Africa="#FF7F0ECC",Uncertain="#7F7F7FCC",
         `Central America`="#D43F3ACC",`South America`="#E18727CC",`North America`="#EEA236CC",
         Europe="#8F7700CC",
         Oceania="#17BECFCC",
         Asia="#008B45CC",
         CHN="#E64B35CC",
         Others="#7F7F7FCC")
dt = read.table("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/Phy/PIE/A08_FH.txt",header=T,sep = '\t',comment.char = "")
dt$Source = sub(" (CHN)","",dt$Source,fixed=T)


hapType = unique(dt$Haplotype)
for(i in 1:length(hapType)){
  hap = hapType[i]
  temp = dt[which(dt$Haplotype == hap),]
  sourceConsistution = table(temp$Source)/sum(table(temp$Source))
  temp$Source[which(!is.na(match(temp$Source,names(sourceConsistution)[which(sourceConsistution<0.05)])))]="Others"
  sourceConsistution = table(temp$Source)/sum(table(temp$Source))
  sourceConsistution = data.frame(Source = names(sourceConsistution),Proportion = as.numeric(sourceConsistution))
  sourceConsistution$Source[which(sourceConsistution$Source=='Uncertain')] = "Others"
  sourceConsistution = aggregate(sourceConsistution$Proportion,by = list(sourceConsistution$Source),sum)
  if(length(unique(sourceConsistution$Source))==1){
    sourceConsistution = data.frame(Source="Others",Proportion=1)
  }
  print(sourceConsistution)
  p1 =
    ggplot(sourceConsistution, aes(x = 0, y = x, fill = Group.1)) +
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
  ggsave(plot = p1,filename = paste0("./A08_FH_Hap",i,".Source_pie.pdf"),width = 1.5,height = 1.5,bg = "transparent")
}

# FH-block pie
library(ggplot2)
setwd("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/Phy/PIE/")
cols = c(`SWC`="#990080CC",`YZR`="#0099CCCC",`YER`="#BA6338CC",
         `NWC`="#339900CC",`SC`="#8A4198CC",`NC`="#003399CC",
         USA="#020249CC",
         `AL`="#8C564BCC",
         GB="#EFC000CC",
         `EC`="#5DB1DDCC",
         FSU="#DC0000CC",Africa="#FF7F0ECC",Uncertain="#7F7F7FCC",
         `Central America`="#D43F3ACC",`South America`="#E18727CC",`North America`="#EEA236CC",
         Europe="#8F7700CC",
         Oceania="#17BECFCC",
         Asia="#008B45CC",
         CHN="#E64B35CC",
         Others="#7F7F7FCC")
dt = read.table("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/Phy/PIE/A08_FH_Block.txt",header=T,sep = '\t',comment.char = "")
dt$Source = sub(" (CHN)","",dt$Source,fixed=T)


hapType = unique(dt$Haplotyp)
for(i in 1:length(hapType)){
  hap = hapType[i]
  temp = dt[which(dt$Haplotype == hap),]
  sourceConsistution = table(temp$Source)/sum(table(temp$Source))
  temp$Source[which(!is.na(match(temp$Source,names(sourceConsistution)[which(sourceConsistution<0.05)])))]="Others"
  sourceConsistution = table(temp$Source)/sum(table(temp$Source))
  sourceConsistution = data.frame(Source = names(sourceConsistution),Proportion = as.numeric(sourceConsistution))
  sourceConsistution$Source[which(sourceConsistution$Source=='Uncertain')] = "Others"
  sourceConsistution = aggregate(sourceConsistution$Proportion,by = list(sourceConsistution$Source),sum)
  if(length(unique(sourceConsistution$Source))==1){
    sourceConsistution = data.frame(Source="Others",Proportion=1)
  }
  print(sourceConsistution)
  p1 =
    ggplot(sourceConsistution, aes(x = 0, y = x, fill = Group.1)) +
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
  ggsave(plot = p1,filename = paste0("./A08_Block_Hap_",hap,".Source_pie.pdf"),width = 1.5,height = 1.5,bg = "transparent")
}