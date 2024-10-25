# Sankey plot
library(ggalluvial)
library(dplyr)
library(ggplot2)

snp_DT = data.frame()
FH_DT = data.frame()
for(i in 1:6){
  snp_group = paste0("E:/QiGA/document/CottonLab/XLD/Imputation/Phy/ALL/group/G", i, ".New.sample.list")
  FH_group = paste0("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/Phy/Group/G", i, ".FH.sample.list.txt")
  
  snp_dt = read.table(snp_group,header=F)
  snp_dt = data.frame(Group = paste0("G",i),Sample = snp_dt$V2)
  snp_DT = rbind(snp_DT, snp_dt)
  
  FH_dt = read.table(FH_group,header=F)
  FH_dt = data.frame(Group = paste0("G",i),Sample = FH_dt$V1)
  FH_DT = rbind(FH_DT, FH_dt)
}
DT = merge(snp_DT,FH_DT,'Sample')
DT$Direction = paste0(DT$Group.x,"_",DT$Group.y)
plotDT = as.data.frame(table(DT$Direction))

plotDT = plotDT[-which(plotDT$Freq<30),]
plotDT$Var1 = as.character(plotDT$Var1)
plotDT$FromSNP = sapply(strsplit(plotDT$Var1,"_",fixed=T),'[[',1)
plotDT$ToFH = sapply(strsplit(plotDT$Var1,"_",fixed=T),'[[',2)
cols = c(G1="#EFC000CC",G2="#0099CCCC",G3="#BA6338CC",G4="#631879CC",G5="#339900CC",G6="#D43F3ACC")
p = 
  ggplot(plotDT,aes(y = Freq, axis1 = ToFH, axis2 = FromSNP)) +
  geom_alluvium(aes(fill = ToFH), width = 1/18) +
  geom_stratum(width = 1/18, fill = 'transparent', color = "grey") +
  #geom_label(stat = "stratum", aes(label = after_stat(stratum)),family='sans',size = 8) +
  scale_fill_manual(values = cols)+
  theme(legend.position = 'none',
        text = element_text(family = 'sans'),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = 'transparent',colour = 'transparent'),
        panel.background = element_rect(fill = 'transparent',colour = 'transparent'),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())
ggsave(plot = p,filename = "E:/QiGA/document/CottonLab/FH/ResultsAllVariants/Phy/Group/Group_compare.pdf",width = 3,height = 8)
