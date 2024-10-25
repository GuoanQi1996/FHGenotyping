# plot
library(ggplot2)
library(ggsci)
dt = read.table("./HapStat/Plot_data.FHNumber.txt",header = T,sep='\t')
dt = dt[order(dt$Group),]
gsum = aggregate(dt$Number,list(dt$Group),sum)
gsum = rbind(gsum,gsum)
gsum = gsum[order(gsum$Group.1),]
dt$Proportion = dt$Number/gsum$x*100
dt$Number = dt$Number/10000
dt$Type = factor(dt$Type,levels = c("A-subgenome","D-subgenome"))
dt$Group = factor(dt$Group,levels = c("G6","G5","G4","G3","G2"))

col = pal_npg(alpha = 0.8)(2)
names(col) = c("A-subgenome","D-subgenome")
p = 
  ggplot(dt, aes(fill = factor(Type,levels = c("D-subgenome","A-subgenome")),y = Group, x = Number))+
  geom_bar(position = "stack", stat = "identity")+
  geom_text(aes(label = paste(round(Proportion, 1),"%")),
            position = position_stack(vjust = 0.5))+
  xlab(expression("FH number"~(phantom() %*% 10^4)))+
  ylab("")+
  scale_fill_manual(values = col)+
  scale_x_continuous(expand = c(0.02,0.02))+
  theme_bw()+
  theme(legend.position = c(0.8,0.15),
        legend.title = element_blank(),
        legend.background = element_rect(fill='transparent',colour = NA),
        text = element_text(color='black',family='sans'),
        axis.title.x = element_text(color='black',family='sans',size = 10),
        axis.text = element_text(color='black',family='sans',size = 10),
        axis.title.y = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = 'grey90'),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent',colour = NA))
ggsave("FH_Number.pdf",p,width = 5,height = 4)

# Eh and FH number
setwd("E:/QiGA/document/CottonLab/FH")
library(data.table)
library(minpack.lm)
FHap_All = fread("./ResultsAllVariants/FHAP_whole.GeneSTAT.txt",header=T)

model_lm <- lm(HapN ~ I(EH^2), data = FHap_All)
summary(model_lm)

library(ggplot2)
library(ggpmisc)
p = 
  ggplot(data = FHap_All,mapping = aes(x = EH, y = HapN))+
  geom_point()+
  xlab(expression(paste("Shannon's equitability (",italic(E[H]),")",sep = "")))+
  ylab("FH number")+
  stat_poly_line(formula = y ~ poly(x, 2, raw = TRUE)) +
  stat_poly_eq(formula = y ~ poly(x, 2, raw = TRUE), use_label(c("eq","adj.R2")))+
  theme_bw()+
  theme(text = element_text(family = 'sans',size = 10),
        axis.text = element_text(family = 'sans',size = 8))
ggsave("FH.vs.EH.pdf",p,width = 5,height = 4.9)

# chromosomal - level diversity
setwd("E:/QiGA/document/CottonLab/FH")
library(data.table)
library(ggplot2)
library(ggsci)
library(scales)
library(agricolae)
DT = read.table("./ResultsAllVariants/FHAP_whole.GeneSTAT.txt",header=T,sep = '\t')

gff = fread("../Cotton_GO/TM-1_V2.1.gene.CDS.gff",header=F,stringsAsFactors = F,data.table = F)
allGene = gff[which(gff$V3=='gene'),]
allGene$V9 = sub("ID=","",sub(";","",allGene$V9,fixed = T),fixed = T)
allGene = data.frame(Gene = allGene$V9)
houseKeeping = read.table("HouseKeeping_gene.list",header=F)
allVariedGene = data.frame(Gene = setdiff(allGene$Gene,houseKeeping$V1))

DT1 = merge(DT,allGene,'Gene',all.y=T)
DT1$HapN[which(is.na(DT1$HapN))] = 0
DT1$EH[which(is.na(DT1$EH))] = 0
DT1$SubGenome = "A"
DT1$SubGenome[grep("GH_A",DT1$Gene)] = "A"
DT1$SubGenome[grep("GH_D",DT1$Gene)] = "D"
DT1$CHR = sub("H_","",sapply(strsplit(DT1$Gene,"G",fixed = T),'[[',2),fixed = T)
DT1$HapN = DT1$HapN/1000
cols = c("#EFC000CC","#0073c2CC")
names(cols) = c("A","D")

SSDT1 = aggregate(list(DT1$EH),list(DT1$CHR),mean)
SSDT2 = aggregate(list(DT1$HapN),list(DT1$CHR),sum)
SSDT = merge(SSDT1,SSDT2,'Group.1')
colnames(SSDT) = c("CHR","EH","FH")
#SSDT$FH = SSDT$FH/1000
SSDT$SubGenome = substr(SSDT$CHR,1,1)

variance = aov(EH~CHR,data = DT1)
MC = LSD.test(variance,"CHR", p.adj="bonferroni")#结果显示：标记字母法out$group
MC = MC[["groups"]]
MC$CHR = rownames(MC)
MC = MC[order(MC$CHR,decreasing = F),]


FH_min_max = range(SSDT$FH)
EH_min_max = range(SSDT$EH)
x_min_max = 1.05*range(SSDT$FH)
p=
  ggplot(SSDT,aes(x = CHR))+
  geom_col(aes(y = FH, fill = SubGenome))+
  scale_fill_manual(name="",values = cols)+
  geom_line(aes(y = rescale(EH,x_min_max), group=1),linewidth=1.5,color='red') +
  geom_point(aes(y = rescale(EH,x_min_max)),shape=21,fill="white",size=2)+
  lims(y=x_min_max)+
  scale_y_continuous(breaks=breaks_pretty(5),sec.axis = sec_axis( ~rescale(.,EH_min_max),name = expression(paste("Avg. ",italic(E[H]),sep = ""))))+
  labs(y = expression(paste("FH number (", 10^3, ")"),sep=""))+
  theme_bw() +
  geom_text(data=MC,
            aes(x=CHR,y=rescale(EH,x_min_max)*1.1,label=groups))+
  theme(text = element_text(color='black',family='sans',size = 16),
        axis.title.x =  element_blank(),
        axis.text = element_text(color='black',family='sans',size = 14),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.title = element_blank(),
        legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent',colour = NA))
p
ggsave("FH_EH_chrLevel_SumFH.pdf",plot = p,width = 10,height = 4)

SSDT1 = aggregate(list(DT1$EH),list(DT1$CHR),mean)
SSDT2 = aggregate(list(DT1$HapN),list(DT1$CHR),mean)
SSDT = merge(SSDT1,SSDT2,'Group.1')
colnames(SSDT) = c("CHR","EH","FH")
SSDT$SubGenome = substr(SSDT$CHR,1,1)
FH_min_max = range(SSDT$FH)
EH_min_max = range(SSDT$EH)
x_min_max = 1.05*range(SSDT$FH)
p=
  ggplot(SSDT,aes(x = CHR))+
  geom_col(aes(y = FH, fill = SubGenome))+
  scale_fill_manual(name="",values = cols)+
  geom_line(aes(y = rescale(EH,x_min_max), group=1),linewidth=1.5,color='red') +
  geom_point(aes(y = rescale(EH,x_min_max)),shape=21,fill="white",size=2)+
  lims(y=x_min_max)+
  scale_y_continuous(breaks=breaks_pretty(5),sec.axis = sec_axis( ~rescale(.,EH_min_max),name = expression(paste("Avg. ", italic(E[H]),sep = ""))))+
  labs(y = paste("Avg. FH number"))+
  theme_bw() +
  geom_text(data=MC,
            aes(x=CHR,y=rescale(EH,x_min_max)*1.1,label=groups))+
  theme(text = element_text(color='black',family='sans',size = 16),
        axis.title.x =  element_blank(),
        axis.text = element_text(color='black',family='sans',size = 14),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.title = element_blank(),
        legend.position = 'none',
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent',colour = NA))
p
ggsave("FH_EH_chrLevel_AveFH.pdf",plot = p,width = 10,height = 4)

# EH and chromosome length 
SSDT1 = aggregate(list(DT1$EH),list(DT1$CHR),mean)
SSDT2 = aggregate(list(DT1$HapN),list(DT1$CHR),sum)
SSDT = merge(SSDT1,SSDT2,'Group.1')
colnames(SSDT) = c("CHR","EH","FH")
SSDT$SubGenome = substr(SSDT$CHR,1,1)
chrLength = read.table("E:/QiGA/document/CottonLab/Introgression/IntroDetect/TM-1_V2.1.CHR.size",header=F)
colnames(chrLength) = c("CHR","Start","End")
SSDT = merge(SSDT,chrLength,"CHR")
SSDT$End = SSDT$End/1000000

as.data.frame(table(DT1$CHR)) -> temp
colnames(temp) = c('CHR','Genes')
SSDT= merge(SSDT,temp,'CHR')
#SSDT$FH = SSDT$FH/1000
library(ggplot2)
library(ggpmisc)
p1 = 
  ggplot(data = SSDT,mapping = aes(x = End, y = FH))+
  geom_point()+
  xlab("Chromosome length (Mb)")+
  ylab(expression(paste("FH number (",10^3,")"),sep=""))+
  stat_poly_line(formula = y ~ poly(x, raw = TRUE)) +
  stat_poly_eq(formula = y ~ poly(x, raw = TRUE), use_label(c("eq","adj.R2")))+
  theme_bw()+
  theme(text = element_text(family = 'sans',size = 10,color = 'black'),
        axis.text = element_text(family = 'sans',size = 9),
        plot.margin = unit(c(1,1,1,1),'lines'),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent',colour = NA))
p1
ggsave("FH_vs_chrLength.pdf",plot = p1,width = 4.4,height = 4)

p2 = 
  ggplot(data = SSDT,mapping = aes(x = Genes, y = FH))+
  geom_point()+
  xlab("The number of genes")+
  ylab(expression(paste("FH number (",10^3,")"),sep=""))+
  stat_poly_line(formula = y ~ poly(x, raw = TRUE)) +
  stat_poly_eq(formula = y ~ poly(x, raw = TRUE), use_label(c("eq","adj.R2")))+
  theme_bw()+
  theme(text = element_text(family = 'sans',size = 10,color = 'black'),
        axis.text = element_text(family = 'sans',size = 9),
        plot.margin = unit(c(1,1,1,1),'lines'),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent',colour = NA))
p2
ggsave("FH_vs_geneCount.pdf",plot = p2,width = 4.4,height = 4)

p2 = 
  ggplot(data = SSDT,mapping = aes(x = End, y = EH))+
  geom_point()+
  ylab(expression(paste(italic(E[H]),sep = "")))+
  xlab("Chromosome length (Mb)")+
  stat_poly_line(formula = y ~ poly(x, raw = TRUE)) +
  stat_poly_eq(formula = y ~ poly(x, raw = TRUE), use_label(c("eq","adj.R2")))+
  theme_bw()+
  theme(text = element_text(family = 'sans',size = 10,color = 'black'),
        plot.margin = unit(c(1,1,1,1),'lines'),
        axis.text = element_text(family = 'sans',size = 9),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent',colour = NA))
ggsave("EH_vs_chrLength.pdf",plot = p2,width = 4.4,height = 4)

PL = list(p1,p2)
P = cowplot::plot_grid(plotlist = PL,align = 'hv',nrow = 1)
ggsave("EH_FH_vs_chrLength.pdf",width = 10,height = 4.5)

p = 
  ggplot(data = SSDT,mapping = aes(x = FH, y = EH))+
  geom_point()+
  ylab(expression(paste(italic(E[H]),sep = "")))+
  xlab(expression(paste("FH number (", 10^3, ")"),sep=""))+
  stat_poly_line(formula = y ~ poly(x, raw = TRUE)) +
  stat_poly_eq(formula = y ~ poly(x, raw = TRUE), use_label(c("eq","R2")))+
  theme_bw()+
  theme(text = element_text(family = 'sans',size = 16,color = 'black'),
        axis.text = element_text(family = 'sans',size = 14),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent',colour = NA))
ggsave("FH_vs_EH.pdf",plot = p,width = 4.2,height = 4)

DT1$SubGenome = factor(DT1$SubGenome,levels = c("A","D"))
#DT1$HapN = DT1$HapN/1000
cols = c("#EFC000CC","#0073c2CC")
names(cols) = c("A","D")
p3 = 
  ggplot(DT1, aes(x = SubGenome, y = EH, fill = SubGenome)) +
  geom_violin(trim = TRUE, scale = "width") +
  scale_fill_manual(name="",values = cols)+
  #geom_boxplot(width=0.05,position=position_dodge(0.9),fill = 'white', color = "black",outliers = F)+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  scale_x_discrete(name="",labels=c("At","Dt"))+
  #scale_fill_npg(alpha = 0.8)+
  #geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
  labs(y = expression(paste(italic(E[H]),sep = ""))) +
  theme_bw() +
  theme(text = element_text(color='black',family='sans',size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text = element_text(color='black',family='sans',size = 14),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.title = element_blank(),
        legend.position = 'none',
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent',colour = NA))
ggsave("EH_SubGenomeLevel.pdf",plot = p3,width = 5,height = 2)

p4 = 
  ggplot(DT1, aes(x = SubGenome, y = HapN, fill = SubGenome)) +
  geom_violin(trim = TRUE, scale = "width") +
  scale_fill_manual(name="",values = cols)+
  #geom_boxplot(width=0.05,position=position_dodge(0.9),fill = 'white', color = "black",outliers = F)+
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.2))+
  scale_x_discrete(name="",labels=c("At","Dt"))+
  #scale_fill_npg(alpha = 0.8)+
  #geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
  labs(y = expression(paste("FH number (", 10^3, ")"),sep="")) +
  theme_bw() +
  theme(text = element_text(color='black',family='sans',size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(color='black',family='sans',size = 14),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.title = element_blank(),
        legend.position = 'none',
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent',colour = NA))
ggsave("FH_SubGenomeLevel.pdf",plot = p4,width = 5,height = 2)

PL = list(p3,NULL,p4)
P = cowplot::plot_grid(plotlist = PL,align = 'hv',nrow = 3,rel_heights = c(2,-0.2,2))
P
ggsave("EH_FH_SubGenomeLevel.pdf",width = 5,height = 4.5)
