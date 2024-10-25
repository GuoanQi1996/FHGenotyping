library(agricolae)
library(ggsci)
dt1 = read.table("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/Phy/PIE/A06_FH_Block.txt",header=T,sep='\t',comment.char = "")
dt2 = read.table("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/Phy/PIE/A08_FH_Block.txt",header=T,sep='\t',comment.char = "")

DT = rbind(dt1,dt2)
DT$TYPE = paste0(DT$Type,"-",DT$HaplotypeBlock)
DT$Sample = sub("0_","",DT$Sample,fixed=T)
pheList = list.files("E:/QiGA/document/CottonLab/XLD/Phe/","^EMMAX*.*txt",full.names = T)

corresDT = data.frame(HapBlock = unique(DT$TYPE))
corresDT$HapBlockAlt = c("A06-Hapblock 1","A06-Hapblock 2","A06-Hapblock 4","A06-Hapblock 3",
                         "A08-Hapblock 1","A08-Hapblock 2","A08-Hapblock 4","A08-Hapblock 3","A08-Hapblock 5","A08-Hapblock 6")
col = pal_npg(alpha = 0.8)(10)
names(col) = corresDT$HapBlockAlt
SUMDT = data.frame()
for(i in 1:length(pheList)){
  phe = pheList[i]
  pheDT = read.table(phe,header=F,sep = '\t',comment.char = "")
  pheShort = strsplit(phe,"/",fixed=T)[[1]][7]
  pheShort = strsplit(pheShort,".",fixed=T)[[1]][1]
  pheShort = sub("EMMAX_","",pheShort,fixed=T)
  
  pheDT = pheDT[,c(2,3)]
  pheDT$PHE = pheShort
  colnames(pheDT) = c("Sample","Value","Phenotype")
  
  Final = merge(pheDT,DT,'Sample')
  NA_Line = which(is.na(Final$Value))
  if(length(NA_Line)>0){
    Final = Final[-NA_Line,]
  }
  
  cat(pheShort)
  cat("\n")
  
  DTAve = aggregate(Final$Value,list(Final$TYPE),mean,na.rm=T)
  DTMedi = aggregate(Final$Value,list(Final$TYPE),median,na.rm=T)
  DTQuantile = aggregate(Final$Value,list(Final$TYPE),function(x){
    iqr = IQR(x)
    Q3 = quantile(x,probs=0.75)
    min(max(x), Q3 + 1.5 * iqr)
  })
  colnames(DTAve)[2] = "average"
  colnames(DTMedi)[2] = "median"
  colnames(DTQuantile)[2] = "quantile"
  DTStat = as.data.frame(table(Final$TYPE))
  DTStat = merge(DTStat,DTMedi,by.x = 'Var1',by.y = 'Group.1')
  DTStat = merge(DTStat,DTAve,by.x = 'Var1',by.y = 'Group.1')
  DTStat = merge(DTStat,DTQuantile,by.x = 'Var1',by.y = 'Group.1')
  DTStat = DTStat[order(DTStat$average,decreasing = F),]
  DTStat = merge(DTStat,corresDT,by.x = 'Var1',by.y = 'HapBlock',sort=F)
  xlabel = paste0(DTStat$HapBlockAlt,"\n(n=",DTStat$Freq,")")
  Final$TYPE = factor(Final$TYPE,levels = as.character(DTStat$Var1),labels = xlabel)
  DTStat$Var1 = factor(DTStat$Var1,levels = as.character(DTStat$Var1),labels = xlabel)
  
  variance = aov(Value~TYPE,data = Final)
  MC = LSD.test(variance,"TYPE", p.adj="bonferroni")#结果显示：标记字母法out$group
  MC = MC[["groups"]]
  MC$HapBlock = rownames(MC)
  MC = MC[order(MC$Value,decreasing = F),]
  print(MC)
  
  MC = merge(MC,DTStat,by.x = 'HapBlock',by.y = 'Var1',sort=F)
  MC$PHE = pheShort
  SUMDT = rbind(SUMDT,MC)
  colnames(MC)[1] = 'TYPE'

  cols_in_use = col[match(MC$HapBlockAlt,names(col))]
  names(cols_in_use) = as.character(MC$TYPE)
  p = ggplot(Final,aes(x=TYPE, y = Value, fill = TYPE))+
    geom_boxplot(width = 0.5,outlier.colour = NA)+
    scale_x_discrete(breaks = MC$TYPE[order(MC$Value,decreasing = T)])+
    scale_y_continuous(limits = c(quantile(Final$Value,0.01),max(Final$Value)*1.05))+
    xlab("")+
    ylab(pheShort)+
    theme_bw()+
    geom_text(data=MC,size = 3,
              aes(x=factor(TYPE),y=quantile*1.02,label=groups))+
    scale_fill_manual(values = cols_in_use)+
    coord_flip()+
    theme(legend.position = 'none',
          plot.margin = unit(c(0.5,1,0,1),'line'),
          plot.background = element_rect(fill = 'transparent',color = 'transparent'),
          panel.background = element_rect(fill = 'transparent',color = 'transparent'),
          text = element_text(family = 'sans',color = 'black'),
          axis.title.y = element_text(family = 'sans',size = 18,color = 'black'),
          axis.title.x = element_text(family = 'sans',size = 18,color = 'black'),
          axis.text.y = element_text(size=16,color = 'black'),
          axis.text.x = element_text(size=14,color = 'black'))
  ggsave(paste0("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/Haplotype/pheDEV/Phe_",pheShort,".hapDev.pdf"),p,width = 4,height = 7)
}
SUMDT = SUMDT[,-1]
SUMDT = merge(SUMDT,corresDT,by = 'HapBlockAlt', sort = F)
write.csv(SUMDT,"E:/QiGA/document/CottonLab/FH/ResultsAllVariants/Haplotype/pheDEV/Summary.csv",quote=F,col.names = T,row.names = F)
