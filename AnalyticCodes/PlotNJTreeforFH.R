# plot nj tree
library(ape)
library(aplot)
library(readxl)
library(ggplot2)
library(ggtree)
library(treeio)
library(ggsci)
library(ggstance)
library(ggnewscale)

TreeType = "NJ"

setwd("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/")
tree = read.newick(paste0("./Phy/FH.Var2.Factor.NJ.tree"))
tree$tip.label = sub("0_","",tree$tip.label)
sample.info = read.table("E:/QiGA/document/CottonLab/XLD/population_info/Detailed_Information.txt",header=T,quote="",comment.char = "@",sep='\t')[,c(1,11,14)]
tree.node = data.frame(VCF_ID=tree$tip.label)
sample.info = merge(sample.info,tree.node,sort=F,all.y=T)

#sample.info$`Data interval`[match(c("SC (CHN)","Central America","South America","North America"),sample.info$Zone1)]="Uncertain"
sample.info$`Data.interval`[which(sample.info$`Data.interval`=='')]='Uncertain'

sample_info = data.frame(Zone=sample.info$Zone1)
rownames(sample_info) = sample.info$VCF_ID

date_info = data.frame(Date=sample.info$`Data.interval`)
rownames(date_info) = sample.info$VCF_ID
date_info$Date[which(date_info$Date=="1950-1959" | date_info$Date=="1960-1969")] = "1950-1969"
date_info$Date[which(date_info$Date=="1970-1979" | date_info$Date=="1980-1989")] = "1970-1989"
date_info$Date[which(date_info$Date=="1990-1999" | date_info$Date=="2000-2009")] = "1990-2009"
date_info$Date = factor(date_info$Date,levels = c("Before 1900","1900-1949","1950-1969","1970-1989","1990-2009","2010-2019"))
#date_info$Date = factor(date_info$Date,levels = c("Before 1900","1900-1949","1950-1959","1960-1969","1970-1979","1980-1989","1990-1999","2000-2009","2010-2019"))
#date_info$Date = as.numeric(date_info$Date)

cols = c(`SWC`="#990080FF",`YZR`="#0099CCFF",`YER`="#BA6338FF",
         `NWC`="#339900FF",`SCL`="#8A419899",`NC`="#003399FF",
         USA="#02024999", Africa="#FF7F0EFF",Europe="#8F770099",
         `AL`="#8C564BFF",Asia="#008B4599",FSU="#B2474599",
         GB="#EFC00099")

# rm low X
rmList = tree[['tip.label']][grep("FL",tree[['tip.label']])]
rmList = c(rmList,tree[['tip.label']][grep("LGZ",tree[['tip.label']])])
rmList = c(rmList,tree[['tip.label']][grep("ZXL",tree[['tip.label']])])
tree = drop.tip(tree,rmList)
tree1 <- root(tree, outgroup = "DXMGB_GH1707", edgelabel = TRUE)
p1 = ggtree(tree1, branch.length="none", size=0.001)
p_data = p1$data

p_data = p_data[which(p_data$isTip==TRUE),]
p_data = p_data[order(p_data$y,decreasing = F),]
rmOutGB = p_data$label[which(p_data$y>87 & grepl("GB",p_data$label))]
tree1 = drop.tip(tree1,rmOutGB)
p1 = ggtree(tree1, branch.length="none", size=0.001)
p_data = p1$data

p_data = p_data[which(p_data$isTip==TRUE),]
p_data = p_data[order(p_data$y,decreasing = F),]

p_data = merge(p_data, sample.info,by.x = 'label',by.y = "VCF_ID",sort=F)
#colnames(p_data)[11] = "Zone"
#write.table(p_data,"Circular.CORE.rotate.setmissing.phy.structure.data.table",quote=F,col.names = T,row.names = F,sep='\t')
#3266
#tree1 <- root(tree, outgroup = "DXMGB_GH1707", edgelabel = TRUE)

pt = ggtree(tree1, size=0.01,layout='circular',color = 'grey40') +
  #layout_dendrogram()+
  theme_tree(strip.background = element_blank(),
             strip.text.x = element_blank(),bgcolor = NA)

#pt = rotate(pt,3721)
#pt = rotate(pt,4503)
#pt = rotate(pt,4503)
p_data = pt$data

p_data = p_data[which(p_data$isTip==TRUE),]
p_data = p_data[order(p_data$y,decreasing = F),]

p_data = merge(p_data, sample.info,by.x = 'label',by.y = "VCF_ID",sort=F)
colnames(p_data)[11] = "Zone"
write.table(p_data,"./Phy/Circular.CORE.rotate.FH.Var2.Factor.phy.structure.data.table",quote=F,col.names = T,row.names = F,sep='\t')

# p1 = pt %<+% sample.info +
#   aes(color=I(Zone1))+
#   geom_text(aes(label=node), hjust=-.3, size = 0.5)+
#   scale_color_manual(values=c(cols),name = "")
# ggsave(paste0("./Phy/CORE.circular.internalNode.pdf"),p1,width = 30,height = 30)

pt = gheatmap(pt, sample_info, offset=0, width=.06,color = NA,colnames = F,font.size = 1) +
  scale_fill_manual(values=c(cols),name = "",na.value = "grey70")
pt = pt + new_scale_fill()
pal_continuous = pal_gsea(palette = c("default"), n = length(unique(date_info$Date)), alpha = 1, reverse = FALSE)(length(unique(date_info$Date)))
#pal_continuous = pal_material("blue",n = length(unique(date_info$Date)),alpha = 0.8)(length(unique(date_info$Date)))
names(pal_continuous) = c("Before 1900","1900-1949","1950-1969","1970-1989","1990-2009","2010-2019")
pt = gheatmap(pt, date_info, offset=1800, width=.06,color = NA,colnames = F,font.size = 1) +
  scale_fill_manual(values=c(pal_continuous),name = "",na.value = "grey70")

pt = pt +  
  theme(plot.margin=margin(0,0,0,0),
        legend.position = 'none',
        panel.grid  = element_blank(),
        panel.border=element_blank(),
        strip.background = element_rect(fill='transparent'),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent',colour = NA),
        legend.background = element_rect(fill='transparent')#,
        # legend.spacing.y=unit(0.02, "cm"), 
        # legend.title=element_text(size=7.5), 
        # legend.text=element_text(size=5.5), 
        # legend.box.spacing=unit(0.02,"cm")
  )

pt = open_tree(pt,10)
ggsave(paste0("./Phy/Circular.CORE.noLegend.FH.Var2.Factor.phy.structure.pdf"),pt,width = 8,height = 8)
