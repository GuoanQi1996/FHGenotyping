# plot
# only 9 genes do not have 00000 haplotype
setwd("E:/QiGA/document/CottonLab/FH")
library(data.table)
library(pheatmap)
library(RColorBrewer)
library(ape)

group1 = data.frame(label=read.table("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/Phy/Group/G1.FH.sample.list.txt",comment.char = "")[,1])
group2 = data.frame(label=read.table("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/Phy/Group/G2.FH.sample.list.txt",comment.char = "")[,1])
group3 = data.frame(label=read.table("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/Phy/Group/G3.FH.sample.list.txt",comment.char = "")[,1])
group4 = data.frame(label=read.table("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/Phy/Group/G4.FH.sample.list.txt",comment.char = "")[,1])
group5 = data.frame(label=read.table("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/Phy/Group/G5.FH.sample.list.txt",comment.char = "")[,1])
group6 = data.frame(label=read.table("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/Phy/Group/G6.FH.sample.list.txt",comment.char = "")[,1])

group1$Group = "G1"
group2$Group = "G2"
group3$Group = "G3"
group4$Group = "G4"
group5$Group = "G5"
group6$Group = "G6"

group = rbind(group1,group2,group3,group4,group5,group6)
group$label = paste0("0_",group$label)

pop_info = read.table("E:/QiGA/document/CottonLab/XLD/population_info/Detailed_Information.txt",header=T,sep='\t',quote = "",comment.char = "@",encoding = 'utf_8')
pop_info = pop_info[,c("VCF_ID","Data.interval","Zone2")]
pop_info$VCF_ID = paste0("0_",pop_info$VCF_ID)
colnames(pop_info) = c("label","Cultivated Date","Source Region")

group = merge(group,pop_info,"label",sort=F)

anno.dt = data.frame(group[,c(2,4,3)])
rownames(anno.dt) = group$label
anno.dt$Group = factor(anno.dt$Group)
anno.dt$Cultivated.Date = factor(anno.dt$Cultivated.Date)
anno.dt$Source.Region = factor(anno.dt$Source.Region)
anno.dt = anno.dt[,c(1,2)]

hapStat = read.table("./ResultsAllVariants/FHAP_whole.HapSTAT.txt",header=T,colClasses = c("character","character","numeric"),stringsAsFactors = F)
hapStat$Gene = sapply(strsplit(hapStat$Name,"_HAP",fixed=T),'[[',1)

gff = fread("../Cotton_GO/TM-1_V2.1.gene.CDS.gff",header=F,stringsAsFactors = F,data.table = F)
allGene = gff[which(gff$V3=='gene'),]
allGene$V9 = sub("ID=","",sub(";","",allGene$V9,fixed = T),fixed = T)
allGene$V10 = floor((allGene$V4+allGene$V5)/2)
allGene = allGene[,c(9,1,10)]

hf_filter <- function(x){
  keep = names(which(sort(table(x),decreasing = T)>10))
  x[which(is.na(match(x,keep)))]=1
  as.integer(factor(x))
}
GT = fread("./ResultsAllVariants/FHAP_whole.type.2.txt", sep = '\t', data.table = F, stringsAsFactors = F)

geneIndex = GT[,1]
includeChr = sapply(strsplit(sapply(strsplit(geneIndex,"_",fixed=T),'[[',2),"G",fixed=T),'[[',1)
includeChrUniq = unique(includeChr)

#pal_continuous = c("grey",ggsci::pal_gsea(palette = c("default"), n = length(na.omit(unique(anno.dt$Cultivated.Date))), alpha = 0.8, reverse = FALSE)(length(na.omit(unique(anno.dt$Cultivated.Date)))))
#names(pal_continuous) = rev(c(as.character(na.omit(unique(anno.dt$Cultivated.Date))),"NA"))
cols = list(Group = c(G1="#EFC000CC",G2="#0099CCCC",G3="#BA6338CC",G4="#631879CC",G5="#339900CC",G6="#D43F3ACC"),
            Source.Region = c(`SWC (CHN)`="#990080CC",`YZR (CHN)`="#0099CCCC",`YER (CHN)`="#BA6338CC",`NWC (CHN)`="#339900CC",
                              `EC (CHN)`="#5DB1DDCC",`SC (CHN)`="#8A4198CC",`NC (CHN)`="#003399CC",CHN="#E66B65CC",
                              `Central America`="#D43F3ACC",`Wild`="#5CB85CCC", USA="#020249CC",
                              FSU="#DC0000CC",Africa="#FF7F0ECC",Uncertain="grey",`South America`="#E18727CC",`North America`="#EEA236CC",
                              Europe="#8F7700CC",Oceania="#17BECFCC",Asia="#008B45CC",GB="#EFC000CC"))

sumDT = data.frame()
for(i in 1:length(includeChrUniq)){
  chr = includeChrUniq[i]
  index = which(includeChr==chr)
  
  #orderS = read.table(paste0("E:/QiGA/document/bioinformatics/HuYan/splitCHR/CHR_", chr, ".majRef.cluster_order.txt"),header=F)
  #orderS$V1 = paste0("0_",orderS$V1)
  #GT1 = GT[index,c(1,as.numeric(na.omit(rev(match(orderS$V1,colnames(GT))))))]
  GT1 = GT[index,]
  GTP = as.matrix(GT1[,-1])
  
  numeric.var.vec = apply(GTP,1,var,na.rm=T)
  #quant.numeric.var = quantile(numeric.var.vec,probs = seq(0,1,0.01))
  GTP = GTP[which(numeric.var.vec>=2),]
  GTP = apply(GTP,2,as.integer)
  GT.C = apply(GTP,2,as.character)
  disMat = dist.gene(t(GT.C), method = 'pairwise')
  
  GTP = t(apply(GTP,1,hf_filter))
  
  # nj.tree3 = nj(disMat)
  # nj.tree3[["tip.label"]] = rownames(GT)
  
  colnames(GTP) = colnames(GT1)[-1]
  rownames(GTP) = GT1$Gene[which(numeric.var.vec>=2)]
  #rownames(GTP) = GT1$Gene
  maxLevel = max(GTP,na.rm = T)
  if(maxLevel>52){
    colours = c('#FBF7D8CC',ggsci::pal_igv(alpha = 0.8)(51),colorRampPalette(colors = brewer.pal(8,"Dark2")[1:8])(maxLevel-52))
  } else {
    colours = c('#FBF7D8CC',ggsci::pal_igv(alpha = 0.8)(maxLevel-1))
  }
  
  x_lab = rep("",nrow(GTP))
  pos = allGene$V10[match(rownames(GTP),allGene$V9)]
  ticks1 = floor(c(seq(1,nrow(GTP),nrow(GTP)/10),nrow(GTP)))
  ticks2 = floor(c(seq(1,nrow(GTP),nrow(GTP)/10),nrow(GTP)-8))
  x_lab[ticks2]=floor(as.numeric(c(0,pos[ticks1[-1]]))/1000000)
  
  # GTP = imputeTS::na_replace(GTP,1)
  # anno.dt1 = anno.dt[which(!is.na(match(colnames(GTP),rownames(anno.dt)))),c(3,2,1)]
  # column_anno = HeatmapAnnotation(df = anno.dt1,col = cols, na_col = 'grey80',
  #                                 show_annotation_name = F,show_legend = F,
  #                                 which = 'row',simple_anno_size = unit(1,"cm"))
  # p=Heatmap(t(GTP),col = colours,na_col = "grey",border = F,
  #         cluster_rows = T,cluster_columns = F,show_row_dend = F,show_row_names = F,
  #         column_labels = x_lab,column_names_rot = 0,column_names_gp = grid::gpar(fontsize = 26),
  #         show_heatmap_legend = F,left_annotation = column_anno,
  #         raster_quality = 8,raster_by_magick = T)
  # png("ComplexHeatmap.CHRA6.Q8.png",width = 16*300,height = 8*300,res = 300)
  # pp = draw(p)
  # dev.off()
  # orderS = colnames(GTP)[row_order(pp)]
  
  p=
    pheatmap(t(GTP),color = colours,abreaks=c(1:maxLevel)-0.5,na_col = "#d8d8d899",font.family= "Arial",clustering_method = 'average',
             cluster_cols = F,clustering_distance_rows = disMat, annotation_names_row = F, scale = 'none',treeheight_row = 0,annotation_colors = cols,labels_col = x_lab,
             show_rownames = F,show_colnames = T,angle_col = 0,fontsize_col = 14,border_color = NA,legend = F,annotation_row = anno.dt, annotation_legend = F,
             filename = paste0("./ResultsAllVariants/Haplotype/CHR_", chr, ".FHAP.pheatmap.reCluster.var2.tiff"),width = 8,height = 3)
  
  x_lab = rep("",nrow(GTP))
  orderS = p$tree_row[['labels']][p$tree_row[["order"]]]
  write.table(orderS,paste0("./ResultsAllVariants/Haplotype/CHR_",chr,".FHAP.cluster_order.var2.txt"),quote=F,col.names = F,row.names = F)
  
  GTP = GTP[,orderS]
  p=
    pheatmap(t(GTP),color = colours,abreaks=c(1:maxLevel)-0.5,na_col = "#d8d8d899",font.family= "Arial",clustering_method = 'average',
             cluster_cols = F,cluster_rows = F,annotation_names_row = F,scale = 'none',treeheight_row = 0,annotation_colors = cols,labels_col = x_lab,
             show_rownames = F,show_colnames = T,angle_col = 0,fontsize_col = 14,border_color = NA,legend = F,annotation_row = anno.dt, annotation_legend = F,
             filename = paste0("./ResultsAllVariants/Haplotype/CHR_", chr, ".FHAP.pheatmap.reCluster.var2.noX.tiff"),width = 8,height = 3)
    
  orderS = data.frame(label = orderS)
  orderS = merge(orderS,group,'label',sort=F)
  orderS = cbind(No = 1:nrow(orderS),Chromosome = chr,orderS)
  sumDT = rbind(sumDT,orderS)
}
colnames(sumDT)[3] = "Label"
write.csv(sumDT,"./ResultsAllVariants/Haplotype/pop_info_with_groupG_perCHR.csv",quote = F,row.names = F)

# cha zhi
setwd("E:/QiGA/document/CottonLab/FH")
library(ggplot2)
library(reshape2)
library(data.table)
library(imputeTS)
library(pheatmap)
library(RColorBrewer)

group1 = data.frame(label=read.table("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/Phy/Group/G1.FH.sample.list.txt",comment.char = "")[,1])
group2 = data.frame(label=read.table("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/Phy/Group/G2.FH.sample.list.txt",comment.char = "")[,1])
group3 = data.frame(label=read.table("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/Phy/Group/G3.FH.sample.list.txt",comment.char = "")[,1])
group4 = data.frame(label=read.table("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/Phy/Group/G4.FH.sample.list.txt",comment.char = "")[,1])
group5 = data.frame(label=read.table("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/Phy/Group/G5.FH.sample.list.txt",comment.char = "")[,1])
group6 = data.frame(label=read.table("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/Phy/Group/G6.FH.sample.list.txt",comment.char = "")[,1])

group1$Group = "G1"
group2$Group = "G2"
group3$Group = "G3"
group4$Group = "G4"
group5$Group = "G5"
group6$Group = "G6"

group = rbind(group1,group2,group3,group4,group5,group6)
group$label = paste0("0_",group$label)

pop_info = read.table("E:/QiGA/document/CottonLab/XLD/population_info/Detailed_Information.txt",header=T,sep='\t',quote = "",comment.char = "@",encoding = 'utf_8')
pop_info = pop_info[,c("VCF_ID","Data.interval","Zone2")]
pop_info$VCF_ID = paste0("0_",pop_info$VCF_ID)
colnames(pop_info) = c("label","Cultivated Date","Source Region")

group = merge(group,pop_info,"label",sort=F)

anno.dt = data.frame(group[,c(1,2,4)])
rownames(anno.dt) = group$label
anno.dt$Group = factor(anno.dt$Group)
#anno.dt$Cultivated.Date = factor(anno.dt$Cultivated.Date)
anno.dt$Source.Region = factor(anno.dt$Source.Region)
anno.dt = anno.dt[,c(2,3)]

gff = fread("../Cotton_GO/TM-1_V2.1.gene.CDS.gff",header=F,stringsAsFactors = F,data.table = F)
allGene = gff[which(gff$V3=='gene'),]
allGene$V9 = sub("ID=","",sub(";","",allGene$V9,fixed = T),fixed = T)
allGene$V10 = floor((allGene$V4+allGene$V5)/2)
allGene = allGene[,c(9,4,5)]
colnames(allGene) = c("Gene","XMIN","XMAX")

hf_filter <- function(x){
  keep = names(which(sort(table(x),decreasing = T)>10))
  x[which(is.na(match(x,keep)))]=1
  as.integer(factor(x))
}
GT_FH = fread("./ResultsAllVariants/FHAP_whole.type.2.txt", sep = '\t', data.table = F, stringsAsFactors = F)

geneIndex = GT_FH[,1]
includeChr = sapply(strsplit(sapply(strsplit(geneIndex,"_",fixed=T),'[[',2),"G",fixed=T),'[[',1)
includeChrUniq = unique(includeChr)

#pal_continuous = c("grey",ggsci::pal_gsea(palette = c("default"), n = length(na.omit(unique(anno.dt$Cultivated.Date))), alpha = 0.8, reverse = FALSE)(length(na.omit(unique(anno.dt$Cultivated.Date)))))
#names(pal_continuous) = rev(c(as.character(na.omit(unique(anno.dt$Cultivated.Date))),"NA"))
cols = list(Group = c(G1="#EFC000CC",G2="#0099CCCC",G3="#BA6338CC",G4="#631879CC",G5="#339900CC",G6="#D43F3ACC"),
            Source.Region = c(`SWC (CHN)`="#990080CC",`YZR (CHN)`="#0099CCCC",`YER (CHN)`="#BA6338CC",`NWC (CHN)`="#339900CC",
                              `EC (CHN)`="#5DB1DDCC",`SC (CHN)`="#8A4198CC",`NC (CHN)`="#003399CC",CHN="#E66B65CC",
                              `Central America`="#D43F3ACC",`Wild`="#5CB85CCC", USA="#020249CC",
                              FSU="#DC0000CC",Africa="#FF7F0ECC",Uncertain="grey",`South America`="#E18727CC",`North America`="#EEA236CC",
                              Europe="#8F7700CC",Oceania="#17BECFCC",Asia="#008B45CC",GB="#EFC000CC"))

sumDT1 = data.frame()
sumDT2 = data.frame()
type = 'Ref'
for(i in 1:length(includeChrUniq)){
  chr = includeChrUniq[i]
  index = which(includeChr==chr)
  
  GT_SNP = fread(paste0("../../bioinformatics/HuYan/splitCHR/",chr, ".", type, ".thin05.traw"),header=T)
  rmList = c(grep("FL",colnames(GT_SNP)),
             grep("LGZ",colnames(GT_SNP)),
             grep("ZXL",colnames(GT_SNP)))
  GT_SNP = as.data.frame(as.data.frame(GT_SNP)[,-c(1,3:6,rmList)])
  
  # orderS1 = read.table(paste0("E:/QiGA/document/bioinformatics/HuYan/splitCHR/CHR_", chr, ".", type, ".cluster_order.txt"), header=F)
  # orderS1$V1 = paste0("0_",orderS1$V1)
  # GT_SNP = GT_SNP[,c(1,as.numeric(na.omit(rev(match(orderS1$V1,colnames(GT_SNP))))))]
  highHet = apply(GT_SNP[,2:ncol(GT_SNP)],1,function(x){
    length(which(x==1))/length(which(!is.na(x)))
  })
  lowVar = apply(GT_SNP[,2:ncol(GT_SNP)],1,function(x){
    var(x,na.rm=T)
  })
  GT_SNP = GT_SNP[which(highHet<=0.2 | lowVar<=0),]
  GT_SNP = na_replace(GT_SNP,2)
  SNP_POS = data.frame(TYPE = "SNP",ID = GT_SNP$SNP, POS = sapply(strsplit(GT_SNP$SNP,"_",fixed=T),'[[',2))
  
  orderS2 = read.table(paste0("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/Haplotype/CHR_", chr, ".FHAP.cluster_order.var2.txt"), header=F,comment.char = "$")
  GT_FH_CHR = GT_FH[index,c(1,as.numeric(na.omit(rev(match(orderS2$V1,colnames(GT_FH))))))]
  rownames(GT_FH_CHR) = GT_FH_CHR$Gene
  nameRow = rownames(GT_FH_CHR)
  nameCol = colnames(GT_FH_CHR)
  GT_FH_CHR_temp = t(apply(GT_FH_CHR[,-1],1,hf_filter))
  GT_FH_CHR = data.frame(Gene = nameRow, GT_FH_CHR_temp)
  colnames(GT_FH_CHR) = nameCol
  rownames(GT_FH_CHR) = nameRow
  
  FH_POS = merge(data.frame(Gene = GT_FH_CHR$Gene),allGene,"Gene",sort=F)
  FH_POS = data.frame(TYPE="FH",ID = FH_POS$Gene,POS = FH_POS$XMIN)
  combinedDT = rbind(SNP_POS,FH_POS)
  combinedDT$POS = as.numeric(combinedDT$POS)
  combinedDT = combinedDT[order(combinedDT$POS),]
  control = data.frame(ID = combinedDT$ID)
  
  GT_SNP1 = merge(control,GT_SNP,by.x="ID",by.y="SNP",sort=F,all = T)
  GT_SNP1 = GT_SNP1[match(control$ID,GT_SNP1$ID),]
  
  GT_SNP_FINAL = na_locf(GT_SNP1,option = 'locf')
  
  GTP = as.matrix(GT_SNP_FINAL[,-1])
  nameCol = colnames(GTP)
  GTP = GTP[,rev(nameCol)]
  
  x_lab = rep("",nrow(GTP))
  pos = combinedDT$POS
  ticks1 = floor(c(seq(1,nrow(GTP),nrow(GTP)/10),nrow(GTP)))
  ticks2 = floor(c(seq(1,nrow(GTP),nrow(GTP)/10),nrow(GTP)-60))
  x_lab[ticks2]=floor(as.numeric(c(0,pos[ticks1[-1]]))/1000000)
  
  p=
    pheatmap(t(GTP),color = c("#9013fe99","#f8e71c99","#d8d8d899"),abreaks=c(-0.5,0.5,1.5,2.5),na_col = "#d8d8d899",
             font.family= "Arial",clustering_method = 'average',
             cluster_cols = F,cluster_rows = T,annotation_names_row = F,scale = 'none',treeheight_row = 0,annotation_colors = cols,annotation_legend = F,
             show_rownames = F,show_colnames = T,angle_col = 0,fontsize = 14, border_color = NA, legend = F, annotation_row = anno.dt, labels_col = x_lab,
             filename = paste0("./ResultsAllVariants/Haplotype/CHR_", chr, ".", type, "_SNP.pheatmap.tiff"),width = 8,height = 3)
  
  write.table(p[["tree_row"]][["labels"]][p[["tree_row"]][["order"]]],
              paste0("./ResultsAllVariants/Haplotype/CHR_",chr,".SNP.cluster_order.txt"),
              quote=F,col.names = F,row.names = F)
  orderS1 = read.table(paste0("./ResultsAllVariants/Haplotype/CHR_",chr,".SNP.cluster_order.txt"),header=F,comment.char = "$")
  GTP = GTP[,orderS1$V1]
  p=
    pheatmap(t(GTP),color = c("#9013fe99","#f8e71c99","#d8d8d899"),abreaks=c(-0.5,0.5,1.5,2.5),na_col = "#d8d8d899",
             font.family= "Arial",clustering_method = 'average',
             cluster_cols = F,cluster_rows = F,annotation_names_row = F,scale = 'none',treeheight_row = 0,annotation_colors = cols,annotation_legend = F,
             show_rownames = F,show_colnames = F,angle_col = 0,fontsize = 14, border_color = NA, legend = F, annotation_row = anno.dt,
             filename = paste0("./ResultsAllVariants/Haplotype/CHR_", chr, ".", type, "_SNP.pheatmap.noX.tiff"),width = 8,height = 3)
  
  #orderS1 = colnames(GTP)

  orderS1 = data.frame(label = orderS1$V1)
  orderS1 = merge(orderS1,group,'label',sort=F)
  orderS1 = cbind(No = 1:nrow(orderS1),Chromosome = chr,orderS1)
  sumDT1 = rbind(sumDT1,orderS1)
  
  
  GT_FH_CHR1 = merge(control,GT_FH_CHR,by.x="ID",by.y="Gene",sort=F,all = T)
  GT_FH_CHR1 = GT_FH_CHR1[match(control$ID,GT_FH_CHR1$ID),]
  
  GT_FH_CHR_FINAL = na_locf(GT_FH_CHR1,option = 'locf')
  GTP = as.matrix(GT_FH_CHR_FINAL[,-1])
  nameCol = colnames(GTP)
  GTP = GTP[,rev(nameCol)]
  
  maxLevel = max(GTP,na.rm = T)
  if(maxLevel>52){
    colours = c('#FBF7D8CC',ggsci::pal_igv(alpha = 0.8)(51),colorRampPalette(colors = brewer.pal(8,"Dark2")[1:8])(maxLevel-52))
  } else {
    colours = c('#FBF7D8CC',ggsci::pal_igv(alpha = 0.8)(maxLevel-1))
  }
  
  x_lab = rep("",nrow(GTP))
  pos = combinedDT$POS
  ticks1 = floor(c(seq(1,nrow(GTP),nrow(GTP)/10),nrow(GTP)))
  ticks2 = floor(c(seq(1,nrow(GTP),nrow(GTP)/10),nrow(GTP)-60))
  x_lab[ticks2]=floor(as.numeric(c(0,pos[ticks1[-1]]))/1000000)
  
  # p=
  #   pheatmap(t(GTP),color = colours,na_col = "#d8d8d899",
  #            font.family= "Arial",clustering_method = 'average',
  #            cluster_cols = F,cluster_rows = T,annotation_names_row = F,scale = 'none',treeheight_row = 0,annotation_colors = cols,annotation_legend = F,
  #            show_rownames = F,show_colnames = T,angle_col = 0,fontsize = 14, border_color = NA, legend = F, annotation_row = anno.dt, labels_col = x_lab,
  #            filename = paste0("./Results/splitCHR/CHR_", chr, ".FHAP.pheatmap.NEW.tiff"),width = 8,height = 3)
  # write.table(p[["tree_row"]][["labels"]][p[["tree_row"]][["order"]]],
  #             paste0("./Results/splitCHR/CHR_",chr,".FHAP.cluster_order.NEW.txt"),
  #             quote=F,col.names = F,row.names = F)
  orderS2 = read.table(paste0("./ResultsAllVariants/Haplotype/CHR_",chr,".FHAP.cluster_order.var2.txt"),header=F,comment.char = "$")
  GTP = GTP[,orderS2$V1]
  p=
    pheatmap(t(GTP),color = colours,na_col = "#d8d8d899",
             font.family= "Arial",clustering_method = 'average',
             cluster_cols = F,cluster_rows = F,annotation_names_row = F,scale = 'none',treeheight_row = 0,annotation_colors = cols,annotation_legend = F,
             show_rownames = F,show_colnames = F,angle_col = 0,fontsize = 14, border_color = NA, legend = F, annotation_row = anno.dt,
             filename = paste0("./ResultsAllVariants/Haplotype/CHR_", chr, ".FHAP.CutX.pheatmap.InsertVar.tiff"),width = 8,height = 3)
  
  #orderS2 = colnames(GTP)
  
  orderS2 = data.frame(label = orderS2$V1)
  orderS2 = merge(orderS2,group,'label',sort=F)
  orderS2 = cbind(No = 1:nrow(orderS2),Chromosome = chr,orderS2)
  sumDT2 = rbind(sumDT2,orderS2)
}
colnames(sumDT1)[3] = "Label"
colnames(sumDT2)[3] = "Label"
write.csv(sumDT1,"./ResultsAllVariants/Haplotype/pop_info_with_groupG_perCHR.SNP.csv",quote = F,row.names = F)
write.csv(sumDT2,"./ResultsAllVariants/Haplotype/pop_info_with_groupG_perCHR.FHAP.csv",quote = F,row.names = F)
