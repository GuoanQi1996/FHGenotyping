# varied genes / extremely varied genes / house-keeping genes
gff = fread("../Cotton_GO/TM-1_V2.1.gene.CDS.gff",header=F,stringsAsFactors = F,data.table = F)
allGene = gff[which(gff$V3=='gene'),]
allGene$V9 = sub("ID=","",sub(";","",allGene$V9,fixed = T),fixed = T)
houseKeeping = setdiff(allGene$V9,rownames(GT))
ExtremelyVaried = setdiff(setdiff(allGene$V9,houseKeeping),rownames(GT_major))
varied = setdiff(allGene$V9,c(houseKeeping,ExtremelyVaried))
write.table(ExtremelyVaried,"ExtremelyVaried_gene.list",quote=F,col.names = F,row.names = F)
write.table(varied,"Varied_gene.list",quote=F,col.names = F,row.names = F)
library(clusterProfiler)
library(enrichplot)
go_gene = read.table("E:/QiGA/document/CottonLab/Cotton_GO/TM-1_final.agriGO.txt")
go_desc = read.csv("E:/QiGA/document/CottonLab/Cotton_GO/GO_desc_all.csv")
go_desc$nameform = substr(go_desc$name, 1, 100)
term2gene = data.frame(TERM = go_gene$V2, gene = go_gene$V1)
term2dec = data.frame(TERM = go_desc$Goid, NAME = go_desc$nameform)
en_go=enricher(houseKeeping,TERM2GENE = term2gene,TERM2NAME = term2dec,pvalueCutoff = 0.05,pAdjustMethod = "BH")
write.table(en_go@result,paste0("./ResultsAllVariants/GO/GO.houseKeeping_gene.txt"),sep='\t',quote=F,col.names = T,row.names = F)
pdf(paste0("./ResultsAllVariants/GO/GO.houseKeeping_gene.pdf"),width = 10,height = 6)
print(dotplot(en_go))
dev.off()

en_go=enricher(ExtremelyVaried,TERM2GENE = term2gene,TERM2NAME = term2dec,pvalueCutoff = 0.05,pAdjustMethod = "BH")
write.table(en_go@result,paste0("./ResultsAllVariants/GO/GO.ExtremelyVaried_gene.txt"),sep='\t',quote=F,col.names = T,row.names = F)
pdf(paste0("./ResultsAllVariants/GO/GO.ExtremelyVaried_gene.pdf"),width = 10,height = 6)
print(dotplot(en_go))
dev.off()

en_go=enricher(varied,TERM2GENE = term2gene,TERM2NAME = term2dec,pvalueCutoff = 0.05,pAdjustMethod = "BH")
write.table(en_go@result,paste0("./ResultsAllVariants/GO/GO.varied_gene.txt"),sep='\t',quote=F,col.names = T,row.names = F)
pdf(paste0("./ResultsAllVariants/GO/GO.varied_gene.pdf"),width = 10,height = 6)
print(dotplot(en_go))
dev.off()

# varied between populations
setwd("E:/QiGA/document/CottonLab/FH")
library(clusterProfiler)
library(enrichplot)
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

houseKeeping = fread("./HouseKeeping_gene.list",header=F)
go_gene = read.table("E:/QiGA/document/CottonLab/Cotton_GO/TM-1_final.agriGO.txt")
go_desc = read.csv("E:/QiGA/document/CottonLab/Cotton_GO/GO_desc_all.csv")
go_desc$nameform = substr(go_desc$name, 1, 100)
term2gene = data.frame(TERM = go_gene$V2, gene = go_gene$V1)
term2dec = data.frame(TERM = go_desc$Goid, NAME = go_desc$nameform)

source.unique = unique(group$`Source Region`)
for(i in 1:length(source.unique)){
  sourceR = source.unique[i]
  out_source = sub(" ","_",sourceR,fixed = T)
  dt = fread(paste0("./ResultsAllVariants//FHAP_", out_source, ".GeneSTAT.txt"),header=T,data.table = F,stringsAsFactors = F)
  
  
  houseKeeping.Gene = union(houseKeeping$V1, dt$Gene[which(dt$EH==0)])
  write.table(houseKeeping.Gene,paste0("./ResultsAllVariants/GO/GO.", out_source, ".houseKeeping_Gene.list"),sep='\t',quote=F,col.names = F,row.names = F)
  if(length(houseKeeping.Gene)>10){
    en_go=enricher(houseKeeping.Gene,TERM2GENE = term2gene,TERM2NAME = term2dec,pvalueCutoff = 0.05,pAdjustMethod = "BH")
    go.dt = en_go@result
    if(length(which(go.dt$p.adjust<=0.05))>0){
      write.table(go.dt,paste0("./ResultsAllVariants/GO/GO.", out_source, ".houseKeeping_gene.txt"),sep='\t',quote=F,col.names = T,row.names = F)
      pdf(paste0("./ResultsAllVariants/GO/GO.", out_source, ".houseKeeping_gene.pdf"),width = 10,height = 6)
      print(dotplot(en_go))
      dev.off()
    }
  }
  
  lowDensity.Gene = dt$Gene[which(dt$EH>0 & dt$EH<=0.05)]
  write.table(lowDensity.Gene,paste0("./ResultsAllVariants/GO/GO.", out_source, ".lowDensity_Gene.list"),sep='\t',quote=F,col.names = F,row.names = F)
  if(length(lowDensity.Gene)>10){
    en_go=enricher(lowDensity.Gene,TERM2GENE = term2gene,TERM2NAME = term2dec,pvalueCutoff = 0.05,pAdjustMethod = "BH")
    go.dt = en_go@result
    if(length(which(go.dt$p.adjust<=0.05))>0){
      write.table(go.dt,paste0("./ResultsAllVariants/GO/GO.", out_source, ".lowDensity_Gene.txt"),sep='\t',quote=F,col.names = T,row.names = F)
      pdf(paste0("./ResultsAllVariants/GO/GO.", out_source, ".lowDensity_Gene.pdf"),width = 10,height = 6)
      print(dotplot(en_go))
      dev.off()
    }
  }
  
  
  moderateDensity.Gene = dt$Gene[which(dt$EH>0.05 & dt$EH<=0.3)]
  write.table(moderateDensity.Gene,paste0("./ResultsAllVariants/GO/GO.", out_source, ".moderateDensity_Gene.list"),sep='\t',quote=F,col.names = F,row.names = F)
  if(length(moderateDensity.Gene)>10){
    en_go=enricher(moderateDensity.Gene,TERM2GENE = term2gene,TERM2NAME = term2dec,pvalueCutoff = 0.05,pAdjustMethod = "BH")
    go.dt = en_go@result
    if(length(which(go.dt$p.adjust<=0.05))>0){
      write.table(go.dt,paste0("./ResultsAllVariants/GO/GO.", out_source, ".moderateDensity_Gene.txt"),sep='\t',quote=F,col.names = T,row.names = F)
      pdf(paste0("./ResultsAllVariants/GO/GO.", out_source, ".moderateDensity_Gene.pdf"),width = 10,height = 6)
      print(dotplot(en_go))
      dev.off()
    }
  }
  
  highDensity.Gene = dt$Gene[which(dt$EH>0.3)]
  write.table(highDensity.Gene,paste0("./ResultsAllVariants/GO/GO.", out_source, ".highDensity_Gene.list"),sep='\t',quote=F,col.names = F,row.names = F)
  if(length(highDensity.Gene)>10){
    en_go=enricher(highDensity.Gene,TERM2GENE = term2gene,TERM2NAME = term2dec,pvalueCutoff = 0.05,pAdjustMethod = "BH")
    go.dt = en_go@result
    if(length(which(go.dt$p.adjust<=0.05))>0){
      write.table(go.dt,paste0("./ResultsAllVariants/GO/GO.", out_source, ".highDensity_Gene.txt"),sep='\t',quote=F,col.names = T,row.names = F)
      pdf(paste0("./ResultsAllVariants/GO/GO.", out_source, ".highDensity_Gene.pdf"),width = 10,height = 6)
      print(dotplot(en_go))
      dev.off()
    }
  }
}

# enrichment for varied genes
library(clusterProfiler)
library(enrichplot)

earlyAcc = c("South America","North America","Central America","SC (CHN)")
improvedAcc = c("USA","YER (CHN)","YZR (CHN)","NWC (CHN)")

DT = data.frame()
for(i in 1:length(c(earlyAcc,improvedAcc))){
  name = sub(" ","_",c(earlyAcc,improvedAcc)[i],fixed=T)
  dt = fread(paste0("./ResultsAllVariants/FHAP_", name, ".GeneSTAT.txt"),header=T,data.table = F,stringsAsFactors = F)
  
  dt = dt[,c(1,4)]
  colnames(dt)[2] = name
  if(i == 1){
    DT = dt
  } else {
    DT = merge(DT,dt,"Gene",all = T)
  }
}
DT = na_replace(DT,0)

test.sig = apply(DT,1,function(x){
  x = as.numeric(x[-1])
  #t.test(x[c(1:4)],x[c(5:8)],paired=F)$p.value
  wilcox.test(x[c(1:4)],x[c(5:8)],paired=F,exact = F)$p.value
})

dev = apply(DT,1,function(x){
  x = as.numeric(x[-1])
  mean(x[c(1:4)])-mean(x[c(5:8)])
  #wilcox.test(x[c(1:4)],x[c(5:8)],paired=F,exact = F)$p.value
})

DT$Sig = test.sig
DT$Sig = na_replace(DT$Sig,1)

DT$DEV = dev
write.table(DT,"./ResultsAllVariants/GO/Varied_Eh_genes.summary.txt",quote=F,col.names = T,row.names = F,sep='\t')
genes1 = DT$Gene[which(DT$DEV >= 0.1 & DT$DEV<=0.3)]
genes2 = DT$Gene[which(DT$DEV >= 0.3)]
en_go = enricher(c(genes1,genes2),TERM2GENE = term2gene,TERM2NAME = term2dec,pvalueCutoff = 0.05,pAdjustMethod = "BH")
en_go1 = enricher(genes1,TERM2GENE = term2gene,TERM2NAME = term2dec,pvalueCutoff = 0.05,pAdjustMethod = "BH")
en_go2 = enricher(genes2,TERM2GENE = term2gene,TERM2NAME = term2dec,pvalueCutoff = 0.05,pAdjustMethod = "BH")

write.table(en_go@result,"./ResultsAllVariants/GO/Varied_Eh_genes.GO.txt",quote=F,col.names = T,row.names = F,sep = '\t')
write.table(en_go1@result,"./ResultsAllVariants/GO/Varied_moderate_Eh_genes.GO.txt",quote=F,col.names = T,row.names = F,sep = '\t')
write.table(en_go2@result,"./ResultsAllVariants/GO/Varied_rapid_Eh_genes.GO.txt",quote=F,col.names = T,row.names = F,sep = '\t')

# plot new GO
SSDT = read.table("./ResultsAllVariants/GO/Varied_Eh_genes.summary.txt", header = T, sep = '\t')
GO1 = read.table("./ResultsAllVariants/GO/Varied_moderate_Eh_genes.GO.txt", header = T, sep = '\t')
GO2 = read.table("./ResultsAllVariants/GO/Varied_rapid_Eh_genes.GO.txt", header = T, sep = '\t',quote="")

GO1 = GO1[which(GO1$p.adjust<=0.01),]
GO2 = GO2[which(GO1$p.adjust<=0.01),]

## construct chord DT
library(circlize)
df = data.frame()
for(i in 1:nrow(GO1)){
  geneset = GO1$geneID[i]
  field = GO1$Description[i]
  
  genes = strsplit(geneset,"/",fixed=T)[[1]]
  temp = data.frame(Gene = genes, Function = field)
  df = rbind(df,temp)
}
for(i in 1:nrow(GO2)){
  geneset = GO2$geneID[i]
  field = GO2$Description[i]
  
  genes = strsplit(geneset,"/",fixed=T)[[1]]
  temp = data.frame(Gene = genes, Function = field)
  df = rbind(df,temp)
}
dupCheck = paste0(df$Gene,"_",df$Function)
dupIndex = which(duplicated(dupCheck))

if(length(dupIndex)>0){
  df = df[-dupIndex,]
}

funcCheck = table(df$Function)
funcRemove = names(funcCheck)[which(funcCheck<10)]
df$Function[which(!is.na(match(df$Function,funcRemove)))] = 'Others'

chrod = merge(df,SSDT[,c(1,11)],"Gene",all.x = T)
chrod$GeneSymbol = "G"
chrod$GeneSymbol[which(chrod$DEV>0.1 & chrod$DEV<=0.3)] = "Moderate-varied"
chrod$GeneSymbol[which(chrod$DEV>0.3)] = "Rapid-varied"
chrod$Function[which(chrod$Function=='ATPase activity, coupled to transmembrane movement of substances')] = "ATPase activity"
chrod$Index = paste0(chrod$GeneSymbol,"_",chrod$Function)
chrod = as.data.frame(table(chrod$Index))
chrod$Var1 = as.character(chrod$Var1)
chrod$From = sapply(strsplit(chrod$Var1,"_",fixed=T),'[[',1)
chrod$To = sapply(strsplit(chrod$Var1,"_",fixed=T),'[[',2)
chrod = chrod[,c("From","To","Freq")]

gridCol = pal_igv()(17)
names(gridCol) = c("Moderate-varied","Rapid-varied",
                   "carbohydrate binding",
                   "polysaccharide binding",
                   "microtubule binding",
                   "microtubule motor activity",
                   "microtubule-based movement",
                   "ADP binding",
                   "ATPase activity",
                   "ubiquitin-protein transferase activity",
                   "protein ubiquitination",
                   "cytokinin metabolic process",
                   "cytokinin dehydrogenase activity",
                   "cytokinesis by cell plate formation",
                   "recognition of pollen",
                   "intracellular signal transduction",
                   "Others")

pdf("./ResultsAllVariants/GO/Varied_Eh_genes.GO.plot2.pdf",width = 6,height = 6)
par(mar=c(0,0,0,0))
circos.par(start.degree = -95,
           gap.after=c(rep(1,length(unique(chrod$GeneSymbol))-1),10,rep(1,length(unique(chrod$Function))-1),10))
chordDiagram(chrod,big.gap = 20,grid.col = gridCol,
             directional = 1,
             link.lwd = 0.1,
             link.lty = 1,    
             link.border = 0,
             order = c("Moderate-varied","Rapid-varied",
                       "carbohydrate binding",
                       "polysaccharide binding",
                       "microtubule binding",
                       "microtubule motor activity",
                       "microtubule-based movement",
                       "ADP binding",
                       "ATPase activity",
                       "ubiquitin-protein transferase activity",
                       "protein ubiquitination",
                       "cytokinin metabolic process",
                       "cytokinin dehydrogenase activity",
                       "cytokinesis by cell plate formation",
                       "recognition of pollen",
                       "intracellular signal transduction",
                       "Others"),
             transparency = 0.5,
             annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.3,0.5))
circos.clear()
dev.off()
