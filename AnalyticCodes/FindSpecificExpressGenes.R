library(data.table)
exp = fread("E:/QiGA/document/CottonLab/Cotton_GO/Expression_TM1_all_organisms_New/Expression_Level_FPKM.csv",header=T,data.table = F)
DT.wide = exp

expr.universal.counts = apply(DT.wide[,-1],1,function(x){
  length(which(x>=1))
})
houseKeeping_gene = DT.wide$Gene[which(expr.universal.counts==ncol(DT.wide)-1)]

Fiber_expression_index = grep("Fiber",colnames(DT.wide))
expr.Fiber.counts = apply(DT.wide[,Fiber_expression_index],1,function(x){
  length(which(x>=1))
})
Fiber_expression_gene = DT.wide$Gene[which(expr.Fiber.counts>0)]

Ovule_expression_index = grep("Ovule",colnames(DT.wide))
expr.Ovule.counts = apply(DT.wide[,Ovule_expression_index],1,function(x){
  length(which(x>=1))
})
Ovule_expression_gene = DT.wide$Gene[which(expr.Ovule.counts>0)]

cotyledon_expression_index = grep("cotyledon",colnames(DT.wide))
expr.cotyledon.counts = apply(DT.wide[,cotyledon_expression_index],1,function(x){
  length(which(x>=1))
})
cotyledon_expression_gene = DT.wide$Gene[which(expr.cotyledon.counts>0)]

root_expression_index = grep("root",colnames(DT.wide))
expr.root.counts = apply(DT.wide[,root_expression_index],1,function(x){
  length(which(x>=1))
})
root_expression_gene = DT.wide$Gene[which(expr.root.counts>0)]

seed_expression_index = grep("seed",colnames(DT.wide))
expr.seed.counts = apply(DT.wide[,seed_expression_index],1,function(x){
  length(which(x>=1))
})
seed_expression_gene = DT.wide$Gene[which(expr.seed.counts>0)]

calycle_expression_index = grep("calycle",colnames(DT.wide))
calycle_expression_gene = DT.wide$Gene[which(DT.wide[,calycle_expression_index]>=1)]

stamen_expression_index = grep("stamen",colnames(DT.wide))
stamen_expression_gene = DT.wide$Gene[which(DT.wide[,stamen_expression_index]>=1)]

stem_expression_index = grep("stem",colnames(DT.wide))
stem_expression_gene = DT.wide$Gene[which(DT.wide[,stem_expression_index]>=1)]

leaf_expression_index = grep("leaf",colnames(DT.wide))
leaf_expression_gene = DT.wide$Gene[which(DT.wide[,leaf_expression_index]>=1)]

petal_expression_index = grep("petal",colnames(DT.wide))
petal_expression_gene = DT.wide$Gene[which(DT.wide[,petal_expression_index]>=1)]

pistil_expression_index = grep("pistil",colnames(DT.wide))
pistil_expression_gene = DT.wide$Gene[which(DT.wide[,pistil_expression_index]>=1)]

torus_expression_index = grep("torus",colnames(DT.wide))
torus_expression_gene = DT.wide$Gene[which(DT.wide[,torus_expression_index]>=1)]

# find organism-specific genes 
organisms = unique(sapply(strsplit(colnames(DT.wide)[-1],'_',fixed=T),'[[',1))

for(i in 1:length(organisms)){
  organism = organisms[i]
  
  others = setdiff(organisms,organism)
  myExpressionGene = get(paste0(organism,"_expression_gene"))
  myExpressionGene = setdiff(myExpressionGene,houseKeeping_gene)
  for(j in 1:length(others)){
    myExpressionGene = setdiff(myExpressionGene,get(paste0(others[j],"_expression_gene")))
  }
  assign(paste0(organism,"_specific_gene"),myExpressionGene)
  write.table(myExpressionGene,paste0("E:/QiGA/document/CottonLab/Cotton_GO/Expression_TM1_all_organisms_New/Specific_express.",organism,".genelist"),quote=F,col.names = F,row.names = F)
}

for(i in 1:(length(organisms)-1)){
  organism1 = organisms[i]
  for(j in (i+1):length(organisms)){
    organism2 = organisms[j]
    others = setdiff(organisms,organism1)
    others = setdiff(others,organism2)
    myExpressionGene = get(paste0(organism,"_expression_gene"))
    myExpressionGene = setdiff(myExpressionGene,houseKeeping_gene)
    for(j in 1:length(others)){
      myExpressionGene = setdiff(myExpressionGene,get(paste0(others[j],"_expression_gene")))
    }
    assign(paste0(organism1,"_",organism2,"_specific_gene"),myExpressionGene)
    write.table(myExpressionGene,paste0("E:/QiGA/document/CottonLab/Cotton_GO/Expression_TM1_all_organisms_New/Specific_express.",organism1,"_", organism2, ".genelist"),quote=F,col.names = F,row.names = F)
  }
}

SCG = read.table("./HouseKeeping_gene.list",header=F)
organisms = unique(sapply(strsplit(colnames(DT.wide)[-1],'_',fixed=T),'[[',1))
for(i in 1:length(organisms)){
  organism = organisms[i]
  write.table(intersect(get(paste0(organism,"_specific_gene")),SCG$V1),paste0("./ResultsAllVariants/GO/SCG.",organism,"_specific_gene.txt"),quote=F,col.names = F,row.names = F,sep='\t')
  cat(paste0("SCG have ",length(intersect(get(paste0(organism,"_specific_gene")),SCG$V1)),"/", length(get(paste0(organism,"_specific_gene"))), ", proportion ", 100*length(intersect(get(paste0(organism,"_specific_gene")),SCG$V1))/length(get(paste0(organism,"_specific_gene")))," with ",organism,".\n"))
}
slient_genes = DT.wide$Gene[which(expr.universal.counts==0)]
SCG_slient = intersect(SCG$V1,slient_genes)
write.table(SCG_slient,paste0("./ResultsAllVariants/GO/SCG.Silent.txt"),quote=F,col.names = F,row.names = F,sep='\t')


SCG_HKG1 = intersect(houseKeeping_gene,SCG$V1)
write.table(SCG_HKG1,paste0("./ResultsAllVariants/GO/SCG.HKG.txt"),quote=F,col.names = F,row.names = F,sep='\t')

SCG_Others = setdiff(SCG$V1,c(houseKeeping_gene,Ovule_specific_gene,Fiber_specific_gene,
                              calycle_specific_gene,stamen_specific_gene,
                              leaf_specific_gene,petal_specific_gene,pistil_specific_gene,
                              root_specific_gene,cotyledon_specific_gene,seed_specific_gene,
                              stem_specific_gene,torus_specific_gene,slient_genes))
write.table(SCG_Others,paste0("./ResultsAllVariants/GO/SCG.Others.txt"),quote=F,col.names = F,row.names = F,sep='\t')

go_gene = read.table("E:/QiGA/document/CottonLab/Cotton_GO/TM-1_final.agriGO.txt")
go_desc = read.csv("E:/QiGA/document/CottonLab/Cotton_GO/GO_desc_all.csv")
go_desc$nameform = substr(go_desc$name, 1, 100)
term2gene = data.frame(TERM = go_gene$V2, gene = go_gene$V1)
term2dec = data.frame(TERM = go_desc$Goid, NAME = go_desc$nameform)
library(clusterProfiler)
en_go=enricher(SCG_Others,TERM2GENE = term2gene,TERM2NAME = term2dec,pvalueCutoff = 0.05,pAdjustMethod = "BH")
write.table(en_go@result,"./ResultsAllVariants/GO/GO.All.SCG_except_organism_specific_and_silent.txt",quote=F,sep ='\t',col.names = T,row.names = F)
pdf("./ResultsAllVariants/GO/SCG_Others.pdf",width = 10,height = 6)
print(dotplot(en_go))
dev.off()

x1 = c(5027,3542,2602,669)
colors = c("grey80", pal_d3('category20',alpha = 0.9)(3))
pdf("Annotation.SCG.pie.pdf",width = 3,height = 3)
par(mar = c(0, 0, 0, 0))
pie(x1,labels = "",init.angle = -80,col = colors)
dev.off()

HVG = read.table("./ResultsAllVariants/FHAP_whole.GeneSTAT.txt",header=T)
HVG = HVG[which(HVG$EH>=0.3),] # original 0.3
write.table(HVG,"./ResultsAllVariants/GO/HVG.data.frame.txt",quote=F,col.names = T,row.names = F,sep ='\t')
HVG = HVG[,c(1:3)]
classification = readxl::read_excel("Gene_symbol_classification.xlsx",sheet = 2)
anno.gene = read.table("E:/QiGA/document/CottonLab/Cotton_GO/TM-1_v2.1_final.annotation.txt",header=T,sep = '\t',quote="")
colnames(anno.gene)[1] = 'Gene'
#HVG = merge(HVG,classification,'Gene',all.x = T)
HVG = merge(HVG,anno.gene,'Gene',all.x = T)
write.table(HVG,"./HVG.ALL.Annotation.txt",quote=F,col.names = T,row.names = F,sep = '\t')

d1 = as.data.frame(sort(table(HVG$ArabDesc)/sum(table(HVG$ArabDesc)),decreasing = T))
d2 = as.data.frame(sort(table(HVG$ArabDesc),decreasing = T))
d = cbind(d1,d2$Freq)
colnames(d) = c("Gene description","Numbers (in HVGs)","Proportion")
write.table(d,"./HVG.ALL.Annotation.Summary.txt",quote=F,col.names = T,row.names = F,sep = '\t')

x2 = c(69,225,309)
color2 = c("#8C564BE5","#D62728E5","grey80")
pdf("Annotation.HVG.pie.pdf",width = 3,height = 3)
par(mar = c(0, 0, 0, 0))
pie(x2,labels = "",init.angle = -0,col = color2)
dev.off()