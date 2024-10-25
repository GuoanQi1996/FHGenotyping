setwd("E:/QiGA/document/CottonLab/FH")
# process big annotation 
library(data.table)
dt = fread("data/FH_recode.ANN.vcf",select = c(3,8),sep='\t',header = T,stringsAsFactors = F)

info = strsplit(sapply(strsplit(dt$INFO,";",fixed=T),'[[',2),"|",fixed=T)
dt$Gene = sapply(info,'[[',4)
dt$Type = sapply(info,'[[',2)
rm(info)

dt = dt[,c(1,3,4)]

gc()
colnames(dt)[1] = "SNP"
write.table(dt,"FH.annotation.txt",quote=F,col.names = T,row.names = F)

# add gene information
dt = fread("FH.annotation.txt",header=T,stringsAsFactors = F)
gff = fread("../Cotton_GO/TM-1_V2.1.gene.CDS.gff",header=F,stringsAsFactors = F,data.table = F)
allGene = gff[which(gff$V3=='gene'),]
allGene$V9 = sub("ID=","",sub(";","",allGene$V9,fixed = T),fixed = T)
allGene = data.frame(Gene = allGene$V9)
# gene body SNP dt
dt = dt[-which(dt$Type=='intergenic_region' | dt$Type=='downstream_gene_variant' | dt$Type=='upstream_gene_variant'),]
snpTypePerGene = aggregate(dt$Type,list(dt$Gene),function(x){
  table(x)
})
TypePerGene = data.table(Gene = snpTypePerGene$Group.1, VariantsCount=sapply(snpTypePerGene$x,sum))
TypePerGene = merge(allGene,TypePerGene,"Gene",all.x=T)
TypePerGene$VariantsCount[which(is.na(TypePerGene$VariantsCount))] = 0
write.table(TypePerGene,"GeneBodyVariantsSumPerGene.txt",quote=F,col.names = T,row.names = F)
write.table(dt,"GeneBodyVariantsPerGene.txt",quote=F,col.names = T,row.names = F)

# gene body synonymous and missense 
dt = dt[-which(grepl("intron_variant",dt$Type) | grepl("5_prime",dt$Type)),]
as.data.frame(table(dt$Type)) -> test
write.table(dt,"Synonymous_nonsynonymousPerGene.txt",quote=F,col.names = F,row.names = F)
write.table(dt$SNP,"Synonymous_nonsynonymous.variant",quote=F,col.names = F,row.names = F)

# gene body high-effect SNP
dt = dt[-which(grepl("synonymous_variant",dt$Type) | grepl("intron_variant",dt$Type) | grepl("5_prime",dt$Type)),]

# house-keeping gene 
houseKeeping = setdiff(allGene$Gene,unique(dt$Gene))
write.table(houseKeeping,"HouseKeeping_gene.list",quote=F,col.names = F,row.names = F)


missenseTypePerGene = aggregate(dt$Type,list(dt$Gene),function(x){
  table(x)
})
missenseTypePerGene = data.table(Gene = missenseTypePerGene$Group.1, MissenseCount=sapply(missenseTypePerGene$x,sum))
missenseTypePerGene = merge(allGene,missenseTypePerGene,"Gene",all.x=T)
missenseTypePerGene$MissenseCount[which(is.na(missenseTypePerGene$MissenseCount))] = 0
write.table(missenseTypePerGene,"missenseSumPerGene.txt",quote=F,col.names = T,row.names = F)
write.table(dt,"missensePerGene.txt",quote=F,col.names = T,row.names = F)
write.table(dt$SNP,"missenseVariants.txt",quote=F,col.names = T,row.names = F)
