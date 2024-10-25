# gene flow for dicMigrate
# convert to genepop format
setwd("E:/QiGA/document/CottonLab/FH")
library(data.table)
GT = fread("./ResultsAllVariants/FHAP_whole.type.2.txt", sep = '\t', data.table = F, stringsAsFactors = F)
pop_info = read.table("E:/QiGA/document/CottonLab/XLD/population_info/Detailed_Information.txt",header=T,sep='\t',quote = "",comment.char = "@",encoding = 'utf_8')
pop_info = pop_info[,c("VCF_ID","Data.interval","Zone2")]
pop_info$VCF_ID = paste0("0_",pop_info$VCF_ID)
colnames(pop_info) = c("label","Cultivated Date","Source Region")
GTM = as.matrix(GT[,-1])
rownames(GTM) = GT[,1]

total_haplotypes = apply(GTM,1,function(x){
  length(unique(x))
})
indexHighVar = which(total_haplotypes>300)
if(length(indexHighVar)>0){
  GTM = GTM[-indexHighVar,]
}

# write genepop header
write("3KCotton Project Function Haplotype GenePop",file = "FHAP_whole.type.3.genepop.txt",append=F)

pops = c("USA","NC (CHN)","SWC (CHN)","YZR (CHN)","FSU","YER (CHN)", "Africa",
         "South America","North America","Europe","Asia","Oceania","SC (CHN)","NWC (CHN)","GB","Central America")
geneNames = rownames(GTM)
indiNames = colnames(GTM)
write.table(geneNames,file = "FHAP_whole.type.3.genepop.txt",append=T,col.names = F,row.names = F,quote=F)
for(i in 1:length(pops)){
  pop = pops[i]
  pop1 = sub(" ","_",pop,fixed=T)
  
  targetSamples = pop_info$label[which(pop_info$`Source Region`==pop)]
  targetSamples = intersect(indiNames,targetSamples)
  if(length(targetSamples)==0){
    next
  }
  write(paste0("POP"),file = "FHAP_whole.type.3.genepop.txt",append=T)
  targetMat = t(GTM[,targetSamples])
  targetMat = apply(targetMat,2,function(x){
    paste0(sprintf("%03d",x),sprintf("%03d",x))
  })
  rownames(targetMat) = targetSamples
  targetMat = cbind(paste0(pop1,","),targetMat)
  write.table(targetMat,"E:/QiGA/document/CottonLab/FH/ResultsAllVariants/ResFHAP_whole.type.3.thres300.genepop.txt",append=T,row.names = F,col.names = F,sep=' ',quote=F)
}

library(diveRsity)
infile = "E:/QiGA/document/CottonLab/FH/ResultsAllVariants/FHAP_whole.type.3.genepop.txt"
dat <- rgp(infile)

res = 
  divMigrate(infile = "E:/QiGA/document/CottonLab/FH/ResultsAllVariants/FHAP_whole.type.3.genepop.txt",
           outfile = "divMargin",
           boots = 100,
           stat = 'all',
           para = F,
           plot_network = T)
save.image("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GeneFlow/DivMigration.res.RData")
dev.off()


infile = "E:/QiGA/document/CottonLab/FH/FHAP_whole.type.3.genepop.txt"
outfile = "divMargin"
boots = 10
stat = 'all'
para = T
plot_network = T
filter_threshold = 0
plot_col = "darkblue"

pops = c("USA","NC (CHN)","SWC (CHN)","YZR (CHN)","FSU","YER (CHN)", "Africa",
         "South America","North America","Europe","Asia","Oceania","SC (CHN)","NWC (CHN)","GB","Central America")
gRelPlt <- res[["gRelMig"]]
filter_threshold = quantile(gRelPlt,na.rm=T,probs = 0.9)
gRelPlt[gRelPlt < filter_threshold] <- 0
qgraph::qgraph(gRelPlt, nodeNames = sapply(pops, 
                                           "[", 1), legend = TRUE, posCol = plot_col, edge.labels = F, 
               mar = c(2, 2, 5, 5), curve = 2.5)
title(paste("\n Relative migration network \n (Filter threshold = ", 
            filter_threshold, "; Gst)", sep = ""))