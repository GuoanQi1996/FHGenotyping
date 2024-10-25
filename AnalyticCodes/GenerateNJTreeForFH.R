# generate tree
setwd("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/")
library(data.table)
library(distances)
library(ape)
GT = fread("FHAP_whole.type.2.txt", sep = '\t', data.table = F, stringsAsFactors = F)
rownames(GT) = GT$Gene
GT = as.data.frame(t(GT[,-1]))

# filter high varied object/variants
numeric.var.vec = apply(GT,2,var,na.rm=T)
quant.numeric.var = quantile(numeric.var.vec,probs = seq(0,1,0.01))

factor.var.vec = apply(GT,2,function(x){
  frq = table(x)
  base = names(frq)[which(frq==max(frq))]
  length(which(x!=base))
})
quant.factor.var = quantile(factor.var.vec,probs = seq(0,1,0.01))

# plot(density(factor.var.vec[which(numeric.var.vec<2)]))
# quantile(factor.var.vec[which(numeric.var.vec<2)],probs = seq(0,1,0.01))

GT.V = GT[,which(numeric.var.vec>=2)]

# Euclidean distances
# dist = distances(GT.V,id_variable = rownames(GT.V))
# distMat = distance_matrix(dist)
# nj.tree = nj(distMat)

# dist.gene for character matrix
GT.V = apply(GT.V,2,as.character)
disMat = dist.gene(GT.V, method = 'pairwise')
nj.tree3 = nj(disMat)
nj.tree3[["tip.label"]] = rownames(GT)


#clustMat = hclust(distMat, method = 'complete')
#tree = as.phylo(clustMat)

# 输出newick格式文件
write.tree(nj.tree3,file = "./Phy/FH.Var2.Factor.NJ.tree")