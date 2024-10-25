dist2df <- function(inDist) {
  if (class(inDist) != "dist") stop("wrong input type")
  A <- attr(inDist, "Size")
  B <- if (is.null(attr(inDist, "Labels"))) sequence(A) else attr(inDist, "Labels")
  if (isTRUE(attr(inDist, "Diag"))) attr(inDist, "Diag") <- FALSE
  if (isTRUE(attr(inDist, "Upper"))) attr(inDist, "Upper") <- FALSE
  data.frame(
    row = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
    col = rep(B[-length(B)], (length(B)-1):1),
    value = as.vector(inDist))
}
setwd("E:/QiGA/document/CottonLab/FH")
library(data.table)
library(ape)
GT = fread("./ResultsAllVariants/FHAP_whole.type.2.txt", sep = '\t', data.table = F, stringsAsFactors = F)
allSample = colnames(GT)[-1]

GT.C = t(GT[,-1])
GT.C = apply(GT.C,2,'as.character')
rownames(GT.C) = allSample
GT.dis = dist.gene(GT.C,method = 'per')
GT.dis.DT = dist2df(GT.dis)
GT.dis.DT$row = sub("0_","",GT.dis.DT$row,fixed=T)
GT.dis.DT$col = sub("0_","",GT.dis.DT$col,fixed=T)
GT.dis.DT$value = 1-GT.dis.DT$value
colnames(GT.dis.DT) = c("IID1","IID2","FHKinship")

fwrite(GT.dis.DT,"FH.GeneticDis.txt.gz",sep = '\t', quote=F,col.names=T,row.names = T)

# compare SNPs-based and FH-based kinship
FHKinship = fread("FH.GeneticDis.txt.gz",header=T,sep = '\t',data.table=F,stringsAsFactors = F)
SNPKinship = fread("FH_include_SNP.king.kin0.gz",header=T,data.table=F,stringsAsFactors = F)
SNPKinship = SNPKinship[,c(2,4,8)]
FHKinship = FHKinship[,c(2,3,4)]

SNPKinship = SNPKinship[order(SNPKinship$IID1,SNPKinship$IID2,decreasing = F),]
FHKinship = FHKinship[order(FHKinship$IID1,FHKinship$IID2,decreasing = F),]
which(SNPKinship$IID1!=FHKinship$IID1)
which(SNPKinship$IID2!=FHKinship$IID2)

cor(SNPKinship$KINSHIP,FHKinship$FHKinship) # 0.7093
plot(scale(SNPKinship$KINSHIP),scale(FHKinship$FHKinship))

# prepare data for th regression
fill_missing_symmetric <- function(matrix) {
  n <- nrow(matrix)
  if(n != ncol(matrix)) {
    stop("The matrix is not square.")
  }
  
  for (i in 2:n) {
    for (j in (i-1):1) {
      if(is.na(matrix[i, j]) && !is.na(matrix[j, i])) {
        matrix[i, j] <- matrix[j, i]
      } else if (!is.na(matrix[i, j]) && is.na(matrix[j, i])) {
        matrix[j, i] <- matrix[i, j]
      }
    }
  }
  return(matrix)
}
setwd("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM")
library(data.table)
library(lme4qtl)
library(lme4)
library(doSNOW)

cl = makeCluster(16)
registerDoSNOW(cl)

GT = fread("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/FHAP_whole.type.2.txt", sep = '\t', data.table = F, stringsAsFactors = F)
colnames(GT) = sub("0_","",colnames(GT))
GT1 = t(GT[,-1])
colnames(GT1) = GT$Gene
GT1 = data.frame(ID = rownames(GT1),GT1)

FHKinship = fread("E:/QiGA/document/CottonLab/FH/FH.GeneticDis.txt",header=T,sep = '\t',data.table=F,stringsAsFactors = F)
#SNPKinship = fread("E:/QiGA/document/CottonLab/FH/FH_include_SNP.king.kin0.gz",header=T,data.table=F,stringsAsFactors = F)
#SNPKinship = SNPKinship[,c(2,4,8)]
FHKinship = FHKinship[,c(2,3,4)]
FHKinship = rbind(FHKinship,data.frame(IID1 = FHKinship$IID2,IID2 = FHKinship$IID1, FHKinship = FHKinship$FHKinship))

FH.Mat = reshape2::dcast(data = FHKinship, IID1~IID2, value.var = 'FHKinship')
#SNP.Mat = reshape2::dcast(data = SNPKinship, IID1~IID2, value.var = 'KINSHIP')
rownames(FH.Mat) = FH.Mat$IID1
#rownames(SNP.Mat) = SNP.Mat$IID1
FH.Mat = as.matrix(FH.Mat[,-1])
diag(FH.Mat) = 1
#write.table(FH.Mat,"E:/QiGA/document/CottonLab/FH/FH.GeneticDis.Mat.txt",quote=F,col.names = T,row.names = F,sep='\t')
#SNP.Mat = SNP.Mat[,-1]

pheList = list.files("E:/QiGA/document/CottonLab/XLD/Phe/EMMAX/Seprated/Phe/","EMMAX.*.txt")
pheList = c(paste0("EMMAX.HZG486.",c("PH_AVE","PH_2017_SHZ","PH_2016_SHZ","FBN_2017_SHZ","FBN_SHZ_AVE","FBN_AVE","BN_2017_SHZ","BN_AVE"),".txt"))

#PHE = sub("EMMAX_","",sapply(strsplit(pheList,".",fixed=T),'[[',1))
phe_field1 = sapply(strsplit(pheList,".",fixed=T),'[[',2)
phe_field2 = sapply(strsplit(pheList,".",fixed=T),'[[',3)
PHE = paste0(phe_field1,"_",phe_field2)

PCs = read.table("E:/QiGA/document/CottonLab/XLD/Phe/PC.Cov.txt",header=F,sep='\t',comment.char = "")
PCs = PCs[,c(2:8)]
colnames(PCs) = c("ID","Intecept",paste0("PC",1:5))

for(i in 1:length(PHE)){
  pheName = PHE[i]
  if(file.exists(paste0("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Seprate/Phe_",pheName,".GLMM.csv"))){
    cat(paste0("Trait ",pheName," is already done GLMM, skip...\n"))
    next
  }
  cat(paste0("\nProceeding trait ",pheName,". "))
  #phe = read.table(paste0("E:/QiGA/document/CottonLab/XLD/Phe/EMMAX_", pheName, ".txt"),header = F,sep='\t',comment.char = "")
  phe = read.table(paste0("E:/QiGA/document/CottonLab/XLD/Phe/EMMAX/Seprated/Phe/EMMAX.", phe_field1[i], ".", phe_field2[i], ".txt"),header = F,sep='\t',comment.char = "")
  phe = phe[,c(2,3)]
  colnames(phe) = c("ID","PHE")
  phe_cov = merge(phe,PCs,"ID")
  phe_cov = na.omit(phe_cov)
  phe_cov$ID = as.factor(phe_cov$ID)
  sample_include = as.character(phe_cov$ID)
  
  FH.Kin = FH.Mat[sample_include,sample_include]
  #SNP.Kin = SNP.Mat[sample_include,sample_include]
  
  data = merge(phe_cov,GT1,by='ID')
  colVar = apply(data[,-1],2,var)
  var0 = as.numeric(which(colVar==0))+1
  data = data[,-var0]
  
  FHKModel = relmatLmer(PHE ~ PC1+PC2+PC3+PC4+PC5+(1|ID), 
                        data, 
                        relmat = list(ID = FH.Kin),
                        REML=F,
                        control = lmerControl(optimizer = "nloptwrap",
                                              optCtrl=list(xtol_abs=1e-6, ftol_abs=1e-6),
                                              calc.derivs = F))
  KinEff = ranef(FHKModel)$ID
  KinEff$ID = rownames(KinEff)
  colnames(KinEff)[1] = 'KinEff'
  data = merge(data,KinEff,"ID")
  data$PHENK = data$PHE-data$KinEff
  width = ncol(data)
  data = data[,c(1:7,(width-1):width,8:(width-2))]
  #SNPKModel = relmatLmer(PHE ~ PC1+PC2+PC3+PC4+PC5+(1|ID), phe_cov, relmat = list(ID = SNP.Kin),REML=F)
  
  iterations = ncol(data)-9
  cat(paste0(nrow(data)," samples remains and ",iterations," polymorphic genes detected.\n"))
  pb = txtProgressBar(max = iterations, style = 3)
  progress = function(n) setTxtProgressBar(pb, n)
  opts = list(progress = progress)
  
  genes_include = colnames(data)
  DT = foreach(j = 1:iterations, .combine = rbind, .options.snow = opts, .packages = c('lme4qtl','lme4')) %dopar% {
    col = j+9
    gene = genes_include[col]
    colnames(data)[col] = 'FH'
    
    data$FH = paste0("FH",data$FH)
    data$ID = as.factor(data$ID)
    data$FH = as.factor(data$FH)
    
    FHType = paste0(unique(as.character(data$FH)),collapse = "/")
    FHNumber = length(unique(as.character(data$FH)))
    FHLength = paste0(as.numeric(table(as.character(data$FH))),collapse = "/")
    
    if(FHNumber>=nrow(data)){
      colnames(data)[col] = gene
      c(pheName,gene,NA,FHType,FHNumber,FHLength,rep("HighVaried",8))
    } else {
      #FHKModel.Full = update(FHKModel,.~.+(1|FH))
      #SNPKModel = update(FHKModel,.~.+data$FH)
      model1 = lm(PHE ~ PC1+PC2+PC3+PC4+PC5,data)
      model2 = lmer(PHE ~ PC1+PC2+PC3+PC4+PC5+(1|FH),data,REML = F)
      #FHMpdel.Test = anova(FHKModel.Full,FHKModel)
      #SNPKModel.Test = anova(SNPKModel.Full,SNPKModel)
      FHMpdel.Test = anova(model2,model1)
      
      vf = as.data.frame(VarCorr(model2))[, c("grp", "vcov")]
      h2 = vf$vcov[which(vf$grp=='FH')]/sum(vf$vcov,na.rm = T)
      
      colnames(data)[col] = gene
      as.character(na.omit(c(pheName,gene,h2,FHType,FHNumber,FHLength,FHMpdel.Test$AIC,FHMpdel.Test$BIC,FHMpdel.Test$logLik,FHMpdel.Test$Chisq,FHMpdel.Test$`Pr(>Chisq)`)))
      #DT = rbind(DT,record)
    }
    
  }
  DT = as.data.frame(DT)
  colnames(DT) = c("Phe","Gene","Est_h2","FH_Type","FH_Number","FH_Length","AIC_Reduced","AIC_Full","BIC_Reduced","BIC_Full","LogLik_Reduced","LogLik_Full","Chisq","P")
  write.csv(DT,paste0("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/Seprate/Phe_",pheName,".GLMM.csv"),quote=F,row.names = F)
}

stopCluster(cl)


# pick up pleio genes
setwd("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/")

GLMs = list.files("./Seprate/","*csv",full.names = T)
phe = sapply(strsplit(GLMs,"/",fixed=T),'[[',3)
phe = sapply(strsplit(phe,".",fixed=T),'[[',1)
phe = sub("Phe_","",phe,fixed=T)
DT = data.frame()
meanChisq = c()
GCChisq = c()
for(i in 1:length(GLMs)){
  cat(paste0(i,"/",length(GLMs)," ...\n"))
  glm = GLMs[i]
  file = read.csv(glm,header = T)
  file$P = as.numeric(file$P)
  file$Chisq = as.numeric(file$Chisq)
  
  mean.chisq = mean(file$Chisq,na.rm=T)
  gc.chisq = median(file$Chisq,na.rm=T)/0.455
  meanChisq = c(meanChisq,mean.chisq)
  GCChisq = c(GCChisq,gc.chisq)
  # if(grepl("DXM",phe[i])){
  #   file$P = pchisq(file$Chisq/mean(file$Chisq,na.rm=T),1,0,lower.tail = F)
  # }
  
  # file2 = merge(file,gff,'Gene')
  # #file2 = file2[,c(1,16,17,15)]
  # file2 = file2[,c(1,16,17,14)]
  # colnames(file2) = c("SNP","Chromosome","Position",phe[i])
  # CMplot(file2,plot.type="q",box=FALSE,file="jpg",file.name="",dpi=300,
  #        conf.int=TRUE,conf.int.col=NULL,threshold.col="red",threshold.lty=2,
  #        file.output=TRUE,verbose=TRUE,width=5,height=5)
  
  #threshold = 0.05/length(which(!is.na(file$P)))
  threshold = 1e-5
  index = which(file$P<=threshold)
  target = file[index,]
  DT = rbind(DT,target)
}
statsDT = data.frame(Phe = phe,Mean = meanChisq, GC = GCChisq)
write.csv(statsDT,file = "GLMM.TestStatistic.Mean_LambdaGC.csv",quote=F,row.names = F)
write.csv(DT,file = "SignificantGene.20240909.csv",quote=F,row.names = F)

sumDT = as.data.frame(table(DT$Gene))
colnames(sumDT) = c("Gene","SigCount")
sumDT = sumDT[order(sumDT$SigCount,decreasing = T),]
write.csv(sumDT,file = "SignificantGene.summary.20240909.csv",quote=F,row.names = F)

# merged with annotation
sumDT = read.csv("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/SignificantGene.summary.20240909.csv",header=T)
anno = read.table("E:/QiGA/document/CottonLab/Cotton_GO/TM-1_v2.1_final.annotation.txt",header=T,sep = '\t',quote="")

dt = merge(sumDT,anno,by.x = 'Gene', by.y = 'ID',sort = F)
write.csv(dt,"E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/SignificantGene.summary.20240909.withAnnotation.csv",quote=F,row.names = F)

TF_home = c("AP2",'WRKY','ARF','ARR-B','B3','BBR-BPC','BES1','C2H2','C3H','CAMTA','CO-like','CPP',
            'DBB','Dof','E2F','EIL','ERF','FAR1','G2-like','GATA','GRAS','GRF','GeBP','HB-PHD','HB-other',
            'HD-ZIP','HRT','HSF','LBD','LFY','LSD','MADS','MYB','NAC','NF-X','NF-Y','NZZ','SPL','Nin-like',
            'RAV','S1Fa','SAP','SBP','SRS','STAT','TALE','TCP','Trihelix','VOZ','ZF-HD','WOX','Whirly',
            'YABBY','bHLH','bZIP')
DT = data.frame()
for(i in 1:length(TF_home)){
  TF = TF_home[i]
  target = grep(TF,dt$ArabDesc)
  targetDT = dt[target,]
  DT = rbind(DT,targetDT)
}
if(length(which(duplicated(DT$Gene)))>0){
  DT = DT[-which(duplicated(DT$Gene)),]
}
write.csv(DT,"E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GLMM/SignificantGene.summary.20240909.TF.csv",quote=F,row.names = F)
