library(data.table)
GT = fread("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/FHAP_whole.type.2.txt", sep = '\t', data.table = F, stringsAsFactors = F)
colnames(GT) = sub("0_","",colnames(GT))
GT1 = t(GT[,-1])
colnames(GT1) = GT$Gene
GT1 = data.frame(ID = rownames(GT1),GT1)
geneList = colnames(GT1)[-1]

FH.group = data.frame()
for(i in 1:6){
  temp = read.table(paste0("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/Phy/Group/G", i, ".FH.sample.list.txt"),header=F)
  temp$FHGroup = paste0("G",i)
  
  FH.group = rbind(FH.group,temp)
}
colnames(FH.group)[1] = 'ID'

Group1 = fread("E:/QiGA/document/CottonLab/XLD/population_info/Detailed_Information.txt",sep = '\t',quote="",header=T)
Group1 = Group1[,c(1,19)]
colnames(Group1)[1] = 'ID'
Group1 = merge(Group1,FH.group,'ID',all.y = T)
colnames(Group1) = c('ID','Source Region','Phylogenetic Group')

GeneFlow = data.frame()
for(i in 1:length(geneList)){
  gene = geneList[i]
  dt = GT1[,c("ID",gene)]
  dt = merge(dt,Group1,'ID',all.x=T)
  
  stats1 = as.matrix(table(dt[,gene],dt$`Source Region`))
  stats1 = stats1>1
  countF1 = apply(stats1,1,function(x){
    length(which(!x))
  })
  stats1 = stats1[which(countF1!=ncol(stats1)),,drop=F]
  if(nrow(stats1)==0){
    next
  }
  hapList1 = paste0(gene,"_",sprintf("HAP%04d", as.numeric(rownames(stats1))))
  flow1 = apply(stats1,1,function(x){
    c(length(which(x)),paste0(colnames(stats1)[which(x)],collapse = ','))
  })
  
  stats2 = as.matrix(table(dt[,gene],dt$`Phylogenetic Group`))
  stats2 = stats2>=10
  countF2 = apply(stats2,1,function(x){
    length(which(!x))
  })
  stats2 = stats2[which(countF2!=ncol(stats2)),,drop=F]
  if(nrow(stats2)==0){
    next
  }
  hapList2 = paste0(gene,"_",sprintf("HAP%04d", as.numeric(rownames(stats2))))
  flow2 = apply(stats2,1,function(x){
    c(length(which(x)),paste0(colnames(stats2)[which(x)],collapse = ','))
  })
  
  flow1 = cbind(hapList1, t(flow1))
  flow2 = cbind(hapList2, t(flow2))
  colnames(flow1) = c("HapName","NumbersAppearRegion","SourceRegion")
  colnames(flow2) = c("HapName","NumbersAppearPhylogenetic","PhylogeneticGroup")
  flow = merge(flow1,flow2,'HapName',all=T)
  GeneFlow = rbind(GeneFlow,flow)
}
openxlsx::write.xlsx(GeneFlow,'E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GeneFlow/GeneFlow.Summary.2support.xlsx')

DT = openxlsx::read.xlsx("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GeneFlow/GeneFlow.Summary.2support.xlsx")
DT$NumbersAppearRegion = as.numeric(DT$NumbersAppearRegion)
DT1 = DT[which(DT$NumbersAppearRegion>1),]

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
weight = function(x){
  exp(-1*((x-2)^2)/10)
}
summaryDT = data.frame()
regions = c("GB","Central America","North America","South America","SC (CHN)",
         "Europe","Asia","Africa",
         "NC (CHN)","FSU","SWC (CHN)","NWC (CHN)", "USA","YZR (CHN)","YER (CHN)")
bootstrap = 10
for(i in 1:length(regions)){
  region1 = regions[i]
  re1 = region1
  if(re1 == "CHN"){
    re1 = ",CHN,"
  }
  hits1 = grep(re1,DT1$SourceRegion,fixed=T)
  for(j in i:length(regions)){
    region2 = regions[j]
    re2 = region2
    if(re2 == 'CHN'){
      re2 = ",CHN,"
    }
    hits2 = grep(re2,DT1$SourceRegion,fixed=T)
    hits = intersect(hits1,hits2)
    
    dt = DT1[hits,]
    countStat = as.data.frame(table(dt$NumbersAppearRegion))
    countStat$Var1 = as.numeric(as.character(countStat$Var1))
    
    #score = sum(weight(countStat$Var1)*countStat$Freq)/(sum(countStat$Freq))
    score = sum(weight(countStat$Var1)*countStat$Freq)/(nrow(DT1))
    print(sum(countStat$Freq))
    
    bootValue = c()
    for(l in 1:bootstrap){
      temp = dt[sample(1:nrow(dt),floor(0.1*nrow(dt)),replace = T),]
      temp.countStat = as.data.frame(table(temp$NumbersAppearRegion))
      temp.countStat$Var1 = as.numeric(as.character(temp.countStat$Var1))
      
      temp.score = sum(weight(temp.countStat$Var1)*temp.countStat$Freq)/(sum(temp.countStat$Freq)*2)
      bootValue = c(bootValue,temp.score)
    }
    stats = t.test(bootValue,alternative = 'greater',mu = 0)
    meanV = as.numeric(stats$estimate)
    meanP = stats$p.value
    record = c(region1,region2,score,meanV,meanP)
    summaryDT = rbind(summaryDT,record)
  }
}
colnames(summaryDT) = c("Region1","Region2","Score","BootScore","BootP")
summaryDT$Score = as.numeric(summaryDT$Score)
summaryDT$BootP = as.numeric(summaryDT$BootP)
write.table(summaryDT,"E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GeneFlow/FlowIndex.fullTable.txt",quote=F,col.names = T,row.names = F,sep = '\t')

summaryMat = reshape2::acast(summaryDT,Region1~Region2,value.var = 'Score')
summaryMat = summaryMat[regions,regions]
summaryMat = apply(summaryMat,1,function(x){
  x/max(x,na.rm=T)
})
summaryMat = fill_missing_symmetric(summaryMat)

pMat = summaryMat
pMat = pMat>=0.15
pMat[which(pMat)] = 1e-10
pMat[which(!pMat)] = 1
summaryMat1 = summaryMat
summaryMat1[which(summaryMat1<0.15)] = 0
corrplot(summaryMat1,p.mat = pMat,tl.srt = 45,tl.cex = 1,tl.col = 'black',number.cex = 0.8,addCoef.col = 'black',insig = 'blank',is.corr = F,method = 'square',type = 'lower',col = COL1('Oranges', 100),diag = T)
pdf("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GeneFlow/FlowIndex.pdf",width = 10,height = 10)
corrplot(summaryMat1,p.mat = pMat,tl.srt = 45,tl.cex = 1,tl.col = 'black',number.cex = 0.8,addCoef.col = 'black',insig = 'blank',is.corr = F,method = 'square',type = 'lower',col = COL1('Oranges', 100),diag = T)
dev.off()
write.table(summaryMat,"E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GeneFlow/FlowIndex.txt",quote=F,col.names = T,row.names = F,sep = '\t')

# plot map 
library(ggplot2)
library(ggthemes)
library(ggspatial)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

world <- ne_countries(scale = "medium", returnclass = "sf") %>% st_make_valid()

degree = 150
target_crs <- st_crs(paste0("+proj=eqc +x_0=0 +y_0=0 +lat_0=0 +lon_0=",degree))

# define a long & slim polygon that overlaps the meridian line & set its CRS to match
# that of world
# Centered in lon 133

offset <- 180 - degree

polygon <- st_polygon(x = list(rbind(
  c(-0.0001 - offset, 90),
  c(0 - offset, 90),
  c(0 - offset, -90),
  c(-0.0001 - offset, -90),
  c(-0.0001 - offset, 90)
))) %>%
  st_sfc() %>%
  st_set_crs(4326)

world2 <- world %>% st_difference(polygon)
# Transform
world3 <- world2 %>% st_transform(crs = target_crs)

sites <- data.frame(longitude = c(0, -102,-100,-70,-68,110,108,95,70,25), 
                    latitude = c(45, 22,40,20,5,20,33,25,47,15),
                    Labels = c("EUR","NAL","USA","CAL","SAL","SCL","CHN","ASL","FSU","AFR")) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs = target_crs)

time1 = 1.05
time2 = 3
arrows1 <- data.frame(XSTART = c(-70,-70,-102,-68,-102,110,-100,70,110,-70,0,95,-100)*time1, 
                     YSTART = c(20,20,22,5,22,20,40,47,20,20,45,25,40)*time1,
                     Weight = c(0.7695,0.3313,0.2928,0.2928,0.4305,0.1633,0.5119,0.0992,0.1,0.1,0.1,0.1,0.1)*time2,
                     Colors = c("#3C5488FF","#3C5488FF","#3C5488FF","#3C5488FF","#3C5488FF","#E64B35FF","#E64B35FF","#E64B35FF","#3C5488FF","#3C5488FF","#E64B35FF","#3C5488FF","#E64B35FF") ) %>%
  st_as_sf(coords = c("XSTART", "YSTART"), crs = 4326) %>%
  st_transform(crs = target_crs)
arrows2 <- data.frame(XEND = c(-102,-68,-68,-102,110,108,108,108,95,0,25,25,70)*time1,
                      YEND = c(22,5,5,22,20,33,33,33,25,45,15,15,47)*time1) %>%
  st_as_sf(coords = c("XEND", "YEND"), crs = 4326) %>%
  st_transform(crs = target_crs)

arrows1$XSTART = st_coordinates(arrows1$geometry)[,1]
arrows1$YSTART = st_coordinates(arrows1$geometry)[,2]
arrows2$XEND = st_coordinates(arrows2$geometry)[,1]
arrows2$YEND = st_coordinates(arrows2$geometry)[,2]
arrows = cbind(arrows1[,c(1,2,4,5)],arrows2[,c(2,3)])
#st_bbox(world3)
p1 = 
  ggplot(data = world3) + 
  geom_sf(fill= 'antiquewhite') +
  #geom_sf(data = sites, aes(), size = 4, shape = 23, fill = "darkred") +
  geom_text(data = sites, aes(x = st_coordinates(geometry)[,1], 
                              y = st_coordinates(geometry)[,2], 
                              label = Labels), size = 4, nudge_y = 5,color='black',fontface='bold') +
  coord_sf(xlim = c(st_bbox(world3)[1]*0.98, st_bbox(world3)[3])*0.98, 
           ylim = c(st_bbox(world3)[2]*0.65, st_bbox(world3)[4])*0.98, expand = FALSE)+
  annotate(geom = 'curve',
           x = arrows$XSTART,
           xend = arrows$XEND,
           y = arrows$YSTART,
           yend = arrows$YEND,
           curvature = .3,
           color = arrows$Colors,
           size = 1
           #arrow = arrow(length = unit(3, 'mm'))
           )+
  xlab('Longitude') + 
  ylab('Latitude') + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = 'dashed', size = 0.25), 
        panel.background = element_rect(fill = 'aliceblue'))
ggsave("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/GeneFlow/WorldMap2.pdf",width = 5,height = 4)
