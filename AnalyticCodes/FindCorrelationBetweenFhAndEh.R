# fit Eh and Fh function
dt = read.table("E:/QiGA/document/CottonLab/FH/ResultsAllVariants/FHAP_whole.GeneSTAT.txt",header=T)

model1 = lm(HapN~EH,data = dt)
summary(model1)

model2 = lm(HapN~poly(EH,2),data = dt)
summary(model2)

model3 = lm(HapN~poly(EH,3,raw = T),data = dt)
summary(model3)

pdf("E:/QiGA/document/CottonLab/FH/EH_vs_FH.pdf",height = 6,width = 6)
plot(dt$EH,dt$HapN,frame.plot = F,ylim = c(0,4000),
     ylab = "Number of functional haplotype (FH)",
     xlab = expression(paste("Shannon's equitability (",italic(E[h]),")")),
     pch = 16)
seq_test = runif(10000)
seq_predict = as.numeric(predict(model3,newdata = data.frame(EH=seq_test)))
points(seq_test,seq_predict,col='red',pch = 16)
dev.off()