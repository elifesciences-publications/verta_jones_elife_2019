library(gtools)
library(gplots)

#--------------------------------
qThreshold = 0.2

#--------------------------------
# DE test results - TYNE
tyneDEtest = read.table('.../expression/cuffdiff/STAR/fw_vs_mar_Tyne_UCSCgtf/gene_exp.diff',header=T,stringsAsFactors=F,sep='\t')
tyneDEsig = tyneDEtest[which(tyneDEtest$q_value<qThreshold),]
tyneDE = tyneDEsig$gene_id
tyneFC = tyneDEtest$log2.fold_change.
names(tyneFC) = tyneDEtest$gene_id
tyneLR = foldchange2logratio(tyneFC)

#--------------------------------
# DE test results - LITC
litcDEtest = read.table('.../expression/cuffdiff/STAR/fw_vs_mar_Litc_UCSCgtf/gene_exp.diff',header=T,stringsAsFactors=F,sep='\t')
litcDEsig = litcDEtest[which(litcDEtest$q_value<qThreshold),]
litcDE = litcDEsig$gene_id
litcFC = litcDEtest$log2.fold_change.
names(litcFC) = litcDEtest$gene_id
litcLR = foldchange2logratio(litcFC)

#--------------------------------
# compare Tyne and Litc
v = venn(list(tyneDE,litcDE))
sharedSig = attr(v,"intersections")$"A:B"
litcSpecific = attr(v,"intersections")$"B"
tyneSpecific = attr(v,"intersections")$"A"

# plot
plot(litcFC,tyneFC,pch=20,cex=0.2)
points(litcFC[litcSpecific],tyneFC[litcSpecific],col='blue')
points(litcFC[tyneSpecific],tyneFC[tyneSpecific],col='green')
points(litcFC[sharedSig],tyneFC[sharedSig],col='red')

# parallel
parallel = sharedSig[which(sign(litcFC[sharedSig]) == sign(tyneFC[sharedSig]))]

# anti-parallel
antiParallel = sharedSig[which(sign(litcFC[sharedSig]) != sign(tyneFC[sharedSig]))]

#--------------------------------
# DE test results - both together
parallelDEtest = read.table('.../expression/cuffdiff/STAR/fw_vs_mar_Litc_Tyne_UCSCgtf/gene_exp.diff',header=T,stringsAsFactors=F,sep='\t')
rownames(parallelDEtest) = parallelDEtest$gene_id
parallelDEsig = parallelDEtest[which(parallelDEtest$q_value<qThreshold),]
parallelDE = parallelDEsig$gene_id

# parallel
parallel = parallelDE[which(sign(litcFC[parallelDE]) == sign(tyneFC[parallelDE]))]

# plot
plot(litcFC,tyneFC,pch=20,cex=0.2)
points(litcFC[litcDE],tyneFC[litcDE],col='blue')
points(litcFC[tyneDE],tyneFC[tyneDE],col='green')
points(litcFC[parallel],tyneFC[parallel],col='red')

# compPC
sameSign = names(litcLR)[which(sign(litcLR) == sign(tyneLR))]
PC2 = read.table('.../ASE/UCSCgtf/Tyne_Litc_PC2_loadings.txt')
PC5 = read.table('.../ASE/UCSCgtf/Tyne_Litc_PC5_loadings.txt')
PC = data.frame(x=0.145*PC2$x+0.063*-1*PC5$x)
rownames(PC) = rownames(PC2)
extremePCglobal = rownames(PC)[which(PC$x<quantile(PC$x,probs=c(0.01)) | PC$x>quantile(PC$x,probs=c(0.99)))]
PC$zPC = scale(PC$x)
extremePCglobalSameSign = extremePCglobal[which(extremePCglobal %in% sameSign)]

# plot
plotData = data.frame(litc=litcFC,tyne=tyneFC,parallel=names(litcLR)%in%parallel,PC=names(litcFC)%in%extremePCglobal)
g = ggplot(plotData[which(plotData$parallel == F),],aes(litc,tyne))
g = g + geom_point(color='grey')
g = g + geom_point(data=plotData[which(plotData$parallel == T),],aes(litc,tyne),color='black')
g = g + theme(axis.text=element_text(hjust=1,size=24,colour='black'),panel.background=element_rect(fill='white',colour='black'),text=element_text(size=22),legend.position="right")
#g = g + geom_point(data=plotData[which(plotData$PC == T),],aes(litc,tyne),color='blue')
g
ggsave(".../ASE/UCSCgtf/LitcTyneParallelDE.ps", device=cairo_ps)





