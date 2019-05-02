library(ggplot2)
library(gplots)
library("pheatmap")
library("RColorBrewer")
library("viridis")

################
# Composite PC #
################

PC2 = read.table('.../Tyne_Litc_PC2_loadings.txt')
PC5 = read.table('.../Tyne_Litc_PC5_loadings.txt')

PC = data.frame(x=0.145*PC2$x+0.063*-1*PC5$x)
rownames(PC) = rownames(PC2)
extremePCglobal = rownames(PC)[which(PC$x<quantile(PC$x,probs=c(0.01)) | PC$x>quantile(PC$x,probs=c(0.99)))]

################################
# compare to salinity DE genes #
################################
library(DESeq2)

counts = read.table('.../expression/cuffnorm/STAR/UCSCgtfSalinity/genes.count_table',row.names=1,header=T,sep='\t')

# take average of parents replicated on two lanes
counts$c172_P_532_F = apply(counts[,c('c172_P_532_F_lane3_0','c172_P_532_F_lane8_0')],1,mean)
counts$c172_P_533_M = apply(counts[,c('c172_P_533_M_lane3_0','c172_P_533_M_lane8_0')],1,mean)

sel = c('c172_P_532_F',
'c172_P_533_M',
'c172_F1_20_F_0',
'c172_F1_20_M_0',
'c172_F1_04_F_0',
'c172_F1_04_M_0',
'c172_F1_01_F_0',
'c172_F1_01_M_0',
'c172_F1_13_F_0',
'c172_F1_13_M_0',
'c172_F1_05_F_0',
'c172_F1_05_M_0',
'c172_F1_10_F_0',
'c172_F1_10_M_0')

sex = c('F',
'M',
'F',
'M','F',
'M','F',
'M','F',
'M','F',
'M','F',
'M')

sal = c('P_mar','P_fw',rep('3.5',4),rep('0.2',4),rep('35',4))

subCounts = counts[,sel]
design = data.frame(row.names=sel,sal=sal, sex=sex)

iCounts = round(subCounts)
iCounts = DESeqDataSetFromMatrix(countData=iCounts, colData=design, design=~sal+sex)
iCounts = DESeq(iCounts)
normCounts = assay(iCounts,normalized=T)

# overall
res = results(iCounts)

# contrasts
res3vs02 = results(iCounts, contrast=c("sal","3.5","0.2"))
summary(res3vs02,alpha=0.01)
sig3vs02 = subset(res3vs02,padj<0.01)

res3vs35 = results(iCounts, contrast=c("sal","3.5","35"))
summary(res3vs35,alpha=0.01)
sig3vs35 = subset(res3vs35,padj<0.01)

res35vs02 = results(iCounts, contrast=c("sal","35","0.2"))
summary(res35vs02,alpha=0.01)
sig35vs02 = subset(res35vs02,padj<0.01)

#postscript('/fml/chones/projects/PJ035_JPV_CisTrans/ASE/salinityVennFDR1.ps')
venn(list(rownames(sig3vs02),rownames(sig3vs35),rownames(sig35vs02)))
salinityDE = unique(c(rownames(sig3vs02),rownames(sig3vs35),rownames(sig35vs02)))
#write.table(c(rownames(sig3vs02),rownames(sig3vs35),rownames(sig35vs02)),'salinityDE.txt',row.names=F,col.names=F,quote=F)
#dev.off()

parallelSalinityDE = attr(venn(list(extremePCglobal,unique(c(rownames(sig3vs02),rownames(sig3vs35),rownames(sig35vs02))))),"intersect")$'A:B'
#write.table(parallelSalinityDE,'parallelSalinityDE.txt',row.names=F,col.names=F,quote=F)

# overlap probability between PC2 outliers and salinity DE genes
total = length(na.omit(res$padj)) # --> tested by DESeq2 and PCA
salN = length(unique(c(rownames(sig3vs02),rownames(sig3vs35),rownames(sig35vs02))))
PC2N = length(extremePCglobal)
overlap = length(parallelSalinityDE)
phyper(overlap - 1, salN, total - salN, PC2N, lower.tail=F)
# sal FDR 1% 3.658975e-14
# sal FDR 5% 5.602629e-13
# sal FDR 10% 1.04518e-13


#########################################
# expression level relative to rowmeans #
#########################################
#---------------------------
# variance stabilized counts
vCounts = varianceStabilizingTransformation(iCounts,blind=FALSE)
normCounts = assay(vCounts,normalized=T)

# exclude parents
sel = c(
'c172_P_532_F',
'c172_P_533_M',
'c172_F1_20_F_0',
'c172_F1_20_M_0',
'c172_F1_04_F_0',
'c172_F1_04_M_0',
'c172_F1_01_F_0',
'c172_F1_01_M_0',
'c172_F1_13_F_0',
'c172_F1_13_M_0',
'c172_F1_05_F_0',
'c172_F1_05_M_0',
#'c172_F1_10_F_0',
'c172_F1_10_M_0')

countMat = data.frame(normCounts[parallelSalinityDE,sel]-rowMeans(normCounts[parallelSalinityDE,sel]))
colnames(countMat) <- design[sel,'sal']

###################
#------------------
# Boxplot
library(reshape2)
countMatPlot = melt(data.frame(countMat,gene=rownames(countMat)))

ggplot(countMatPlot[which(countMatPlot$variable %in% c('P_mar','P_fw') == F),])+
geom_boxplot(aes(x=gene,y=value))+
facet_wrap(~variable)


###################
#------------------
# Heatmap
paletteLength <- 50
myColor <- colorRampPalette(c( "turquoise", "black", "yellow"))(paletteLength)
myBreaks <- c(seq(min(countMat), 0, length.out=ceiling(paletteLength/2) + 1),seq(max(countMat)/paletteLength, max(countMat), length.out=floor(paletteLength/2)))

#pdf('.../salinityHeatmapRelativeNormCountsParallelSalinityDE1FDR.pdf',height=5,width=7)
pheatmap(countMat,clustering_distance_cols='correlation',clustering_distance_rows='manhattan',col=myColor,border_color=NA)
#dev.off()

# candidate heatmap
aqp3a = 'XLOC_020754'
slc12a2 = 'XLOC_022471'
kcnj1a = 'XLOC_000145'
trpv6='XLOC_028455'

candMat = countMat[c(aqp3a,kcnj1a),]
paletteLength <- 50
myColor <- colorRampPalette(c( "turquoise", "black", "yellow"))(paletteLength)
myBreaks <- c(seq(min(candMat), 0, length.out=ceiling(paletteLength/2) + 1),seq(max(candMat)/paletteLength, max(candMat), length.out=floor(paletteLength/2)))
x11(h=2,w=5)
pheatmap(candMat,clustering_distance_cols='correlation',clustering_distance_rows='manhattan',col=myColor,border_color=NA,breaks=myBreaks)+geom_hline(yintercept=0,colour='grey78')
#dev.off()


###############
#--------------
# Correlation

# DESeq2 normalized counts relative to rowmeans - expression profiles
countDist = cor(countMat[parallelSalinityDE,],method='spearman')
countDistMatrix <- as.matrix(countDist)
rownames(countDistMatrix) <- design[sel,'sal']
colnames(countDistMatrix) <- design[sel,'sal']
myBreaks <- c(seq(min(countDistMatrix), 0, length.out=ceiling(paletteLength/2) + 1),seq(max(countDistMatrix)/paletteLength, max(countDistMatrix), length.out=floor(paletteLength/2)))
paletteLength <- 50
colors <- colorRampPalette(c("turquoise", "black", "yellow"))(paletteLength)
#pdf('.../salinityCorrelationNormCountsParallelSalinityDE10FDR.pdf')
pheatmap(countDistMatrix,
clustering_distance_rows='correlation',
clustering_distance_cols='correlation',
col=colors,
breaks=myBreaks)
#dev.off()

# only transcipts where ASE ratio was tested also
xlocSel = read.table('candidateSalinityASExloc.txt',stringsAsFactors=F)
xlocSel = xlocSel[,1]


#---------------------------------------------
# FPKM counts --> can be compared across genes
fpkmCounts = read.table('.../expression/cuffnorm/STAR/UCSCgtfSalinity/genes.fpkm_table',row.names=1,header=T,sep='\t')
# take average of parents replicated on two lanes
fpkmCounts$c172_P_532_F = apply(fpkmCounts[,c('c172_P_532_F_lane3_0','c172_P_532_F_lane8_0')],1,mean)
fpkmCounts$c172_P_533_M = apply(fpkmCounts[,c('c172_P_533_M_lane3_0','c172_P_533_M_lane8_0')],1,mean)

# transform counts on log-scale to reduce the effect of outlier genes with large variation
# substract the log(FPKM) level of each sample from the row mean --> average across all samples
# this gives a "profile" --> whether sample had lower of higher expression relative to others
sel = c(
'c172_P_532_F',
'c172_P_533_M',
'c172_F1_20_F_0',
'c172_F1_20_M_0',
'c172_F1_04_F_0',
'c172_F1_04_M_0',
'c172_F1_01_F_0',
'c172_F1_01_M_0',
'c172_F1_13_F_0',
'c172_F1_13_M_0',
'c172_F1_05_F_0',
'c172_F1_05_M_0',
#'c172_F1_10_F_0',
'c172_F1_10_M_0')

sal = c('P_mar','P_fw',rep('3.5',4),rep('0.2',4),rep('35',3))

countMat = data.frame(log(fpkmCounts[parallelSalinityDE,sel]+0.01)-rowMeans(log(fpkmCounts[parallelSalinityDE,sel]+0.01)))
colnames(countMat) <- design[sel,'sal']

# exclude extreme values for better plotting
#select = function(x){all(x>-4)} # excludes 8 rows
#countMat = countMat[apply(countMat,1,select),]

# plot
countDist = cor(countMat,method='spearman')
countDistMatrix <- as.matrix(countDist)
rownames(countDistMatrix) <- paste(design[sel,'sal'],1:13,sep='_')
colnames(countDistMatrix) <- paste(design[sel,'sal'],1:13,sep='_')
paletteLength <- 50
myBreaks <- c(seq(-1, 0, length.out=ceiling(paletteLength/2) + 1),seq(max(countDistMatrix)/paletteLength, max(countDistMatrix), length.out=floor(paletteLength/2)))
colors <- colorRampPalette(c("turquoise", "black", "yellow"))(paletteLength)
# order rows
countDistMatrix = countDistMatrix[c(8,10,7,9,2,4,5,6,3,13,1,11,12),c(8,10,7,9,2,4,5,6,3,13,1,11,12)]
pdf('.../salinityCorrelationNormCountsParallelSalinityDE1FDR.pdf')
pheatmap(countDistMatrix,
#clustering_distance_rows='euclidean',
#clustering_distance_cols='euclidean',
col=colors,
breaks=myBreaks,
cluster_rows=F, cluster_cols=F
)
dev.off()


#######
# PCA #
#######
library(ggbiplot)

subDesign = design[which(rownames(design) %in% c("c172_P_532_F","c172_P_533_M")==F),]

# PCA based on F1's and PC outlier genes overlapping salinity DE genes
pca = prcomp(t(countMat[salinityDE,which(colnames(countMat) %in% c("P_mar","P_fw")==F)]),scale=F)

subCountsFS = normCounts[salinityDE,which(colnames(normCounts) %in% c("P_mar","P_fw"))]
FSproject = data.frame(scale(t(subCountsFS), pca$center, pca$scale) %*% pca$rotation[,c(1,2)],
eco = c('M','F'))
percentVar <- round(100 * pca$sdev^2/sum(pca$sdev^2))

g = ggbiplot(pca, var.axes=F, choices=c(1,2), obs.scale=1, var.scale=1, groups=as.factor(subDesign$sal), ellipse=FALSE, ellipse.prob=0.95, alpha=0)
g = g + scale_color_manual(values=c('darkblue','darkred','grey','black','black'))
g = g + theme(axis.text=element_text(hjust=1,size=24,colour='black'),panel.background=element_rect(fill='white',colour='black'),panel.grid.major=element_line(colour='grey78'),text=element_text(size=22),legend.position="right")
g = g + geom_point(data=as.data.frame(pca$x), aes(PC1, PC2, color=subDesign$sal),size=8)
g = g + geom_point(data=FSproject,aes(PC1, PC2,shape=eco),size=8)
g
ggsave('.../salinityPCANormCountsSalinityDE10FDR.ps')


# PCA based on relative expression levels
pca = prcomp(t(normCounts[salinityDE,which(colnames(normCounts) %in% c("c172_P_532_F","c172_P_533_M")==F)]),scale=F)

subCountsFS = normCounts[salinityDE,which(colnames(normCounts) %in% c("c172_P_532_F","c172_P_533_M"))]
FSproject = data.frame(scale(t(subCountsFS), pca$center, pca$scale) %*% pca$rotation[,c(1,2)],
eco = c('M','F'))
percentVar <- round(100 * pca$sdev^2/sum(pca$sdev^2))

g = ggbiplot(pca, var.axes=F, choices=c(1,2), obs.scale=1, var.scale=1, groups=as.factor(subDesign$sal), ellipse=FALSE, ellipse.prob=0.95, alpha=0)
g = g + scale_color_manual(values=c('darkblue','darkred','grey','black','black'))
g = g + theme(axis.text=element_text(hjust=1,size=24,colour='black'),panel.background=element_rect(fill='white',colour='black'),panel.grid.major=element_line(colour='grey78'),text=element_text(size=22),legend.position="right")
g = g + geom_point(data=as.data.frame(pca$x), aes(PC1, PC2, color=subDesign$sal),size=8)
g = g + geom_point(data=FSproject,aes(PC1, PC2,shape=eco),size=8)
g
ggsave('.../salinityPCANormCountsSalinityDE10FDR.ps')

#################
# example genes #
#################


# slc12a1/slc12a2 on chrXIV (best example)
# ion secretion -> upregulation in Mar
xloc='XLOC_022471'
plotData = data.frame(value=c(unlist(countMat[xloc,7:10]),
unlist(countMat[xloc,3:6]),
unlist(countMat[xloc,11:13])),
sal=c(rep('0.2',4),rep('3.5',4),rep('35',3)))
plotData$sal = factor(plotData$sal,levels=c("0.2","3.5","35"))
x11(h=2,w=5)
ggplot(plotData,aes(sal,value))+geom_boxplot()+theme(panel.background=element_rect(fill='white',colour='grey78'),legend.position="right",text=element_text(size=16,colour='black'))+labs(title="slc12a1/slc12a2")
ggsave('.../slc12a1salinityBoxplot.ps')

# aqp3a
xloc='XLOC_020754'
plotData = data.frame(value=c(unlist(countMat[xloc,7:10]),
unlist(countMat[xloc,3:6]),
unlist(countMat[xloc,11:13])),
sal=c(rep('0.2',4),rep('3.5',4),rep('35',3)))
plotData$sal = factor(plotData$sal,levels=c("0.2","3.5","35"))
x11(h=2,w=5)
ggplot(plotData,aes(sal,value))+geom_boxplot()+theme(panel.background=element_rect(fill='white',colour='grey78'),legend.position="right",text=element_text(size=16,colour='black'))+labs(title="aqp3a")
ggsave('.../aqp3asalinityBoxplotFDR1.ps')

# trpv6
xloc='XLOC_028455'
plotData = data.frame(value=c(unlist(countMat[xloc,7:10]),
unlist(countMat[xloc,3:6]),
unlist(countMat[xloc,11:13])),
sal=c(rep('0.2',4),rep('3.5',4),rep('35',3)))
plotData$sal = factor(plotData$sal,levels=c("0.2","3.5","35"))
x11(h=2,w=5)
ggplot(plotData,aes(sal,value))+geom_boxplot()+theme(panel.background=element_rect(fill='white',colour='grey78'),legend.position="right",text=element_text(size=16,colour='black'))+labs(title="trpv6")
ggsave('.../trpv6salinityBoxplotFDR1.ps')


# kcnj1a.6
xloc = 'XLOC_000145'
plotData = data.frame(value=c(unlist(countMat[xloc,7:10]),
unlist(countMat[xloc,3:6]),
unlist(countMat[xloc,11:13])),
sal=c(rep('0.2',4),rep('3.5',4),rep('35',3)))
plotData$sal = factor(plotData$sal,levels=c("0.2","3.5","35"))
x11(h=2,w=5)
ggplot(plotData,aes(sal,value))+geom_boxplot()+theme(panel.background=element_rect(fill='white',colour='grey78'),legend.position="right",text=element_text(size=16,colour='black'))+labs(title="kcnj1a.3")
ggsave('.../kcnj1salinityBoxplotFDR1.ps')


# combined
aqp3a = 'XLOC_020754'
slc12a2 = 'XLOC_022471'
kcnj1a = 'XLOC_000145'
trpv6='XLOC_028455'
unknown='XLOC_001256'
NADH5C='XLOC_011835'
mucin='XLOC_010054'
cands = c(unknown,NADH5C,slc12a2,aqp3a)

plotData = melt(data.frame(t(countMat[cands,3:13]),sal=c(rep('3.5',4),rep('0.2',4),rep('35',3))))
plotData$sal = factor(plotData$sal,levels=c("0.2","3.5","35"))
plotData$variable = factor(plotData$variable,levels=c("XLOC_020754","XLOC_001256","XLOC_022471","XLOC_011835"))

x11(h=5,w=2)
ggplot(plotData,aes(sal,value))+geom_boxplot()+theme(panel.background=element_rect(fill='white',colour='grey78'),legend.position="right",text=element_text(size=16,colour='black'))+facet_wrap(~variable,ncol=1,scales='free')+geom_hline(yintercept=0,colour='black',linetype='dotted')
ggsave('.../candidateSalinityBoxplotFDR1.ps')





