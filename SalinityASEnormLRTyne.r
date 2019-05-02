library(DESeq2)
library(GenomicRanges)
library(ggplot2)
library(reshape)

##################################
# allele-specific counts over SNPs
cov = read.table('.../aseReadCounts_c172_allF1_infSites_alleleCov_DESeqNormalized_STAR_duprem.txt',header=T)
ref = cov[,which(grepl('REF',colnames(cov)))]
alt = cov[,which(grepl('ALT',colnames(cov)))]

# polarise counts
genoInf = read.table('.../c172_P_genomic_recalibrated_filtered_SNPs_informative_minCov10_genotypes.txt',header=T,stringsAsFactors=F)
rownames(genoInf) = paste(genoInf[,1],genoInf[,2],sep=':')
fmREF = rownames(genoInf)[which(genoInf$female_P == '0/0')]
mREF = rownames(genoInf)[which(genoInf$male_P == '0/0')]
fmALT = rownames(genoInf)[which(genoInf$female_P == '1/1')]
mALT = rownames(genoInf)[which(genoInf$male_P == '1/1')]

# create new data frames that are polarized allele counts for female and male
fm = data.frame()
fm = ref[which(rownames(ref) %in% fmREF),]
fm[nrow(fm)+1:length(which(rownames(alt) %in% fmALT)),] = alt[which(rownames(alt) %in% fmALT),]
colnames(fm) = gsub('REFCOV','FM',colnames(fm))

m = data.frame()
m = ref[which(rownames(ref) %in% mREF),]
m[nrow(m)+1:length(which(rownames(alt) %in% mALT)),] = alt[which(rownames(alt) %in% mALT),]
colnames(m) = gsub('REFCOV','M',colnames(m))

normCounts = cbind(fm,m)
design = data.frame(row.names=colnames(cbind(fm,m)), salinity=c(rep('3.5',times=4),rep('35',times=2),rep('0.2',times=4),rep('35',times=2)), allele=c(rep('FM',times=12),rep('M',times=12)))

#write.table(normCounts,'.../c172_polarised_allele_counts_allF1_DESeqNormalized_STAR_duprem_allSNP.txt')

library(gtools)
normFC = data.frame()
normLR = data.frame()
for (i in rownames(normCounts)){
    for (y in colnames(fm)){
        normFC[i,y] = foldchange(fm[i,y],m[i,sub('_FM','_M',y)])
        normLR[i,y] = foldchange2logratio(foldchange(fm[i,y],m[i,sub('_FM','_M',y)]))
    }
}
colnames(normFC) = unlist(lapply(colnames(fm),function(x){sub('_FM','',x)}))
colnames(normLR) = unlist(lapply(colnames(fm),function(x){sub('_FM','',x)}))

#write.table(normFC,'.../c172_polarised_allele_foldchange_allF1_DESeqNormalized_STAR_duprem_allSNP.txt')


###########################################################
#----------------------------------------------------------
# ASE in parallel diverged genes / salinity treatment genes

library(gtools)
normFC = read.table('.../c172_polarised_allele_foldchange_allF1_DESeqNormalized_STAR_duprem_allSNP.txt')
# exclude 10 Female sample --> outlier
normFC = normFC[,colnames(normFC) %in% 'c172_F1_10_F' == F]
normFC[mapply(is.infinite,normFC)] = NA
normLR = mapply(foldchange2logratio,normFC)
rownames(normLR) = rownames(normFC)

# pick SNPs that were tested for effects (one SNP per xloc per F1)
effects172 = read.table('.../AseReadCountsEffects172_bestQvalue_cisParentTrans_STAR_duprem_nosexFDR10.txt',header=T,stringsAsFactors=F)
normLR = normLR[which(rownames(normLR) %in% paste(effects172$seqnames,effects172$start,sep=':')),]

# turn into genomicRange
getChr = function(x){unlist(unlist(strsplit(x,':'))[1])}
getPos = function(x){unlist(as.numeric(unlist(strsplit(x,':'))[2]))}
normLRrange = makeGRangesFromDataFrame(data.frame(chrom=unlist(lapply(rownames(normLR),getChr)),start=unlist(lapply(rownames(normLR),getPos)),end=unlist(lapply(rownames(normLR),getPos))+1,normLR[,1:ncol(normLR)]),keep.extra.columns=T)


# xloc annotation
geneFeatures = read.table('.../assembly/merged_assembly_STAR_UCSCgtf/mergedParsed.txt',stringsAsFactors=F,sep='\t',skip=1,fill=NA)
colnames(geneFeatures) = c('chrom','source','type','start','end','V6','strand','V8','V9','gene_id','isoform','exon','name','nearestRef','classCode')
xlocRange = makeGRangesFromDataFrame(geneFeatures,keep.extra.columns=T)

sel = read.table('.../parallelSalinityDE.txt',stringsAsFactors=F)
sel = sel[,1]

extremeRange = xlocRange[which(xlocRange$gene_id %in% sel)]

normLRrange$xloc = NA
normLRrange$gene = NA
hits = findOverlaps(normLRrange,extremeRange)
normLRrange$xloc[queryHits(hits)] = extremeRange$gene_id[subjectHits(hits)]
normLRrange$gene[queryHits(hits)] = extremeRange$name[subjectHits(hits)]
normLRrangePar = subsetByOverlaps(normLRrange,extremeRange)

parNormLR = as.data.frame(mcols(normLRrangePar))
rownames(parNormLR) = names(normLRrangePar)[as.numeric(rownames(parNormLR))]

parNormLR[mapply(is.infinite,parNormLR)] = NA
parNormLR = na.omit(parNormLR)
#rownames(parNormLR) = paste(parNormLR[,13],rownames(parNormLR))
colnames(parNormLR) = c('3.5','3.5','3.5','3.5','35','0.2','0.2','0.2','0.2','35','35','xloc','name')

# pick one SNP per xloc
parNormLR = parNormLR[!duplicated(parNormLR$xloc),]

# plot
library("pheatmap")
library("RColorBrewer")
paletteLength <- 50
myColor <- colorRampPalette(c("turquoise", "black", "yellow"))(paletteLength)
myBreaks <- c(seq(min(parNormLR[,1:11]), 0, length.out=ceiling(paletteLength/2) + 1),seq(max(parNormLR[,1:11])/paletteLength, max(parNormLR[,1:11]), length.out=floor(paletteLength/2)))
#pdf('.../salinityDEaseHeatmap.pdf')
pheatmap(parNormLR[,1:11],col=myColor,breaks=myBreaks,border_color=NA,clustering_distance_rows='euclidean',clustering_distance_cols='euclidean')
#dev.off()


#############
# Correlation
countDist = cor(parNormLR[,1:11],method='spearman')
countDistMatrix <- as.matrix(countDist)
#rownames(countDistMatrix) <- design$sal
#colnames(countDistMatrix) <- design$sal
myBreaks <- c(-1, 0, length.out=ceiling(paletteLength/2) + 1),seq(max(countDistMatrix)/paletteLength, max(countDistMatrix), length.out=floor(paletteLength/2)))
paletteLength <- 50
colors <- colorRampPalette(c("turquoise", "black", "yellow"))(paletteLength)
pdf('.../parallelSalinityDEaseCorrelation.pdf')
pheatmap(countDistMatrix,
clustering_distance_rows='euclidean',
clustering_distance_cols='euclidean',
col=colors,
breaks=myBreaks)
dev.off()

##############
# candidates
aqp3a = 'chrXIII:10682178'
slc12a2 = 'chrXIV:13967752'
cands = c(slc12a2,aqp3a)

plotData = melt(data.frame(t(parNormLR[cands,1:11]),sal=c(rep('3.5',4),rep('0.2',4),rep('35',3))),gene=rep(cands,11))

plotData$sal = factor(plotData$sal,levels=c("0.2","3.5","35"))
plotData$variable = factor(plotData$variable, levels=c("chrUn.53967592","chrXIV.13967752","chrXIII.10682178","chrI:16782412"))

x11(h=5,w=2)
ggplot(plotData,aes(sal,value))+geom_boxplot()+theme(panel.background=element_rect(fill='white',colour='grey78'),legend.position="right",text=element_text(size=16,colour='black'))+facet_wrap(~variable,ncol=1,scales='free')+geom_hline(yintercept=0,colour='black',linetype='dotted')
ggsave('.../candidateSalinityASEBoxplotFDR1.ps')


#####
# PCA
library(ggbiplot)
pca = prcomp(t(parNormLR[,1:11]),scale=F)
percentVar <- round(100 * pca$sdev^2/sum(pca$sdev^2))

g = ggbiplot(pca, var.axes=F, choices=c(1,2), obs.scale=1, var.scale=1, groups=as.factor(colnames(parNormLR)[1:11]), ellipse=TRUE, ellipse.prob=0.95, alpha=0)
g = g + scale_color_manual(values=c('darkblue','darkred','grey','black','black'))
g = g + theme(axis.text=element_text(hjust=1,size=24,colour='black'),panel.background=element_rect(fill='white',colour='black'),panel.grid.major=element_line(colour='grey78'),text=element_text(size=22),legend.position="right")
g = g + geom_point(data=as.data.frame(pca$x), aes(PC1, PC2, color=colnames(parNormLR)[1:11]),size=8)
g
#ggsave('.../salinityPCAaseLogRatioSalinityDE1FDR.ps')

