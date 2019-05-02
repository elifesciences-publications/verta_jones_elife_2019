library(ggplot2)

################
# Composite PC #
################
PC2 = read.table('.../Tyne_Litc_PC2_loadings.txt')
PC5 = read.table('.../Tyne_Litc_PC5_loadings.txt')

PC = data.frame(x=0.145*PC2$x+0.063*-1*PC5$x)
rownames(PC) = rownames(PC2)
extremePCglobal = rownames(PC)[which(PC$x<quantile(PC$x,probs=c(0.01)) | PC$x>quantile(PC$x,probs=c(0.99)))]

#################
# Tyne-Litc CSS #
#################

# CSS peaks versus PC2 loadings -- all transcripts
library(GenomicRanges)

# xloc annotation
geneFeatures = read.table('.../mergedParsed.txt',stringsAsFactors=F,sep='\t',skip=1,fill=NA)
colnames(geneFeatures) = c('chrom','source','type','start','end','V6','strand','V8','V9','gene_id','isoform','exon','name','nearestRef','classCode')
xlocRange = makeGRangesFromDataFrame(geneFeatures,keep.extra.columns=T)

# css
css = read.table('.../Litc_Tyne_genomic_recalibrated_filtered_SNPs_CSS_w10000_s5000.txt',header=T,sep='\t')
# filter out windows with few SNPs
#css = css[which(css$tyne_mar_inds.txt.tyne_fw_inds.txt_nVariants >= 10 & css$litc_mar_inds.txt.litc_fw_inds.txt_nVariants >= 10 & css$tyne_mar_inds.txt.litc_mar_inds.txt_nVariants >= 10 & css$tyne_mar_inds.txt.litc_fw_inds.txt_nVariants >= 10 & css$tyne_fw_inds.txt.litc_mar_inds.txt_nVariants >= 10 & css$tyne_fw_inds.txt.litc_fw_inds.txt_nVariants >= 10),]
# transform into range
cssRange = makeGRangesFromDataFrame(data.frame(chrom=css[,'CHROM'],start=css[,'BIN_START'],end=css[,'BIN_END'],css=css[,'CSS']),keep.extra.columns=T)

# convert to zCSS
cssRange$zCSS = scale(cssRange$css)

# extreme CSS windows
cssExtreme = css[which(css$CSS >= quantile(css$CSS,.995)),]
cssExRange = makeGRangesFromDataFrame(data.frame(chrom=cssExtreme[,'CHROM'],start=cssExtreme[,'BIN_START'],end=cssExtreme[,'BIN_END'],css=cssExtreme[,'CSS']),keep.extra.columns=T)

# distance of XLOC to closest css peak
xlocRange$distToCssPeak = NA
hits = distanceToNearest(xlocRange,cssExRange)
xlocRange$distToCssPeak[queryHits(hits)] = mcols(hits)$distance

# are PC extremes (1%) closer to CSS peaks than other transcripts?
windows = seq(1,1000000,100000)
exPCxloc = xlocRange[which(xlocRange$gene_id %in% extremePCglobal)]
controlXloc = xlocRange[which(xlocRange$gene_id %in% extremePCglobal == F)]

distToCssPeak = data.frame()
for (i in 1:length(windows)){
    distToCssPeak[i,'globalExtremePC'] = length(unique(exPCxloc[which(exPCxloc$distToCssPeak < windows[i])]$gene_id))/length(unique(exPCxloc$gene_id))
    distToCssPeak[i,'control'] = length(unique(controlXloc[which(controlXloc$distToCssPeak < windows[i])]$gene_id))/length(unique(controlXloc$gene_id))
    distToCssPeak[i,'distance'] = windows[i]
}
rownames(distToCssPeak) = as.character(windows)

g = ggplot(distToCssPeak)
g = g + geom_line(aes(distance,control,color='grey'),linetype=2,size=2)
g = g + geom_line(aes(distance,globalExtremePC,color='blue'),size=2)
g = g + scale_color_manual(values=c('blue','darkgrey'))
g = g + theme(panel.background=element_rect(fill='white',colour='grey78'),panel.grid.major=element_line(colour='grey78'),legend.position="right",text=element_text(size=16))
g

# random expetation
distToCssPeakRandom = data.frame()
for (y in 1:1000){
    for (i in 1:length(windows)){
        sampleXloc = controlXloc[which(controlXloc$gene_id %in% sample(unique(controlXloc$gene_id),586)),]
        distToCssPeakRandom[i,y] = length(unique(sampleXloc[which(sampleXloc$distToCssPeak < windows[i])]$gene_id))/length(unique(sampleXloc$gene_id))
    }
}
rownames(distToCssPeakRandom) = as.character(windows)

quant95 = function(x){quantile(x,probs=c(0.05,0.95),na.rm=T)}
randomExp95IC = apply(distToCssPeakRandom,1,quant95)
plotRand = melt(randomExp95IC)

g = ggplot(distToCssPeak)
g = g + geom_line(aes(distance,control,color='grey'),size=2)
g = g + geom_line(data=plotRand[plotRand$Var1=='5%',],aes(Var2,value,color='grey'),linetype=2,size=2)
g = g + geom_line(data=plotRand[plotRand$Var1=='95%',],aes(Var2,value,color='grey'),linetype=2,size=2)
g = g + geom_line(aes(distance,globalExtremePC,color='blue'),size=2)
g = g + scale_color_manual(values=c('blue','darkgrey'))
g = g + theme(panel.background=element_rect(fill='white',colour='grey78'),panel.grid.major=element_line(colour='grey78'),legend.position="right",text=element_text(size=16))
g
ggsave('.../distToCssPeak995GlobalPCoutliers.ps')


####################
# Jones et al loci #
####################

# adaptive loci
adaptiveLoci = read.table('.../243_adaptive_regions.bed',stringsAsFactors=F,header=F)
colnames(adaptiveLoci) = c('chrom','start','end')
adaptiveRange = makeGRangesFromDataFrame(adaptiveLoci,keep.extra.columns=T)

# distance of XLOC to closest adaptive peak
xlocRange$distToAdaptive = NA
hits = distanceToNearest(xlocRange,adaptiveRange)
xlocRange$distToAdaptive[queryHits(hits)] = mcols(hits)$distance

# are PC2 extremes (1% - same sign Tyen & Litc) closer to adaptive peaks than other transcripts?
windows = seq(1,1000000,100000)
exPCxloc = xlocRange[which(xlocRange$gene_id %in% extremePCglobal)]
controlXloc = xlocRange[which(xlocRange$gene_id %in% extremePCglobal == F)]

# are PC extremes (1%) closer to Jones et al peaks than other transcripts?
distToAdaptive = data.frame()
for (i in 1:length(windows)){
    distToAdaptive[i,'globalExtremePC'] = length(unique(exPCxloc[which(exPCxloc$distToAdaptive < windows[i])]$gene_id))/length(unique(exPCxloc$gene_id))
    distToAdaptive[i,'control'] = length(unique(controlXloc[which(controlXloc$distToAdaptive < windows[i])]$gene_id))/length(unique(controlXloc$gene_id))
    distToAdaptive[i,'distance'] = windows[i]
}
rownames(distToAdaptive) = as.character(windows)

g = ggplot(distToAdaptive)
g = g + geom_line(aes(distance,control,color='grey'),linetype=2,size=2)
g = g + geom_line(aes(distance,globalExtremePC,color='blue'),size=2)
g = g + scale_color_manual(values=c('blue','darkgrey'))
g = g + theme(panel.background=element_rect(fill='white',colour='grey78'),panel.grid.major=element_line(colour='grey78'),legend.position="right",text=element_text(size=16))
g

# random expetation
distToAdaptiveRandom = data.frame()
for (y in 1:1000){
    for (i in 1:length(windows)){
        sampleXloc = controlXloc[which(controlXloc$gene_id %in% sample(unique(controlXloc$gene_id),586)),]
        distToAdaptiveRandom[i,y] = length(unique(sampleXloc[which(sampleXloc$distToAdaptive < windows[i])]$gene_id))/length(unique(sampleXloc$gene_id))
    }
}
rownames(distToAdaptiveRandom) = as.character(windows)

quant95 = function(x){quantile(x,probs=c(0.05,0.95))}
randomExp95IC = apply(distToAdaptiveRandom,1,quant95)
plotRand = melt(randomExp95IC)

g = ggplot(distToAdaptive)
g = g + geom_line(aes(distance,control,color='grey'),size=2)
g = g + geom_line(data=plotRand[plotRand$X1=='5%',],aes(X2,value,color='grey'),linetype=2,size=2)
g = g + geom_line(data=plotRand[plotRand$X1=='95%',],aes(X2,value,color='grey'),linetype=2,size=2)
g = g + geom_line(aes(distance,globalExtremePC,color='blue'),size=2)
g = g + scale_color_manual(values=c('blue','darkgrey'))
g = g + theme(panel.background=element_rect(fill='white',colour='grey78'),panel.grid.major=element_line(colour='grey78'),legend.position="right",text=element_text(size=16))
g
ggsave('.../distToAdaptiveJonesPCoutliers.ps')

