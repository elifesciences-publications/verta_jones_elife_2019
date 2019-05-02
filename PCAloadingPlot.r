###########################
# plot PC over chromosomes
library(GenomicRanges)
library(zoo)
library(qqman)
library(GenomicFeatures)
library(ggrepel)

geneFeatures = read.table('.../assembly/merged_assembly_STAR_UCSCgtf/mergedParsed.txt',stringsAsFactors=F,sep='\t',skip=1)
colnames(geneFeatures) = c('chrom','source','type','start','end','V6','strand','V8','V9','gene_id','tcons','exon','name','oID','ensGene')
xlocRange = makeGRangesFromDataFrame(geneFeatures,keep.extra.columns=T)

PC2 = read.table('.../ASE/UCSCgtf/Tyne_Litc_PC2_loadings.txt')
PC5 = read.table('.../ASE/UCSCgtf/Tyne_Litc_PC5_loadings.txt')

PC = data.frame(x=0.145*PC2$x+0.063*-1*PC5$x)
rownames(PC) = rownames(PC2)
extremePCglobal = rownames(PC)[which(PC$x<quantile(PC$x,probs=c(0.01)) | PC$x>quantile(PC$x,probs=c(0.99)))]
PC$zPC = scale(PC$x)

xlocRange$PC = PC[xlocRange$gene_id,'zPC']

PCplot = data.frame(chrom=seqnames(xlocRange),pos=start(xlocRange),PC=xlocRange$PC,xloc=xlocRange$gene_id)
PCplot = PCplot[which(PCplot$chrom != 'MT'),]
PCplot$plotPos = as.numeric(PCplot$pos)
PCplot$chromNum = as.numeric(PCplot$chrom)
PCplot = PCplot[order(PCplot$chromNum),]
PCplot$bg = PCplot$chromNum
PCplot$bg[which(PCplot$chromNum %in% seq(1,21,2))] = 'grey10'
PCplot$bg[which(PCplot$chromNum %in% seq(2,20,2))] = 'grey80'

# reduce range to one line per gene
PCplotRed=data.frame()
for (i in unique(PCplot$xloc)){
    xloc = PCplot[which(PCplot$xloc == i),]
    s = min(xloc$pos)
    e = max(xloc$pos)
    PCplotRed[i,'chrom'] = unique(xloc$chrom)
    PCplotRed[i,'pos'] = (s+e)/2
    PCplotRed[i,'PC'] = unique(xloc$PC)
    PCplotRed[i,'bg'] = unique(xloc$bg)
}
PCplotRed$chromNum = as.numeric(PCplotRed$chrom)

write.table(PCplotRed,'.../ASE/UCSCgtf/compositePCloadingsXlocLitcTyne.txt')

#---------------

PCplotRed = read.table('.../ASE/UCSCgtf/compositePCloadingsXlocLitcTyne.txt')
PCplot = PCplotRed
PCplot$zPC = PC[rownames(PCplot),'zPC']
PCplot$chrom = factor(PCplot$chrom,levels=c('chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrVI', 'chrVII', 'chrVIII', 'chrIX','chrX', 'chrXI','chrXII', 'chrXIII', 'chrXIV', 'chrXV', 'chrXVI', 'chrXVII', 'chrXVIII', 'chrXIX', 'chrXX', 'chrXXI'))

PCplot$plotPos = 0
for (i in PCplot$chromNum){
    if (i == 1){
        PCplot[which(PCplot$chromNum == i),'plotPos'] = PCplot[which(PCplot$chromNum == i),'pos']
    }
    else {
        lastPos = max(PCplot[which(PCplot$chromNum == i-1),'plotPos'])
        PCplot[which(PCplot$chromNum == i),'plotPos'] = PCplot[which(PCplot$chromNum == i),'pos'] + lastPos
    }
}

candPlot = PCplot[which(rownames(PCplot) %in% c('XLOC_000759','XLOC_027088','XLOC_028346','XLOC_019499','XLOC_020689','XLOC_034273','XLOC_000171','XLOC_003523','XLOC_025119','XLOC_034131')),]
candPlot$name = NA
for (i in 1:nrow(candPlot)){
    candPlot$name[i] = unique(geneFeatures[which(geneFeatures$gene_id == rownames(candPlot)[i]),'gene_name'])
}
candPlot['XLOC_034131','name'] = 'XLOC_034131'
candPlot$chrom = factor(candPlot$chrom,levels=c('chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrVI', 'chrVII', 'chrVIII', 'chrIX','chrX', 'chrXI','chrXII', 'chrXIII', 'chrXIV', 'chrXV', 'chrXVI', 'chrXVII', 'chrXVIII', 'chrXIX', 'chrXX', 'chrXXI'))
lowPC = min(PC[extremePCglobal,'zPC'][which(PC[extremePCglobal,1] > 0)])
highPC = max(PC[extremePCglobal,'zPC'][which(PC[extremePCglobal,1] < 0)])

#postscript('.../ASE/UCSCgtf/zPCloadingsAllChrWithAnnotation.ps',height=3,width=10)
x11(height=3,width=10)
ggplot(PCplot[which(PCplot$chromNum<=21),])+
geom_point(aes(as.numeric(plotPos),zPC,colour=bg))+
scale_colour_manual(values=c("grey20","grey60"))+
theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.text.y=element_text(hjust=1,size=24,colour='black'),panel.background=element_rect(fill='white',colour='black'),panel.grid.major=element_blank(),text=element_text(size=22),legend.position="none")+
geom_hline(yintercept=c(lowPC,highPC),linetype='dotted')+
geom_label_repel(data=candPlot,aes(as.numeric(plotPos),zPC,label=name),segment.color='red',nudge_y=-4, size = 5)
#dev.off()

#--------------------
# chromosomes separately

lowPC = min(PC[extremePCglobal,'zPC'][which(PC[extremePCglobal,1] > 0)])
highPC = max(PC[extremePCglobal,'zPC'][which(PC[extremePCglobal,1] < 0)])

invXI=data.frame(chr='chrXI',chromNum=11,start=5440949,end=5852218)
invI = data.frame(chr='chrI',chromNum=1,start=21497667,end=21936224)
invXXI=data.frame(chr='chrXXI',chromNum=21,start=5789748,end=7482685)

#x11(height=10,width=15)
ggplot(PCplot[which(PCplot$chromNum<=21),])+
geom_point(aes(as.numeric(pos),zPC,colour=bg),size=0.5)+
scale_colour_manual(values=c("grey20","grey60"))+
theme(axis.title.x=element_blank(),axis.text.y=element_text(hjust=1,size=24,colour='black'),panel.background=element_rect(fill='white',colour='black'),panel.grid.major=element_blank(),text=element_text(size=22),legend.position="none")+
geom_hline(yintercept=c(lowPC,highPC),linetype='dotted')+
#geom_label_repel(aes(as.numeric(plotPos),zPC,label=geneFeatures[rownames(PCplot),"gene_name"]),segment.color='red',nudge_y=-4, size = 5)+
geom_segment(data=invXI,aes(y=0,yend=0,x=start,xend=end),colour='red')+
geom_segment(data=invI,aes(y=0,yend=0,x=start,xend=end),colour='red')+
geom_segment(data=invXXI,aes(y=0,yend=0,x=start,xend=end),colour='red')+
facet_wrap(~chromNum)+
ylim(c(-3,3))
ggsave('.../ASE/UCSCgtf/zPCloadingsSeparateChr.ps')

#---------------------
# plot labels for 10 most extreme compPC transcripts

compPCoutliers = c(rownames(PC[order(PC$zPC),])[1:30],rownames(PC[order(PC$zPC,decreasing=T),])[1:30])

candPlot = PCplot[which(rownames(PCplot) %in% compPCoutliers),]
candPlot$name = NA
for (i in 1:nrow(candPlot)){
    if (is.na(unique(geneFeatures[which(geneFeatures$gene_id == rownames(candPlot)[i]),'gene_name'])) == F){
        candPlot$name[i] = unique(geneFeatures[which(geneFeatures$gene_id == rownames(candPlot)[i]),'gene_name'])
    }
    else {
        candPlot$name[i] = rownames(candPlot)[i]
    }
}
candPlot = candPlot[which(candPlot$chromNum<=21),]
candPlot$chrom = factor(candPlot$chrom,levels=c('chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrVI', 'chrVII', 'chrVIII', 'chrIX','chrX', 'chrXI','chrXII', 'chrXIII', 'chrXIV', 'chrXV', 'chrXVI', 'chrXVII', 'chrXVIII', 'chrXIX', 'chrXX', 'chrXXI'))
lowPC = min(PC[extremePCglobal,'zPC'][which(PC[extremePCglobal,1] > 0)])
highPC = max(PC[extremePCglobal,'zPC'][which(PC[extremePCglobal,1] < 0)])

#x11(height=10,width=15)
ggplot(PCplot[which(PCplot$chromNum<=21),])+
geom_point(aes(as.numeric(plotPos),zPC,colour=bg))+
scale_colour_manual(values=c("grey20","grey60"))+
theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.text.y=element_text(hjust=1,size=24,colour='black'),panel.background=element_rect(fill='white',colour='black'),panel.grid.major=element_blank(),text=element_text(size=22),legend.position="none")+
geom_hline(yintercept=c(lowPC,highPC),linetype='dotted')+
geom_label_repel(data=candPlot,aes(as.numeric(plotPos),zPC,label=name),segment.color='red',nudge_y=-4, size = 5)


#-----------------------
# plot extreme compPC AND candidate osmoregulatory genes, informative labels

osmoGenes = c('XLOC_000759','XLOC_027088','XLOC_028346','XLOC_019499','XLOC_020689','XLOC_034273','XLOC_000171','XLOC_003523','XLOC_025119','XLOC_034131')
compPCoutliers = c("XLOC_008723", "XLOC_017483", "XLOC_022789", "XLOC_010811", "XLOC_011575", "XLOC_019359", "XLOC_010139", "XLOC_030091", "XLOC_022350", "XLOC_018499", "XLOC_011672", "XLOC_031542", "XLOC_000171", "XLOC_022351", "XLOC_016212", "XLOC_011673", "XLOC_008821", "XLOC_023191","XLOC_001323","XLOC_000420")

cands = c(osmoGenes,compPCoutliers)

candPlot = PCplot[which(rownames(PCplot) %in% cands),]
candPlot$name = NA
for (i in 1:nrow(candPlot)){
    if (is.na(unique(geneFeatures[which(geneFeatures$gene_id == rownames(candPlot)[i]),'gene_name'])) == F){
        candPlot$name[i] = unique(geneFeatures[which(geneFeatures$gene_id == rownames(candPlot)[i]),'gene_name'])
    }
    else {
        candPlot$name[i] = rownames(candPlot)[i]
    }
}
candPlot = candPlot[which(candPlot$chromNum<=21),]
candPlot$chrom = factor(candPlot$chrom,levels=c('chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrVI', 'chrVII', 'chrVIII', 'chrIX','chrX', 'chrXI','chrXII', 'chrXIII', 'chrXIV', 'chrXV', 'chrXVI', 'chrXVII', 'chrXVIII', 'chrXIX', 'chrXX', 'chrXXI'))
lowPC = min(PC[extremePCglobal,'zPC'][which(PC[extremePCglobal,1] > 0)])
highPC = max(PC[extremePCglobal,'zPC'][which(PC[extremePCglobal,1] < 0)])

# change names to more informative
candPlot[which(candPlot$name == 'ENSGACG00000018955'),'name'] = c('slc12a10')
candPlot[which(candPlot$name == 'ENSGACG00000012554'),'name'] = c('IFI44L 1')
candPlot[which(candPlot$name == 'ENSGACG00000012548'),'name'] = c('IFI44L 2')
candPlot[which(candPlot$name == 'ENSGACG00000010910'),'name'] = c('HBB-like 1')
candPlot[which(candPlot$name == 'ENSGACG00000013918'),'name'] = c('HBB-like 2')
candPlot[which(candPlot$name == 'ENSGACG00000002627'),'name'] = c('GTPase-like 1')
candPlot[which(candPlot$name == 'si:ch211-110e4.2 (5 of 6)'),'name'] = c('GTPase-like 2')
candPlot[which(candPlot$name == 'si:ch211-141h20.6'),'name'] = c('SAM9L')
candPlot[which(candPlot$name == 'ENSGACG00000018668'),'name'] = c('FGB-like')
candPlot[which(candPlot$name == 'XLOC_034131'),'name'] = c('Ig-like')
candPlot[which(candPlot$name == 'ENSGACG00000022540'),'name'] = c('slc12-like')

#x11(height=3,width=10)
postscript('.../ASE/UCSCgtf/zPCloadingsAllChrWithAnnotation_v2.ps',height=3,width=10)
ggplot(PCplot[which(PCplot$chromNum<=21),])+
geom_point(aes(as.numeric(plotPos),zPC,colour=bg))+
scale_colour_manual(values=c("grey20","grey60"))+
theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.text.y=element_text(hjust=1,size=24,colour='black'),panel.background=element_rect(fill='white',colour='black'),panel.grid.major=element_blank(),text=element_text(size=22),legend.position="none")+
geom_hline(yintercept=c(lowPC,highPC),linetype='dotted')+
geom_label_repel(data=candPlot,aes(as.numeric(plotPos),zPC,label=name),segment.color='red', size = 5)+
geom_point(data=candPlot,aes(as.numeric(plotPos),zPC),colour='red')
dev.off()









