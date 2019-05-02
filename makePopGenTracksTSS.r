library(GenomicFeatures)
library(GenomicRanges)

# xloc annotation
geneFeatures = read.table('.../assembly/merged_assembly_STAR_UCSCgtf/mergedParsed.txt',stringsAsFactors=F,sep='\t',skip=1,fill=NA)
colnames(geneFeatures) = c('chrom','source','type','start','end','V6','strand','V8','V9','gene_id','tcons','exon','name','nearestRef','classCode')
xlocRange = makeGRangesFromDataFrame(geneFeatures,keep.extra.columns=T)

# load transcript GTF track
txdb = loadDb(".../assembly/merged_assembly_STAR_UCSCgtf/mergedNoUknownStrand.sqlite")

###############################################################
# create a range with bins centered around TSS,-250 kb, +250 kb #
###############################################################

# target genes (all xloc)
target = data.frame(chr=NA,start=NA,end=NA,strand=NA,xloc=NA)
for (i in na.omit(unique(xlocRange$gene_id))){
    rg = xlocRange[which(xlocRange$gene_id == i)]
    target[i,'chr'] = as.character(unique(seqnames(rg)))
    target[i,'start'] = min(start(rg))-250000
    target[i,'end'] = min(start(rg))+250000
    target[i,'strand'] = unique(as.character(strand(rg)))
    target[i,'xloc'] = i
}
target = target[-1,]
write.table(target,'.../xlocTssCordinates500kb.txt')
target = read.table('.../xlocTssCordinates500kb.txt')
targetRange = makeGRangesFromDataFrame(target,keep.extra.columns=T)

###################
# window function #
###################

windows = function(gr,size=1000){
    starts = seq(min(start(gr)),max(end(gr))-size,size)
    ends = seq(min(start(gr))+size-1,max(end(gr)),size)
    ends[length(ends)-1] = max(end(gr))
    chrom = rep(seqnames(gr),length(starts))
    str = rep(strand(gr),length(starts))
    return(makeGRangesFromDataFrame(data.frame(chrom,start=starts,end=ends,strand=str)))
}

####################
# ----- DATA ----- #
####################

#------#
# LITC
#------#

# Fst
litcFst = read.table('.../LITC_genomic_recalibrated_filtered_SNPs_maxMissing4.weir.fst',sep='\t',header=T,stringsAsFactors=F)
fstRange = makeGRangesFromDataFrame(data.frame(chrom=litcFst[,'CHROM'],start=litcFst[,'POS'],end=litcFst[,'POS']+1,Fst=litcFst[,'WEIR_AND_COCKERHAM_FST']),keep.extra.columns=T)

# Pi marine
litcPiMar = read.table('.../LITC_MAR_genomic_recalibrated_filtered_SNPs_maxMissing2.sites.pi',sep='\t',header=T,stringsAsFactors=F)
PiMarRange = makeGRangesFromDataFrame(data.frame(chrom=litcPiMar[,'CHROM'],start=litcPiMar[,'POS'],end=litcPiMar[,'POS']+1,Pi=litcPiMar[,'PI']),keep.extra.columns=T)

#  Pi Freshwater
litcPiFw = read.table('.../LITC_FW_genomic_recalibrated_filtered_SNPs_maxMissing2.sites.pi',sep='\t',header=T,stringsAsFactors=F)
PiFwRange = makeGRangesFromDataFrame(data.frame(chrom=litcPiFw[,'CHROM'],start=litcPiFw[,'POS'],end=litcPiFw[,'POS']+1,Pi=litcPiFw[,'PI']),keep.extra.columns=T)

# css
css = read.table('.../Litc_Tyne_genomic_recalibrated_filtered_SNPs_CSS_w1_s1_informative.txt',header=T,sep='\t')
cssRange = makeGRangesFromDataFrame(data.frame(chrom=css[,'CHROM'],start=css[,'BIN_START'],end=css[,'BIN_END'],css=css[,'CSS']),keep.extra.columns=T)


# calculate average values for windows - ALL XLOC
targetFst = as.data.frame(matrix(ncol=500,nrow=length(unique(targetRange$xloc)),dimnames=list(unique(targetRange$xloc),c(1:500))))
targetPiMar = as.data.frame(matrix(ncol=500,nrow=length(unique(targetRange$xloc)),dimnames=list(unique(targetRange$xloc),c(1:500))))
targetPiFw = as.data.frame(matrix(ncol=500,nrow=length(unique(targetRange$xloc)),dimnames=list(unique(targetRange$xloc),c(1:500))))
targetCSS = as.data.frame(matrix(ncol=500,nrow=length(unique(targetRange$xloc)),dimnames=list(unique(targetRange$xloc),c(1:500))))

for (xloc in targetRange$xloc){
    rg = targetRange[which(targetRange$xloc == xloc)]

    rgFst = fstRange[subjectHits(findOverlaps(rg,fstRange))]
    rgPiMar = PiMarRange[subjectHits(findOverlaps(rg,PiMarRange))]
    rgPiFw = PiFwRange[subjectHits(findOverlaps(rg,PiFwRange))]
    rgCSS = cssRange[subjectHits(findOverlaps(rg,cssRange))]

    w = windows(rg,1000)
    overlapsFst = findOverlaps(rgFst, w, select="arbitrary")
    overlapsPiMar = findOverlaps(rgPiMar, w, select="arbitrary")
    overlapsPiFw = findOverlaps(rgPiFw, w, select="arbitrary")
    overlapsCSS = findOverlaps(rgCSS, w, select="arbitrary")

    if (unique(strand(w)) == '+'){
        averagedFst = aggregate(rgFst$Fst,list(overlapsFst),mean,na.rm=T)
        averagedPiMar = aggregate(rgPiMar$Pi,list(overlapsPiMar),mean,na.rm=T)
        averagedPiFw = aggregate(rgPiFw$Pi,list(overlapsPiFw),mean,na.rm=T)
        averagedCSS = aggregate(rgCSS$css,list(overlapsCSS),mean,na.rm=T)
        targetFst[xloc,averagedFst$Group.1] = averagedFst$x
        targetPiMar[xloc,averagedPiMar$Group.1] = averagedPiMar$x
        targetPiFw[xloc,averagedPiFw$Group.1] = averagedPiFw$x
        targetCSS[xloc,averagedCSS$Group.1] = averagedCSS$x
    }

    if(unique(strand(w)) == '-'){
        averagedFst = aggregate(rgFst$Fst,rev(list(overlapsFst)),mean,na.rm=T)
        averagedPiMar = aggregate(rgPiMar$Pi,rev(list(overlapsPiMar)),mean,na.rm=T)
        averagedPiFw = aggregate(rgPiFw$Pi,rev(list(overlapsPiFw)),mean,na.rm=T)
        averagedCSS = aggregate(rgCSS$css,rev(list(overlapsCSS)),mean,na.rm=T)
        targetFst[xloc,averagedFst$Group.1] = averagedFst$x
        targetPiMar[xloc,averagedPiMar$Group.1] = averagedPiMar$x
        targetPiFw[xloc,averagedPiFw$Group.1] = averagedPiFw$x
        targetCSS[xloc,averagedCSS$Group.1] = averagedCSS$x
    }
}
write.table(targetFst,'.../xlocLitcFstTssmatrix500kb.txt')
write.table(targetPiMar,'.../xlocLitcPiMarTssmatrix500kb.txt')
write.table(targetPiFw,'.../xlocLitcPiFwTssmatrix500kb.txt')
write.table(targetCSS,'.../xlocLitcCssTssmatrix500kb.txt')



#------#
# TYNE
#------#

# Fst
tyneFst = read.table('.../Tyne_genomic_recalibrated_filtered_SNPs_maxMissing4.weir.fst',sep='\t',header=T,stringsAsFactors=F)
fstRange = makeGRangesFromDataFrame(data.frame(chrom=tyneFst[,'CHROM'],start=tyneFst[,'POS'],end=tyneFst[,'POS']+1,Fst=tyneFst[,'WEIR_AND_COCKERHAM_FST']),keep.extra.columns=T)

# Pi Marine
tynePiMar = read.table('.../Tyne_MAR_genomic_recalibrated_filtered_SNPs_maxMissing2.sites.pi',sep='\t',header=T,stringsAsFactors=F)
PiMarRange = makeGRangesFromDataFrame(data.frame(chrom=tynePiMar[,'CHROM'],start=tynePiMar[,'POS'],end=tynePiMar[,'POS']+1,Pi=tynePiMar[,'PI']),keep.extra.columns=T)

# Pi Freshwater
tynePiFw = read.table('.../Tyne_FW_genomic_recalibrated_filtered_SNPs_maxMissing2.sites.pi',sep='\t',header=T,stringsAsFactors=F)
PiFwRange = makeGRangesFromDataFrame(data.frame(chrom=tynePiFw[,'CHROM'],start=tynePiFw[,'POS'],end=tynePiFw[,'POS']+1,Pi=tynePiFw[,'PI']),keep.extra.columns=T)

# calculate average values for windows - ALL XLOC
targetFst = as.data.frame(matrix(ncol=500,nrow=length(unique(targetRange$xloc)),dimnames=list(unique(targetRange$xloc),c(1:500))))
targetPiMar = as.data.frame(matrix(ncol=500,nrow=length(unique(targetRange$xloc)),dimnames=list(unique(targetRange$xloc),c(1:500))))
targetPiFw = as.data.frame(matrix(ncol=500,nrow=length(unique(targetRange$xloc)),dimnames=list(unique(targetRange$xloc),c(1:500))))
targetCSS = as.data.frame(matrix(ncol=500,nrow=length(unique(targetRange$xloc)),dimnames=list(unique(targetRange$xloc),c(1:500))))

for (xloc in targetRange$xloc){
    rg = targetRange[which(targetRange$xloc == xloc)]

    rgFst = fstRange[subjectHits(findOverlaps(rg,fstRange))]
    rgPiMar = PiMarRange[subjectHits(findOverlaps(rg,PiMarRange))]
    rgPiFw = PiFwRange[subjectHits(findOverlaps(rg,PiFwRange))]
    rgCSS = cssRange[subjectHits(findOverlaps(rg,cssRange))]

    w = windows(rg,1000)
    overlapsFst = findOverlaps(rgFst, w, select="arbitrary")
    overlapsPiMar = findOverlaps(rgPiMar, w, select="arbitrary")
    overlapsPiFw = findOverlaps(rgPiFw, w, select="arbitrary")
   overlapsCSS = findOverlaps(rgCSS, w, select="arbitrary")

    if (unique(strand(w)) == '+'){
        averagedFst = aggregate(rgFst$Fst,list(overlapsFst),mean,na.rm=T)
        averagedPiMar = aggregate(rgPiMar$Pi,list(overlapsPiMar),mean,na.rm=T)
        averagedPiFw = aggregate(rgPiFw$Pi,list(overlapsPiFw),mean,na.rm=T)
        averagedCSS = aggregate(rgCSS$css,list(overlapsCSS),mean,na.rm=T)
        targetFst[xloc,averagedFst$Group.1] = averagedFst$x
        targetPiMar[xloc,averagedPiMar$Group.1] = averagedPiMar$x
        targetPiFw[xloc,averagedPiFw$Group.1] = averagedPiFw$x
        targetCSS[xloc,averagedCSS$Group.1] = averagedCSS$x
    }

    if(unique(strand(w)) == '-'){
        averagedFst = aggregate(rgFst$Fst,rev(list(overlapsFst)),mean,na.rm=T)
        averagedPiMar = aggregate(rgPiMar$Pi,rev(list(overlapsPiMar)),mean,na.rm=T)
        averagedPiFw = aggregate(rgPiFw$Pi,rev(list(overlapsPiFw)),mean,na.rm=T)
        averagedCSS = aggregate(rgCSS$css,rev(list(overlapsCSS)),mean,na.rm=T)
        targetFst[xloc,averagedFst$Group.1] = averagedFst$x
        targetPiMar[xloc,averagedPiMar$Group.1] = averagedPiMar$x
        targetPiFw[xloc,averagedPiFw$Group.1] = averagedPiFw$x
        targetCSS[xloc,averagedCSS$Group.1] = averagedCSS$x
    }
}
write.table(targetFst,'.../xlocTyneFstTssmatrix500kb.txt')
write.table(targetPiMar,'.../xlocTynePiMarTssmatrix500kb.txt')
write.table(targetPiFw,'.../xlocTynePiFwTssmatrix500kb.txt')
write.table(targetCSS,'.../xlocTyneCssTssmatrix500kb.txt')


