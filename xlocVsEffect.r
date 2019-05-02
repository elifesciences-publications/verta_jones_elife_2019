#####################
# add transcript ID #
# summarize effects #
#####################
library(GenomicRanges)

# xloc annotation
geneFeatures = read.table('.../assembly/merged_assembly_STAR_UCSCgtf/mergedParsed.txt',stringsAsFactors=F,sep='\t',skip=1,fill=NA)
colnames(geneFeatures) = c('chrom','source','type','start','end','V6','strand','V8','V9','gene_id','isoform','exon','name','nearestRef','classCode')
xlocRange = makeGRangesFromDataFrame(geneFeatures,keep.extra.columns=T)

# load effect table
ind = 'c214_F1_1_M'

effectTable = read.table(paste('.../effectTable_',ind,'_FDR10_minCov10_aseReadCounts_refCorrected_F1normalized_genomicMask_STAR_duprem.txt',sep=''),header=T,stringsAsFactors=F)

eRange = makeGRangesFromDataFrame(data.frame(chrom=effectTable[,1],start=effectTable[,2],end=effectTable[,2]+1,effectTable[,3:ncol(effectTable)]),keep.extra.columns=T)

eRange$xloc = NA
hits = findOverlaps(xlocRange,eRange)
idx = queryHits(hits)
queryXloc = as.character(xlocRange$gene_id[idx])
recipIdx = subjectHits(hits)
eRange$xloc[recipIdx] = queryXloc

eRange$isoform = NA
hits = findOverlaps(xlocRange,eRange)
idx = queryHits(hits)
queryIso = as.character(xlocRange$isoform[idx])
recipIdx = subjectHits(hits)
eRange$isoform[recipIdx] = queryIso

eRange$exon = NA
hits = findOverlaps(xlocRange,eRange)
idx = queryHits(hits)
queryExon = as.character(xlocRange$exon[idx])
recipIdx = subjectHits(hits)
eRange$exon[recipIdx] = queryExon

eRange$name = NA
hits = findOverlaps(xlocRange,eRange)
idx = queryHits(hits)
queryName = as.character(xlocRange$name[idx])
recipIdx = subjectHits(hits)
eRange$name[recipIdx] = queryName

effectTable = as.data.frame(eRange)
head(effectTable)

# save the table with the xloc info
write.table(effectTable,paste('.../effectTable_',ind,'_FDR10_minCov10_aseReadCounts_refCorrected_F1normalized_genomicMask_STAR_duprem.txt',sep=''))
rm(effectTable)
rm(ind)





