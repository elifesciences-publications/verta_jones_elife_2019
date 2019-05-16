library(gplots)
library(ggplot2)
library(gdata)
library(reshape)
library(gtools)
library(RColorBrewer)


#####################################################################
# load effect tables (combined F1's, best effect per xloc in each F1)
# reorder classes
# calculate summary freqs for effect classes
effects172 = read.table('AseReadCountsEffects172_bestQvalue_cisParentTrans_STAR_duprem_nosexFDR10.txt',header=T,stringsAsFactors=F)
effects172$colClass = factor(effects172$colClass,levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous'))

effects363 = read.table('AseReadCountsEffects363_bestQvalue_cisParentTrans_STAR_duprem_nosexFDR10.txt',header=T,stringsAsFactors=F)
effects363$colClass = factor(effects363$colClass,levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous'))

###########################################
# Venn diagrams of effect sharing in F1's #
###########################################
library(colorfulVennPlot)

# ==> consider only xloc that are tested for effects in all F1's
tt = venn(list(
effects363[which(effects363$source=='div01F'),'xloc'],
effects363[which(effects363$source=='div01M'),'xloc'],
effects363[which(effects363$source=='div02F'),'xloc'],
effects363[which(effects363$source=='div02M'),'xloc']
),show.plot=F)

tested = attr(tt,"intersections")$"A:B:C:D"

# cis
cis = venn(list(
effects363[which(effects363$xloc %in% tested & effects363$colClass=='cis' & effects363$source=='div01F'),'xloc'],
effects363[which(effects363$xloc %in% tested & effects363$colClass=='cis' & effects363$source=='div01M'),'xloc'],
effects363[which(effects363$xloc %in% tested & effects363$colClass=='cis' & effects363$source=='div02F'),'xloc'],
effects363[which(effects363$xloc %in% tested & effects363$colClass=='cis' & effects363$source=='div02M'),'xloc']
),show.plot=F)

rgPal <- colorRampPalette(c('black','green'))
col = rgPal(10)[cut(unlist(lapply(attr(cis,"intersections"),length)),breaks=10)]
#pdf('bestEffects363cisEffectCisParentTransSharing.pdf')
plotVenn4d(round(unlist(lapply(attr(cis,"intersections"),length))/sum(unlist(lapply(attr(cis,"intersections"),length)))*100),Colors=col)
#dev.off()

# cis + trans
cisPtrans = venn(list(
effects363[which(effects363$xloc %in% tested & effects363$colClass=='cis+trans' & effects363$source=='div01F'),'xloc'],
effects363[which(effects363$xloc %in% tested & effects363$colClass=='cis+trans' & effects363$source=='div01M'),'xloc'],
effects363[which(effects363$xloc %in% tested & effects363$colClass=='cis+trans' & effects363$source=='div02F'),'xloc'],
effects363[which(effects363$xloc %in% tested & effects363$colClass=='cis+trans' & effects363$source=='div02M'),'xloc']
),show.plot=F)

rgPal <- colorRampPalette(c('black','green'))
col = rgPal(10)[cut(unlist(lapply(attr(cisPtrans,"intersections"),length)),breaks=10)]
#pdf('bestEffects363cisPtransEffectCisParentTransSharing.pdf')
plotVenn4d(round(unlist(lapply(attr(cisPtrans,"intersections"),length))/sum(unlist(lapply(attr(cisPtrans,"intersections"),length)))*100))
#dev.off()


# trans
trans = venn(list(
effects363[which(effects363$xloc %in% tested & effects363$colClass=='trans' & effects363$source=='div01F'),'xloc'],
effects363[which(effects363$xloc %in% tested & effects363$colClass=='trans' & effects363$source=='div01M'),'xloc'],
effects363[which(effects363$xloc %in% tested & effects363$colClass=='trans' & effects363$source=='div02F'),'xloc'],
effects363[which(effects363$xloc %in% tested & effects363$colClass=='trans' & effects363$source=='div02M'),'xloc']
),show.plot=F)

rgPal <- colorRampPalette(c('black','green'))
col = rgPal(10)[cut(unlist(lapply(attr(trans,"intersections"),length)),breaks=10)]
pdf('bestEffects363transEffectCisParentTransSharing.pdf')
plotVenn4d(round(unlist(lapply(attr(trans,"intersections"),length))/sum(unlist(lapply(attr(trans,"intersections"),length)))*100),Colors=col)
dev.off()

# cis - trans
cisMtrans = venn(list(
effects363[which(effects363$xloc %in% tested & effects363$colClass=='cis-trans' & effects363$source=='div01F'),'xloc'],
effects363[which(effects363$xloc %in% tested & effects363$colClass=='cis-trans' & effects363$source=='div01M'),'xloc'],
effects363[which(effects363$xloc %in% tested & effects363$colClass=='cis-trans' & effects363$source=='div02F'),'xloc'],
effects363[which(effects363$xloc %in% tested & effects363$colClass=='cis-trans' & effects363$source=='div02M'),'xloc']
),show.plot=F)

rgPal <- colorRampPalette(c('black','green'))
col = rgPal(10)[cut(unlist(lapply(attr(cisMtrans,"intersections"),length)),breaks=10)]
#pdf('bestEffects363cisMtransEffectCisParentTransSharing.pdf')
plotVenn4d(round(unlist(lapply(attr(cisMtrans,"intersections"),length))/sum(unlist(lapply(attr(cisMtrans,"intersections"),length)))*100))
#dev.off()

# compensatory
compensatory = venn(list(
effects363[which(effects363$xloc %in% tested & effects363$colClass=='compensatory' & effects363$source=='div01F'),'xloc'],
effects363[which(effects363$xloc %in% tested & effects363$colClass=='compensatory' & effects363$source=='div01M'),'xloc'],
effects363[which(effects363$xloc %in% tested & effects363$colClass=='compensatory' & effects363$source=='div02F'),'xloc'],
effects363[which(effects363$xloc %in% tested & effects363$colClass=='compensatory' & effects363$source=='div02M'),'xloc']
),show.plot=F)

rgPal <- colorRampPalette(c('black','green'))
col = rgPal(10)[cut(unlist(lapply(attr(compensatory,"intersections"),length)),breaks=10)]
pdf('bestEffects363compensatoryEffectCisParentTransSharing.pdf')
plotVenn4d(round(unlist(lapply(attr(compensatory,"intersections"),length))/sum(unlist(lapply(attr(compensatory,"intersections"),length)))*100),Colors=col)
dev.off()

################################################################
# Tyne
# ==> consider only xloc that are tested for effects in all F1's
tt = venn(list(
effects172[which(effects172$source=='div20F'),'xloc'],
effects172[which(effects172$source=='div20M'),'xloc'],
effects172[which(effects172$source=='div04F'),'xloc'],
effects172[which(effects172$source=='div04M'),'xloc']
),show.plot=F)

tested = attr(tt,"intersections")$"A:B:C:D"

# cis
cis = venn(list(
effects172[which(effects172$xloc %in% tested & effects172$colClass=='cis' & effects172$source=='div20F'),'xloc'],
effects172[which(effects172$xloc %in% tested & effects172$colClass=='cis' & effects172$source=='div20M'),'xloc'],
effects172[which(effects172$xloc %in% tested & effects172$colClass=='cis' & effects172$source=='div04F'),'xloc'],
effects172[which(effects172$xloc %in% tested & effects172$colClass=='cis' & effects172$source=='div04M'),'xloc']
),show.plot=F)

rgPal <- colorRampPalette(c('black','green'))
col = rgPal(10)[cut(unlist(lapply(attr(cis,"intersections"),length)),breaks=10)]
pdf('bestEffects172cisEffectCisParentTransSharing.pdf')
plotVenn4d(round(unlist(lapply(attr(cis,"intersections"),length))/sum(unlist(lapply(attr(cis,"intersections"),length)))*100),Colors=col)
dev.off()


# trans
trans = venn(list(
effects172[which(effects172$xloc %in% tested & effects172$colClass=='trans' & effects172$source=='div20F'),'xloc'],
effects172[which(effects172$xloc %in% tested & effects172$colClass=='trans' & effects172$source=='div20M'),'xloc'],
effects172[which(effects172$xloc %in% tested & effects172$colClass=='trans' & effects172$source=='div04F'),'xloc'],
effects172[which(effects172$xloc %in% tested & effects172$colClass=='trans' & effects172$source=='div04M'),'xloc']
),show.plot=F)

rgPal <- colorRampPalette(c('black','green'))
col = rgPal(10)[cut(unlist(lapply(attr(trans,"intersections"),length)),breaks=10)]
pdf('bestEffects172transEffectCisParentTransSharing.pdf')
plotVenn4d(round(unlist(lapply(attr(trans,"intersections"),length))/sum(unlist(lapply(attr(trans,"intersections"),length)))*100),Colors=col)
dev.off()

# cis + trans
cisPtrans = venn(list(
effects172[which(effects172$xloc %in% tested & effects172$colClass=='cis+trans' & effects172$source=='div20F'),'xloc'],
effects172[which(effects172$xloc %in% tested & effects172$colClass=='cis+trans' & effects172$source=='div20M'),'xloc'],
effects172[which(effects172$xloc %in% tested & effects172$colClass=='cis+trans' & effects172$source=='div04F'),'xloc'],
effects172[which(effects172$xloc %in% tested & effects172$colClass=='cis+trans' & effects172$source=='div04M'),'xloc']
),show.plot=F)

rgPal <- colorRampPalette(c('black','green'))
# ovelap lengths ==> add zeros to zero overlaps
overlap = round(unlist(lapply(attr(cisPtrans,"intersections"),length))/sum(unlist(lapply(attr(cisPtrans,"intersections"),length)))*100)
overlap$"A:D" = 0
overlap$"B:C" = 0
overlap$"B:D" = 0
overlap$"A:B:C" = 0
overlap$"A:C:D" = 0
overlap = overlap[names(attr(tt,"intersections"))]
col = rgPal(10)[cut(unlist(overlap),breaks=10)]
pdf('bestEffects172cisPtransEffectCisParentTransSharing.pdf')
plotVenn4d(unlist(overlap))
dev.off()

# cis - trans
cisMtrans = venn(list(
effects172[which(effects172$xloc %in% tested & effects172$colClass=='cis-trans' & effects172$source=='div20F'),'xloc'],
effects172[which(effects172$xloc %in% tested & effects172$colClass=='cis-trans' & effects172$source=='div20M'),'xloc'],
effects172[which(effects172$xloc %in% tested & effects172$colClass=='cis-trans' & effects172$source=='div04F'),'xloc'],
effects172[which(effects172$xloc %in% tested & effects172$colClass=='cis-trans' & effects172$source=='div04M'),'xloc']
),show.plot=F)

rgPal <- colorRampPalette(c('black','green'))
# ovelap lengths ==> add zeros to zero overlaps
overlap = round(unlist(lapply(attr(cisMtrans,"intersections"),length))/sum(unlist(lapply(attr(cisMtrans,"intersections"),length)))*100)
overlap$"A:C:D" = 0
overlap$"A:B:C" = 0
overlap = overlap[names(attr(tt,"intersections"))]
col = rgPal(10)[cut(unlist(overlap),breaks=10)]
pdf('bestEffects172cisMtransEffectCisParentTransSharing.pdf')
plotVenn4d(unlist(overlap))
dev.off()


# compensatory
compensatory = venn(list(
effects172[which(effects172$xloc %in% tested & effects172$colClass=='compensatory' & effects172$source=='div20F'),'xloc'],
effects172[which(effects172$xloc %in% tested & effects172$colClass=='compensatory' & effects172$source=='div20M'),'xloc'],
effects172[which(effects172$xloc %in% tested & effects172$colClass=='compensatory' & effects172$source=='div04F'),'xloc'],
effects172[which(effects172$xloc %in% tested & effects172$colClass=='compensatory' & effects172$source=='div04M'),'xloc']
),show.plot=F)

rgPal <- colorRampPalette(c('black','green'))
# ovelap lengths ==> add zeros to zero overlaps
overlap = round(unlist(lapply(attr(compensatory,"intersections"),length))/sum(unlist(lapply(attr(compensatory,"intersections"),length)))*100)
overlap$"A:B:D" = 0
overlap = overlap[names(attr(tt,"intersections"))]
col = rgPal(10)[cut(unlist(overlap),breaks=10)]
pdf('bestEffects172compensatoryEffectCisParentTransSharing.pdf')
plotVenn4d(unlist(overlap),Colors=col)
dev.off()


