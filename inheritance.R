library(hexbin)
library(gplots)
library(ggplot2)
library(gtools)
library(gdata)


##################################################
# c172                                           #
# with SNP counts of fully informative positions #
# normalized counts from DESeq2                  #
##################################################

# ASEreadCounts
normSumCounts = read.table('aseReadCounts_c172_PallF1_infSites_totalCov_DESeqNormalized_STAR_duprem.txt',stringsAsFactors=F,header=T)
normSumCounts = normSumCounts[,c(1:4,13,14)]

colnames(normSumCounts) = c('F20','M20','F04','M04','M','F')
# take only positions w/ more than 10 reads in any parent (same cutoff for ASE test regarding alleles)
normSumCounts = normSumCounts[apply(normSumCounts>=10,1,any),]

ind = 'F20'
div = 'div20F'

mDiff = log2(normSumCounts[,ind]) - log2(normSumCounts[,'M'])
fDiff = log2(normSumCounts[,ind]) - log2(normSumCounts[,'F'])

names(mDiff) = rownames(normSumCounts)
names(fDiff) = rownames(normSumCounts)

inh = data.frame(fDiff=fDiff,mDiff=mDiff)
rownames(inh) = rownames(normSumCounts)
inh[,'fcF'] = abs(foldchange(normSumCounts[,ind],normSumCounts[,'F']))
inh[,'fcM'] = abs(foldchange(normSumCounts[,ind],normSumCounts[,'M']))
inh = inh[is.infinite(inh[,1])==F & is.infinite(inh[,2])==F,]
inh = na.omit(inh)

##################################################
# c212 & c214 & c363                             #
# with SNP counts of fully informative positions #
# normalized counts from DESeq2                  #
##################################################

normSumCounts = read.table('aseReadCounts_c363_PF1_infSites_totalCov_DESeqNormalized_STAR_duprem.txt',stringsAsFactors=F,header=T)

colnames(normSumCounts) = c('F','M','F01','M01','F02','M02')

# take only positions w/ more than 10 reads in any parent (same cutoff for ASE test regarding alleles)
normSumCounts = normSumCounts[apply(normSumCounts>=10,1,any),]

ind = 'F02'
div = 'div02F'

mDiff = log2(normSumCounts[,ind]) - log2(normSumCounts[,'M'])
fDiff = log2(normSumCounts[,ind]) - log2(normSumCounts[,'F'])

names(mDiff) = rownames(normSumCounts)
names(fDiff) = rownames(normSumCounts)

inh = data.frame(fDiff=fDiff,mDiff=mDiff)
rownames(inh) = rownames(normSumCounts)
inh[,'fcF'] = abs(foldchange(normSumCounts[,ind],normSumCounts[,'F']))
inh[,'fcM'] = abs(foldchange(normSumCounts[,ind],normSumCounts[,'M']))
inh = inh[is.infinite(inh[,1])==F & is.infinite(inh[,2])==F,]
inh = na.omit(inh)

##############
# categories #
##############
inh[,'iClass'] = 'unassigned'

# conserved (F1-P expression diff less than log2(0.25) or folchange 1.25
#inh[inh[,1]>c(-1) & inh[,1]<c(1) & inh[,2]>c(-1) & inh[,2]<c(1),'iClass'] = 'conserved'
inh[inh[,'fcM']<1.25 & inh[,'fcF']<1.25,'iClass'] = 'conserved'

# F parent dominant (hyb similar to F)
#inh[inh[,2]>c(-1) & inh[,2]<c(1),'iClass'] = 'F_dominant'
inh[inh[,'fcF']<1.25 & inh[,'fcM']>1.25,'iClass'] = 'F_dominant'

# M parent dominant (hyb similar to M)
#inh[inh[,1]>c(-1) & inh[,1]<c(1),'iClass'] = 'M_dominant'
inh[inh[,'fcM']<1.25 & inh[,'fcF']>1.25,'iClass'] = 'M_dominant'

# additive (in both directions)
#inh[inh[,'iClass']!='F_dominant' & inh[,'iClass']!='M_dominant' & ((inh[,1]<c(-1) & inh[,2]>c(1)) | inh[,1]>c(1) & inh[,2]<c(-1)),'iClass'] = 'additive'
inh[inh[,'iClass']!='conserved' & inh[,'iClass']!='F_dominant' & inh[,'iClass']!='M_dominant' & ((sign(inh[,'fDiff'])==c(-1) & sign(inh[,'mDiff'])==1) | (sign(inh[,'fDiff'])==c(1) & sign(inh[,'mDiff'])==c(-1))),'iClass'] = 'additive'

# transgressive (underdominant)
#inh[inh[,'iClass']!='F_dominant' & inh[,'iClass']!='M_dominant' & inh[,1]<c(-1) & inh[,2]<c(-1),'iClass'] = 'underdominant'
inh[inh[,'iClass']!='conserved' & inh[,'iClass']!='F_dominant' & inh[,'iClass']!='M_dominant' & sign(inh[,'fDiff'])==c(-1) & sign(inh[,'mDiff'])==c(-1),'iClass'] = 'underdominant'

# transgressive (overdominant)
#inh[inh[,'iClass']!='F_dominant' & inh[,'iClass']!='M_dominant' & inh[,1]>c(1) & inh[,2]>c(1),'iClass'] = 'overdominant'
inh[inh[,'iClass']!='conserved' & inh[,'iClass']!='F_dominant' & inh[,'iClass']!='M_dominant' & sign(inh[,'fDiff'])==c(1) & sign(inh[,'mDiff'])==c(1),'iClass'] = 'overdominant'

# add effect class
effectTable = read.table('AseReadCountsEffects172_bestQvalue_cisParentTrans_STAR_duprem_nosexFDR10.txt',header=T,stringsAsFactors=F)
effectTable = effectTable[which(effectTable$source == div),]

# for SNP-counts
# add xloc and select for unique xloc-effect combinations
xloc = c(effectTable[,'xloc'])
names(xloc) = paste(effectTable$seqnames,effectTable$start,sep=':')
inh[,'xloc'] = xloc[rownames(inh)]
effects = c(effectTable[,'colClass'])
names(effects) = paste(effectTable$seqnames,effectTable$start,sep=':')
inh[,'effectClass'] = effects[rownames(inh)]

############
# DA ratio #
############

# SNP counts - all positions
A = abs(normSumCounts[,'F']-normSumCounts[,'M'])/2
midpoint = c(normSumCounts[,'F']+normSumCounts[,'M'])/2
D = normSumCounts[,ind]-midpoint
DAratio = D/A
names(DAratio) = rownames(normSumCounts)
DAratio = DAratio[DAratio<100 & DAratio>c(-100) & is.infinite(DAratio)==F]

# add DAratio to inh data frame
inh[,'DA'] = DAratio[rownames(inh)]

# STORE inh dataframes
#inh04M=inh
#inh04F=inh
#inh20M=inh
#inh20F=inh

#inh01F=inh
#inh01M=inh
#inh02F=inh
#inh02M=inh

# combine

inhAll = combine(inh20F,inh20M,inh04F,inh04M)
#inhAll = gdata::combine(inh01F,inh01M,inh02F,inh02M)

write.table(inhAll,'AseReadCountsEffects172_bestQvalue_cisParentTrans_STAR_duprem_nosexFDR10_inheritance.txt')

#########
# PLOTS #
#########

# effect class ~ inheritance class
eInh = as.data.frame(table(inhAll[,c('iClass','effectClass')]))
colSum=apply(table(inhAll[,c('iClass','effectClass')]),2,sum)
eInh[,'percentage'] = eInh[,'Freq']/rep(colSum,each=6)*100
# plot
g = ggplot(eInh)
g = g + geom_bar(aes(y=percentage,x=iClass,fill=as.factor(iClass)),stat="identity")
g = g + facet_wrap(~effectClass)
g = g + theme(strip.text = element_text(size=14,face='bold'))
g = g + scale_fill_brewer(palette="Set1")
g

# inheritance class ~ effect class
#eInh = as.data.frame(table(inh[,c('effectClass','iClass')]))
eInh = as.data.frame(table(inh[inh[,'effectClass']!='conserved' & inh[,'effectClass']!='ambiguous',c('effectClass','iClass')]))
colSum=apply(table(inh[,c('effectClass','iClass')]),2,sum)
#eInh[,'percentage'] = eInh[,'Freq']/rep(colSum,each=7)*100
eInh[,'percentage'] = eInh[,'Freq']/rep(colSum,each=5)*100
# plot
g = ggplot(eInh)
g = g + geom_bar(aes(y=percentage,x=effectClass,fill=as.factor(effectClass)),stat="identity")
g = g + facet_wrap(~iClass)
g = g + theme(strip.text = element_text(size=14,face='bold'))
g

#####################################
# D/A ratio for continuous analyses #
#####################################

# load data
inhAll = read.table('AseReadCountsEffects172_bestQvalue_cisParentTrans_STAR_duprem_nosexFDR10_inheritance.txt',stringsAsFactors=F)

# verify that this makes sense
inh[which(inh[,'iClass']=='additive'),'DA']
hist(inh[which(inh[,'iClass']=='additive'),'DA'])

inhAll$effectClass = factor(inhAll$effectClass,levels=c('cis','trans','cis+trans','cis-trans','compensatory')) # re-order classes
cols = brewer.pal(6,"Set1")
cols[6] = "#999999"
g = ggplot(na.omit(inhAll[which(inhAll$effectClass %in% c('ambiguous','conserved')==F & is.na(inhAll$xloc) == F),]))
g = g + geom_density(aes(DA,fill=as.factor(effectClass)))
g = g + xlim(-4,4)
g = g + facet_wrap(~effectClass,nrow=1,drop=T)
g = g + coord_flip()
g = g + xlab('Dominance:Additivity')
g = g + ylab('')
g = g + theme_classic()
g = g + theme(legend.position="none",axis.text.x = element_text(size=24,colour = "black"),axis.text.y = element_text(size=24,colour = "black"), text = element_text(size=24,colour = "black"),strip.background = element_blank(), strip.text = element_blank())
g = g + scale_fill_manual(values=cols)
g
# save as DAratioc172BestEffectsCisParentTrans.pdf






