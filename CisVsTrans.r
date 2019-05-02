library(GenomicRanges)
library(gtools)
library(ggplot2)
library(qvalue)
library(grid)

setwd('.../')

# clear memory
rm(list=ls())

# DESeq2 normalized counts for parents
normCounts = read.table('.../aseReadCounts_c172_PallF1_infSites_totalCov_DESeqNormalized_STAR_duprem.txt',stringsAsFactors=F,header=T)

# average coverages of parents ran in two lanes
normCounts$c172_P_532_F_TOTCOV = apply(normCounts[,c('c172_P_532_F_lane3_TOTCOV','c172_P_532_F_lane8_TOTCOV')],1,mean)
normCounts$c172_P_533_M_TOTCOV = apply(normCounts[,c('c172_P_533_M_lane3_TOTCOV','c172_P_533_M_lane8_TOTCOV')],1,mean)

#######
# select the individual to be analyzed
ind='20_M'
#######


# ASE table
aseBed = read.table(paste('.../aseTable_c172_F1',ind,'FDR10_minCov10_aseReadCounts_refCorrected_F1normalized_genomicMask_STAR_duprem.txt',sep='_'),header=T,stringsAsFactors=F)
rownames(aseBed) = paste(aseBed[,1],aseBed[,2],sep=':')

# take overlap (only positions tested for ASE)
library(gplots)
m = venn(list(rownames(normCounts),rownames(aseBed)))
mNames = attr(m,"intersections")$"A:B"
normCounts = normCounts[mNames,]
aseBed = aseBed[mNames,]
all(rownames(normCounts) == rownames(aseBed))

# results data frame
P = data.frame(chr=aseBed[,1],pos=aseBed[,2],'FM'=normCounts[,'c172_P_532_F_TOTCOV'],'M'=normCounts[,'c172_P_533_M_TOTCOV'], 'P.fm.m'=foldchange(normCounts[,'c172_P_532_F_TOTCOV'],normCounts[,'c172_P_533_M_TOTCOV']),'F1.fm.m'=aseBed[,'fm.m'])
rownames(P) = rownames(aseBed)

# test DE between parents with binomial test
# use the same depth requirements as for F1 test!!!
binomialTest = function(x,prob){
    if (x[1]>10 | x[2]>10 & x[1] < c(x[1]+x[2])){
        return(unlist(binom.test(c(as.integer(x[1]),as.integer(x[2])),p=prob)['p.value']))
    }
    else {
        return(NA)
    }
}

Pbinom = apply(X=data.frame(normCounts[,'c172_P_532_F_TOTCOV'],normCounts[,'c172_P_533_M_TOTCOV']),MARGIN=1,FUN=binomialTest,prob=sum(normCounts[,'c172_P_532_F_TOTCOV'])/(sum(normCounts[,'c172_P_532_F_TOTCOV'])+sum(normCounts[,'c172_P_533_M_TOTCOV'])))
names(Pbinom) = rownames(normCounts)
Pbinom = na.omit(Pbinom)
Pq = qvalue(Pbinom)$qvalue
P[names(Pq),'parent_qvalue'] = Pq
P = na.omit(P)

# trans test using Fisher's exact test of F1 vs parental allele rations (as in McManus 2010)
# test is significant if allele ratios differ between parental comparison and within F1 hybrids --> P1:P2 != A1:A2
# to be called trans regulated, other requirements apply (no ASE effect in F1's, parents expressed at different levels)
# counts for parents and F1's corected for total number of parent/allele specific reads (this takes into account differences in library size between the parents and mapping bias towards the reference in F1's)
Ptrans=c()
for (i in rownames(P)){ # trans testing for only the genes that have bee tested for ASE!!!
    if (sign(aseBed[i,'ref.alt']) == sign(aseBed[i,'fm.m'])){
        Ptrans[i] = fisher.test(data.frame(c(normCounts[i,'c172_P_532_F_TOTCOV'],normCounts[i,'c172_P_533_M_TOTCOV']),c(aseBed[i,'refCount'],aseBed[i,'altCount'])))$p.value
    }
    else {
        if (sign(aseBed[i,'ref.alt']) != sign(aseBed[i,'fm.m'])){
            Ptrans[i] = fisher.test(data.frame(c(normCounts[i,'c172_P_532_F_TOTCOV'],normCounts[i,'c172_P_533_M_TOTCOV']),c(aseBed[i,'altCount'],aseBed[i,'refCount'])))$p.value
        }
    }
}

Ptrans[Ptrans>1] = 1
Ptrans[Ptrans<0] = 0
Qtrans = qvalue(Ptrans)$qvalue
transSig = Qtrans < 0.1

# ASE test Q values
Qcis = aseBed[rownames(P),'qvalue']
cisSig = Qcis < 0.1

# add Qvalues to P table
P[,'trans_qvalue'] = Qtrans
P[,'cis_qvalue'] = Qcis

# parse all genes that have not been tested for all comparisons
P = na.omit(P)

# labels
P[,'trans'] = FALSE
P[P[,'trans_qvalue']<0.1, 'trans'] = TRUE

P[,'cis'] = FALSE
P[P[,'cis_qvalue']<0.1,'cis'] = TRUE

P[,'parent'] = FALSE
P[P[,'parent_qvalue']<0.1,'parent'] = TRUE

P[,'cis+trans'] = FALSE
P[P[,'cis'] & P[,'trans'] & P[,'parent'] & sign(P[,'P.fm.m'])==sign(P[,'F1.fm.m']),'cis+trans'] = TRUE

P[,'cis-trans'] = FALSE
P[P[,'cis'] & P[,'trans'] & P[,'parent'] & sign(P[,'P.fm.m'])!=sign(P[,'F1.fm.m']),'cis-trans'] = TRUE

P[,'compensatory'] = FALSE
P[P[,'cis'] & P[,'trans'] & P[,'parent']==F,'compensatory'] = TRUE

#################################
# add column for effect classes #
P[,'colClass'] = NA             #
#################################

#---------------------------------------------------------------#
# only cis --> A1!=A2 & P1!=P2 & A1:A2==P1:P2                   #
P[P[,'cis'] & P[,'parent'] & P[,'trans']==F,'colClass'] = 'cis' #
#---------------------------------------------------------------#

#-----------------------------------------------------------------#
# only trans --> A1==A2 & P1!=P2 & A1:A2!=P1:P2                   #
P[P[,'trans'] & P[,'parent'] & P[,'cis']==F,'colClass'] = 'trans' #
#-----------------------------------------------------------------#

#------------------------------------------------------------------------------#
# cis plus trans --> A1!=A2 & P1!=P2 & A1:A2!=P1:P2 & sign(A1:A2)==sign(P1:P2) #
P[P[,'cis+trans'],'colClass'] = 'cis+trans'                                    #
#------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------#
# cis minus trans --> A1!=A2 & P1!=P2 & A1:A2!=P1:P2 & sign(A1:A2)!=sign(P1:P2) #
P[P[,'cis-trans'],'colClass'] = 'cis-trans'                                     #
#-------------------------------------------------------------------------------#

#-------------------------------------------------#
# compensatory --> A1!=A2 & A1:A2!=P1:P2 & P1==P2 #
P[P[,'compensatory'],'colClass'] = 'compensatory' #
#-------------------------------------------------#

#-----------------------------------------------------------#
# optional to plot non-significant or ambiguous genes       #
#P[is.na(P[,'colClass']),'colClass'] = 'conserved/ambiguous' #
#-----------------------------------------------------------#

#----------------------------------------------------------------------------
# non-significant and ambiguous genes                                       #
P[P[,'cis']==F & P[,'parent']==F & P[,'trans']==F,'colClass'] = 'conserved' #
P[is.na(P[,'colClass']),'colClass'] = 'ambiguous'                           #
#----------------------------------------------------------------------------

# save data frame for future analyses
write.table(P,paste('.../effectTable_c172_F1',ind,'FDR10_minCov10_aseReadCounts_refCorrected_F1normalized_genomicMask_STAR_duprem.txt',sep='_'))

