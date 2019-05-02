
########
# c172 #
########
library(DESeq2)

# all F1s w/ RNAseq masked ref
cov = read.table('aseReadCounts_c172_PallF1_infSites_totalCov_STAR_duprem.txt',row.names='pos')

# verify parent lane correlation
cor(na.omit(cov[,c('c172_P_533_M_lane3_TOTCOV','c172_P_533_M_lane8_TOTCOV')]))
cor(na.omit(cov[,c('c172_P_532_F_lane3_TOTCOV','c172_P_532_F_lane8_TOTCOV')]))

# exclude chrXIX
#cov = cov[grepl('chrXIX',rownames(cov))==F,]

cov[is.na(cov)] = 0
cov = cov[,2:ncol(cov)]
cov = as.numeric(cov)
colSums(cov)

# add parent lanes together
#cov$c172_P_533_M_TOTCOV = cov$c172_P_533_M_lane3_TOTCOV + cov$c172_P_533_M_lane8_TOTCOV
#cov$c172_P_532_F_TOTCOV = cov$c172_P_532_F_lane3_TOTCOV + cov$c172_P_532_F_lane8_TOTCOV

# remove parent lane columns
#cov = cov[,which(colnames(cov) %in% c('c172_P_532_F_lane3_TOTCOV','c172_P_532_F_lane8_TOTCOV','c172_P_533_M_lane3_TOTCOV','c172_P_533_M_lane8_TOTCOV') == F)]

design = data.frame(row.names=colnames(cov),
sex=rep(c('F','M'),8),
generation=c(rep('F1',12),rep('P',4)))

icounts = DESeqDataSetFromMatrix(countData=cov, colData=design, design=~1)
icounts = estimateSizeFactors(icounts)

# plot
plot(colSums(counts(icounts)),sizeFactors(icounts),type='n')
text(colSums(counts(icounts)),sizeFactors(icounts),colnames(counts(icounts)),cex=.75)
abline(lm(sizeFactors(icounts) ~ colSums(counts(icounts))))

write.table(counts(icounts,normalized=T),'aseReadCounts_c172_PallF1_infSites_totalCov_DESeqNormalized_STAR_duprem.txt')


# allele-specific counts
alleleCov = read.table('aseReadCounts_c172_PallF1_infSites_STAR_duprem.txt',row.names='pos')

# exclude chrXIX
#alleleCov = alleleCov[grepl('chrXIX',rownames(alleleCov))==F,]

alleleCov[is.na(alleleCov)] = 0
alleleCov = alleleCov[,2:ncol(alleleCov)]
alleleCov = as.numeric(alleleCov)

#########
# for F1s

# remove parent lane columns
alleleCov = alleleCov[,c('c172_F1_20_F_REFCOV','c172_F1_20_F_ALTCOV','c172_F1_20_M_REFCOV','c172_F1_20_M_ALTCOV','c172_F1_04_F_REFCOV','c172_F1_04_F_ALTCOV','c172_F1_04_M_REFCOV','c172_F1_04_M_ALTCOV')]

aDesign = data.frame(row.names=colnames(alleleCov),sex=rep(c('F','F','M','M'),2),generation=c(rep('F1',8)))
acounts = DESeqDataSetFromMatrix(countData=alleleCov, colData=aDesign, design=~1)
acounts = estimateSizeFactors(acounts)

plot(colSums(counts(acounts)),sizeFactors(acounts),type='n')
text(colSums(counts(acounts)),sizeFactors(acounts),colnames(counts(acounts)),cex=.75)
abline(lm(sizeFactors(acounts) ~ colSums(counts(acounts))))

colSums(counts(acounts,normalized=T))

write.table(counts(acounts,normalized=T),'aseReadCounts_c172_F1_infSites_alleleCov_DESeqNormalized_STAR_duprem.txt')



###################
# for salinity F1s

alleleCov = read.table('aseReadCounts_c172_PallF1_infSites_STAR_duprem.txt',row.names='pos')

# remove parent lane columns
alleleCov = alleleCov[,which(colnames(alleleCov) %in% c('c172_P_532_F_lane3_ALTCOV','c172_P_532_F_lane8_ALTCOV','c172_P_533_M_lane3_ALTCOV','c172_P_533_M_lane8_ALTCOV','c172_P_532_F_lane3_REFCOV','c172_P_532_F_lane8_REFCOV','c172_P_533_M_lane3_REFCOV','c172_P_533_M_lane8_REFCOV') == F)]

alleleCov[is.na(alleleCov)] = 0
alleleCov = alleleCov[,2:ncol(alleleCov)]
alleleCov = as.numeric(alleleCov)

# for F1s
aDesign = data.frame(row.names=colnames(alleleCov),sex=rep(c('F','F','M','M'),6),generation=c(rep('F1',24)))
acounts = DESeqDataSetFromMatrix(countData=alleleCov, colData=aDesign, design=~1)
acounts = estimateSizeFactors(acounts)

plot(colSums(counts(acounts)),sizeFactors(acounts),type='n')
text(colSums(counts(acounts)),sizeFactors(acounts),colnames(counts(acounts)),cex=.75)
abline(lm(sizeFactors(acounts) ~ colSums(counts(acounts))))

colSums(counts(acounts,normalized=T))

write.table(counts(acounts,normalized=T),'.../aseReadCounts_c172_allF1_infSites_alleleCov_DESeqNormalized_STAR_duprem.txt')



