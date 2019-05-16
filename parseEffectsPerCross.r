library(ggplot2)
library(reshape)
library(gtools)
library(gdata)


MF = read.table('gene_exp.diff',stringsAsFactors=F,sep='\t',header=T)
MFgenes = MF$gene_id[which(MF$q_value<0.1)]

##
# Merge all F1 effectTables for a given cross
##

##################################################
# pick one effect per xloc per F1 and then combine
##################################################

########
# c172 #
########
e20F = read.table('effectTable_c172_F1_20_F_FDR10_minCov10_aseReadCounts_refCorrected_F1normalized_genomicMask_STAR_duprem.txt',header=T,stringsAsFactors=F)
e20M = read.table('effectTable_c172_F1_20_M_FDR10_minCov10_aseReadCounts_refCorrected_F1normalized_genomicMask_STAR_duprem.txt',header=T,stringsAsFactors=F)
e04F = read.table('effectTable_c172_F1_04_F_FDR10_minCov10_aseReadCounts_refCorrected_F1normalized_genomicMask_STAR_duprem.txt',header=T,stringsAsFactors=F)
e04M = read.table('effectTable_c172_F1_04_M_FDR10_minCov10_aseReadCounts_refCorrected_F1normalized_genomicMask_STAR_duprem.txt',header=T,stringsAsFactors=F)

# exclude genes showing DE between sexes
e20F = e20F[which(e20F$xloc %in% MFgenes == F),]
e20M = e20M[which(e20M$xloc %in% MFgenes == F),]
e04F = e04F[which(e04F$xloc %in% MFgenes == F),]
e04M = e04M[which(e04M$xloc %in% MFgenes == F),]

# OPTIONAL - filter for mono allelic expression
e20F = e20F[is.infinite(e20F$F1.fm.m)==F & is.infinite(e20F$P.fm.m)==F,]
e20M = e20M[is.infinite(e20M$F1.fm.m)==F & is.infinite(e20M$P.fm.m)==F,]
e04F = e04F[is.infinite(e04F$F1.fm.m)==F & is.infinite(e04F$P.fm.m)==F,]
e04M = e04M[is.infinite(e04M$F1.fm.m)==F & is.infinite(e04M$P.fm.m)==F,]

# reorder classes
e20F$colClass = as.character(factor(e20F$colClass, levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous')))
e20M$colClass = as.character(factor(e20M$colClass, levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous')))
e04F$colClass = as.character(factor(e04F$colClass, levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous')))
e04M$colClass = as.character(factor(e04M$colClass, levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous')))

# exclude SNP loci without xloc annotations - intergenic expression that is not supported by Cufflinks transcript models
e20F = e20F[which(is.na(e20F$xloc)==F),]
e20M = e20M[which(is.na(e20M$xloc)==F),]
e04F = e04F[which(is.na(e04F$xloc)==F),]
e04M = e04M[which(is.na(e04M$xloc)==F),]

# pick best effect per xloc - multiply cis qvalue and trans qvalue
div20F=as.data.frame(matrix(ncol=length(colnames(e20F)),nrow=0))
colnames(div20F) = colnames(e20F)
for (i in na.omit(unique(e20F$xloc))){
    # order F1s based on q values - multiply cis and trans q-values
    # modified 3/5/2017 --> multiply cis, trans and parent q-values
    e172 = e20F[which(e20F$xloc %in% i),]
    e172 = e172[order(e172$cis_qvalue*e172$trans_qvalue*e172$parent_qvalue),]
    # pick the top effect
    div20F[i,] = as.data.frame(e172[1,])
}
div20M=as.data.frame(matrix(ncol=length(colnames(e20M)),nrow=0))
colnames(div20M) = colnames(e20M)
for (i in na.omit(unique(e20M$xloc))){
    # order F1s based on q values - multiply cis and trans q-values
    # modified 3/5/2017 --> multiply cis, trans and parent q-values
    e172 = e20M[which(e20M$xloc %in% i),]
    e172 = e172[order(e172$cis_qvalue*e172$trans_qvalue*e172$parent_qvalue),]
    # pick the top effect
    div20M[i,] = as.data.frame(e172[1,])
}
div04F=as.data.frame(matrix(ncol=length(colnames(e04F)),nrow=0))
colnames(div04F) = colnames(e04F)
for (i in na.omit(unique(e04F$xloc))){
    # order F1s based on q values - multiply cis and trans q-values
    # modified 3/5/2017 --> multiply cis, trans and parent q-values
    e172 = e04F[which(e04F$xloc %in% i),]
    e172 = e172[order(e172$cis_qvalue*e172$trans_qvalue*e172$parent_qvalue),]
    # pick the top effect
    div04F[i,] = as.data.frame(e172[1,])
}
div04M=as.data.frame(matrix(ncol=length(colnames(e04M)),nrow=0))
colnames(div04M) = colnames(e04M)
for (i in na.omit(unique(e04M$xloc))){
    # order F1s based on q values - multiply cis and trans q-values
    # modified 3/5/2017 --> multiply cis, trans and parent q-values
    e172 = e04M[which(e04M$xloc %in% i),]
    e172 = e172[order(e172$cis_qvalue*e172$trans_qvalue*e172$parent_qvalue),]
    # pick the top effect
    div04M[i,] = as.data.frame(e172[1,])
}

# combine
combinedEffects172 = gdata::combine(div20F,div20M,div04F,div04M)
effects172 = cbind(table(div20F$colClass),table(div20M$colClass),table(div04F$colClass),table(div04M$colClass))
colnames(effects172) = c('20F','20M','04F','04M')

# save
write.table(combinedEffects172,'AseReadCountsEffects172_bestQvalue_cisParentTrans_STAR_duprem_nosexFDR10.txt')
write.table(effects172,'effectTable172_bestQvalue_cisParentTrans_STAR_duprem_nosexFDR10.txt')

# plot
effects = gdata::combine(melt(effects172[rownames(effects172) %in% c('ambiguous',"conserved")==F,]),names=c('Tyne'))
colnames(effects) = c('class','ind','frequency','cross')
effects$class = factor(effects$class,levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved'))
#postscript('boxplotAseReadCountsEffects172_bestQvalue_STAR_duprem.ps',width=7)
g = ggplot(effects)
g = g + geom_boxplot(aes(as.factor(class),frequency,fill=as.factor(class)))
g = g + geom_jitter(aes(as.factor(class),frequency),height=0,width=0.3)
g = g + facet_wrap(~cross)
g = g + theme(axis.text.x = element_text(angle=45, hjust=1,size=24,colour='black'),axis.title.x=element_blank(),panel.background=element_rect(fill='white',colour='grey78'),panel.grid.major=element_line(colour='grey78'),text=element_text(size=22),legend.position="none")
g = g + ylab('Frequency')
g = g + scale_fill_brewer(palette="Set1")
g
#dev.off()


########
# c212 #
########
rm(list=ls())

MF = read.table('gene_exp.diff',stringsAsFactors=F,sep='\t',header=T)
MFgenes = MF$gene_id[which(MF$q_value<0.1)]

e01F = read.table('effectTable_c212_F1_1_FM_FDR10_minCov10_aseReadCounts_refCorrected_F1normalized_genomicMask_STAR_duprem.txt',header=T,stringsAsFactors=F)
e01M = read.table('effectTable_c212_F1_1_M_FDR10_minCov10_aseReadCounts_refCorrected_F1normalized_genomicMask_STAR_duprem.txt',header=T,stringsAsFactors=F)
e02F = read.table('effectTable_c212_F1_2_FM_FDR10_minCov10_aseReadCounts_refCorrected_F1normalized_genomicMask_STAR_duprem.txt',header=T,stringsAsFactors=F)
e02M = read.table('effectTable_c212_F1_2_M_FDR10_minCov10_aseReadCounts_refCorrected_F1normalized_genomicMask_STAR_duprem.txt',header=T,stringsAsFactors=F)

# exclude genes showing DE between sexes
e01F = e01F[which(e01F$xloc %in% MFgenes == F),]
e01M = e01M[which(e01M$xloc %in% MFgenes == F),]
e02F = e02F[which(e02F$xloc %in% MFgenes == F),]
e02M = e02M[which(e02M$xloc %in% MFgenes == F),]

# OPTIONAL - filter for mono allelic expression
e01F = e01F[is.infinite(e01F$F1.fm.m)==F & is.infinite(e01F$P.fm.m)==F,]
e01M = e01M[is.infinite(e01M$F1.fm.m)==F & is.infinite(e01M$P.fm.m)==F,]
e02F = e02F[is.infinite(e02F$F1.fm.m)==F & is.infinite(e02F$P.fm.m)==F,]
e02M = e02M[is.infinite(e02M$F1.fm.m)==F & is.infinite(e02M$P.fm.m)==F,]

# reorder classes
e01F$colClass = as.character(factor(e01F$colClass, levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous')))
e01M$colClass = as.character(factor(e01M$colClass, levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous')))
e02F$colClass = as.character(factor(e02F$colClass, levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous')))
e02M$colClass = as.character(factor(e02M$colClass, levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous')))

# exclude SNP loci without xloc annotations - intergenic expression that is not supported by Cufflinks transcript models
e01F = e01F[which(is.na(e01F$xloc)==F),]
e01M = e01M[which(is.na(e01M$xloc)==F),]
e02F = e02F[which(is.na(e02F$xloc)==F),]
e02M = e02M[which(is.na(e02M$xloc)==F),]

# pick best effect per xloc - multiply cis qvalue and trans qvalue
div01F=as.data.frame(matrix(ncol=length(colnames(e01F)),nrow=0))
colnames(div01F) = colnames(e01F)
for (i in na.omit(unique(e01F$xloc))){
  # order F1s based on q values - multiply cis and trans q-values
  e212 = e01F[which(e01F$xloc %in% i),]
  e212 = e212[order(e212$cis_qvalue*e212$trans_qvalue*e212$parent_qvalue),]
  # pick the top effect
  div01F[i,] = as.data.frame(e212[1,])
}
div01M=as.data.frame(matrix(ncol=length(colnames(e01M)),nrow=0))
colnames(div01M) = colnames(e01M)
for (i in na.omit(unique(e01M$xloc))){
  # order F1s based on q values - multiply cis and trans q-values
  # modified 3/5/2017 --> multiply cis, trans and parent q-values
  e212 = e01M[which(e01M$xloc %in% i),]
  e212 = e212[order(e212$cis_qvalue*e212$trans_qvalue*e212$parent_qvalue),]
  # pick the top effect
  div01M[i,] = as.data.frame(e212[1,])
}
div02F=as.data.frame(matrix(ncol=length(colnames(e02F)),nrow=0))
colnames(div02F) = colnames(e02F)
for (i in na.omit(unique(e02F$xloc))){
  # order F1s based on q values - multiply cis and trans q-values
  # modified 3/5/2017 --> multiply cis, trans and parent q-values
  e212 = e02F[which(e02F$xloc %in% i),]
  e212 = e212[order(e212$cis_qvalue*e212$trans_qvalue*e212$parent_qvalue),]
  # pick the top effect
  div02F[i,] = as.data.frame(e212[1,])
}
div02M=as.data.frame(matrix(ncol=length(colnames(e02M)),nrow=0))
colnames(div02M) = colnames(e02M)
for (i in na.omit(unique(e02M$xloc))){
  # order F1s based on q values - multiply cis and trans q-values
  # modified 3/5/2017 --> multiply cis, trans and parent q-values
  e212 = e02M[which(e02M$xloc %in% i),]
  e212 = e212[order(e212$cis_qvalue*e212$trans_qvalue*e212$parent_qvalue),]
  # pick the top effect
  div02M[i,] = as.data.frame(e212[1,])
}

# combine
combinedEffects212 = gdata::combine(div01F,div01M,div02F,div02M)
effects212 = cbind(table(div01F$colClass),table(div01M$colClass),table(div02F$colClass),table(div02M$colClass))
colnames(effects212) = c('01F','01M','02F','02M')

# save
write.table(combinedEffects212,'AseReadCountsEffects212_bestQvalue_cisParentTrans_STAR_duprem_nosexFDR10.txt')
write.table(effects212,'effectTable212_bestQvalue_cisParentTrans_STAR_duprem_nosexFDR10.txt')


# plot
effects = gdata::combine(melt(effects212[rownames(effects212) %in% c('ambiguous')==F,]),names=c('Forss'))
colnames(effects) = c('class','ind','frequency','cross')
effects$class = factor(effects$class,levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved'))
#postscript('/fml/chones/projects/PJ035_JPV_CisTrans/ASE/boxplotAseReadCountsEffects212_bestQvalue_STAR_duprem.ps',width=7)
g = ggplot(effects)
g = g + geom_boxplot(aes(as.factor(class),frequency,fill=as.factor(class)))
g = g + geom_jitter(aes(as.factor(class),frequency),height=0,width=0.3)
g = g + facet_wrap(~cross)
g = g + theme(axis.text.x = element_text(angle=45, hjust=1,size=24,colour='black'),axis.title.x=element_blank(),panel.background=element_rect(fill='white',colour='grey78'),panel.grid.major=element_line(colour='grey78'),text=element_text(size=22),legend.position="none")
g = g + ylab('Frequency')
g = g + scale_fill_brewer(palette="Set1")
g
#dev.off()


########
# c214 #
########
rm(list=ls())

MF = read.table('gene_exp.diff',stringsAsFactors=F,sep='\t',header=T)
MFgenes = MF$gene_id[which(MF$q_value<0.1)]

e01F = read.table('effectTable_c214_F1_1_FM_FDR10_minCov10_aseReadCounts_refCorrected_F1normalized_genomicMask_STAR_duprem.txt',header=T,stringsAsFactors=F)
e01M = read.table('effectTable_c214_F1_1_M_FDR10_minCov10_aseReadCounts_refCorrected_F1normalized_genomicMask_STAR_duprem.txt',header=T,stringsAsFactors=F)
e02F = read.table('effectTable_c214_F1_2_FM_FDR10_minCov10_aseReadCounts_refCorrected_F1normalized_genomicMask_STAR_duprem.txt',header=T,stringsAsFactors=F)
e02M = read.table('effectTable_c214_F1_2_M_FDR10_minCov10_aseReadCounts_refCorrected_F1normalized_genomicMask_STAR_duprem.txt',header=T,stringsAsFactors=F)

# exclude genes showing DE between sexes
e01F = e01F[which(e01F$xloc %in% MFgenes == F),]
e01M = e01M[which(e01M$xloc %in% MFgenes == F),]
e02F = e02F[which(e02F$xloc %in% MFgenes == F),]
e02M = e02M[which(e02M$xloc %in% MFgenes == F),]

# OPTIONAL - filter for mono allelic expression
e01F = e01F[is.infinite(e01F$F1.fm.m)==F & is.infinite(e01F$P.fm.m)==F,]
e01M = e01M[is.infinite(e01M$F1.fm.m)==F & is.infinite(e01M$P.fm.m)==F,]
e02F = e02F[is.infinite(e02F$F1.fm.m)==F & is.infinite(e02F$P.fm.m)==F,]
e02M = e02M[is.infinite(e02M$F1.fm.m)==F & is.infinite(e02M$P.fm.m)==F,]

# reorder classes
e01F$colClass = as.character(factor(e01F$colClass, levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous')))
e01M$colClass = as.character(factor(e01M$colClass, levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous')))
e02F$colClass = as.character(factor(e02F$colClass, levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous')))
e02M$colClass = as.character(factor(e02M$colClass, levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous')))

# pick best effect per xloc - multiply cis qvalue and trans qvalue
div01F=as.data.frame(matrix(ncol=length(colnames(e01F)),nrow=0))
colnames(div01F) = colnames(e01F)
for (i in na.omit(unique(e01F$xloc))){
  # order F1s based on q values - multiply cis and trans q-values
  # modified 3/5/2017 --> multiply cis, trans and parent q-values
  e214 = e01F[which(e01F$xloc %in% i),]
  e214 = e214[order(e214$cis_qvalue*e214$trans_qvalue*e214$parent_qvalue),]
  # pick the top effect
  div01F[i,] = as.data.frame(e214[1,])
}
div01M=as.data.frame(matrix(ncol=length(colnames(e01M)),nrow=0))
colnames(div01M) = colnames(e01M)
for (i in na.omit(unique(e01M$xloc))){
  # order F1s based on q values - multiply cis and trans q-values
  # modified 3/5/2017 --> multiply cis, trans and parent q-values
  e214 = e01M[which(e01M$xloc %in% i),]
  e214 = e214[order(e214$cis_qvalue*e214$trans_qvalue*e214$parent_qvalue),]
  # pick the top effect
  div01M[i,] = as.data.frame(e214[1,])
}
div02F=as.data.frame(matrix(ncol=length(colnames(e02F)),nrow=0))
colnames(div02F) = colnames(e02F)
for (i in na.omit(unique(e02F$xloc))){
  # order F1s based on q values - multiply cis and trans q-values
  # modified 3/5/2017 --> multiply cis, trans and parent q-values
  e214 = e02F[which(e02F$xloc %in% i),]
  e214 = e214[order(e214$cis_qvalue*e214$trans_qvalue*e214$parent_qvalue),]
  # pick the top effect
  div02F[i,] = as.data.frame(e214[1,])
}
div02M=as.data.frame(matrix(ncol=length(colnames(e02M)),nrow=0))
colnames(div02M) = colnames(e02M)
for (i in na.omit(unique(e02M$xloc))){
  # order F1s based on q values - multiply cis and trans q-values
  # modified 3/5/2017 --> multiply cis, trans and parent q-values
  e214 = e02M[which(e02M$xloc %in% i),]
  e214 = e214[order(e214$cis_qvalue*e214$trans_qvalue*e214$parent_qvalue),]
  # pick the top effect
  div02M[i,] = as.data.frame(e214[1,])
}


# combine
combinedEffects214 = gdata::combine(div01F,div01M,div02F,div02M)
effects214 = cbind(table(div01F$colClass),table(div01M$colClass),table(div02F$colClass),table(div02M$colClass))
colnames(effects214) = c('01F','01M','02F','02M')

# save
write.table(combinedEffects214,'AseReadCountsEffects214_bestQvalue_cisParentTrans_STAR_duprem_nosexFDR10.txt')
write.table(effects214,'effectTable214_bestQvalue_cisParentTrans_STAR_duprem_nosexFDR10.txt')


# plot
effects = gdata::combine(melt(effects214[rownames(effects214) %in% c('ambiguous')==F,]),names=c('Shiel'))
colnames(effects) = c('class','ind','frequency','cross')
effects$class = factor(effects$class,levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved'))
#postscript('/fml/chones/projects/PJ035_JPV_CisTrans/ASE/boxplotAseReadCountsEffects214_bestQvalue_STAR_duprem.ps',width=7)
g = ggplot(effects)
g = g + geom_boxplot(aes(as.factor(class),frequency,fill=as.factor(class)))
g = g + geom_jitter(aes(as.factor(class),frequency),height=0,width=0.3)
g = g + facet_wrap(~cross)
g = g + theme(axis.text.x = element_text(angle=45, hjust=1,size=24,colour='black'),axis.title.x=element_blank(),panel.background=element_rect(fill='white',colour='grey78'),panel.grid.major=element_line(colour='grey78'),text=element_text(size=22),legend.position="none")
g = g + ylab('Frequency')
g = g + scale_fill_brewer(palette="Set1")
g
#dev.off()


########
# c363 #
########
rm(list=ls())

MF = read.table('gene_exp.diff',stringsAsFactors=F,sep='\t',header=T)
MFgenes = MF$gene_id[which(MF$q_value<0.1)]

e01F = read.table('effectTable_c363_F1_1_F_FDR10_minCov10_aseReadCounts_refCorrected_F1normalized_genomicMask_STAR_duprem.txt',header=T,stringsAsFactors=F)
e01M = read.table('effectTable_c363_F1_1_M_FDR10_minCov10_aseReadCounts_refCorrected_F1normalized_genomicMask_STAR_duprem.txt',header=T,stringsAsFactors=F)
e02F = read.table('effectTable_c363_F1_2_F_FDR10_minCov10_aseReadCounts_refCorrected_F1normalized_genomicMask_STAR_duprem.txt',header=T,stringsAsFactors=F)
e02M = read.table('effectTable_c363_F1_2_M_FDR10_minCov10_aseReadCounts_refCorrected_F1normalized_genomicMask_STAR_duprem.txt',header=T,stringsAsFactors=F)

# exclude genes showing DE between sexes
e01F = e01F[which(e01F$xloc %in% MFgenes == F),]
e01M = e01M[which(e01M$xloc %in% MFgenes == F),]
e02F = e02F[which(e02F$xloc %in% MFgenes == F),]
e02M = e02M[which(e02M$xloc %in% MFgenes == F),]

# OPTIONAL - filter for mono allelic expression
e01F = e01F[is.infinite(e01F$F1.fm.m)==F & is.infinite(e01F$P.fm.m)==F,]
e01M = e01M[is.infinite(e01M$F1.fm.m)==F & is.infinite(e01M$P.fm.m)==F,]
e02F = e02F[is.infinite(e02F$F1.fm.m)==F & is.infinite(e02F$P.fm.m)==F,]
e02M = e02M[is.infinite(e02M$F1.fm.m)==F & is.infinite(e02M$P.fm.m)==F,]

# reorder classes
e01F$colClass = as.character(factor(e01F$colClass, levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous')))
e01M$colClass = as.character(factor(e01M$colClass, levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous')))
e02F$colClass = as.character(factor(e02F$colClass, levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous')))
e02M$colClass = as.character(factor(e02M$colClass, levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous')))

# pick best effect per xloc - multiply cis qvalue and trans qvalue
div01F=as.data.frame(matrix(ncol=length(colnames(e01F)),nrow=0))
colnames(div01F) = colnames(e01F)
for (i in na.omit(unique(e01F$xloc))){
  # order F1s based on q values - multiply cis and trans q-values
  # modified 3/5/2017 --> multiply cis, trans and parent q-values
  e363 = e01F[which(e01F$xloc %in% i),]
  e363 = e363[order(e363$cis_qvalue*e363$trans_qvalue*e363$parent_qvalue),]
  # pick the top effect
  div01F[i,] = as.data.frame(e363[1,])
}
div01M=as.data.frame(matrix(ncol=length(colnames(e01M)),nrow=0))
colnames(div01M) = colnames(e01M)
for (i in na.omit(unique(e01M$xloc))){
  # order F1s based on q values - multiply cis and trans q-values
  # modified 3/5/2017 --> multiply cis, trans and parent q-values
  e363 = e01M[which(e01M$xloc %in% i),]
  e363 = e363[order(e363$cis_qvalue*e363$trans_qvalue*e363$parent_qvalue),]
  # pick the top effect
  div01M[i,] = as.data.frame(e363[1,])
}
div02F=as.data.frame(matrix(ncol=length(colnames(e02F)),nrow=0))
colnames(div02F) = colnames(e02F)
for (i in na.omit(unique(e02F$xloc))){
  # order F1s based on q values - multiply cis and trans q-values
  # modified 3/5/2017 --> multiply cis, trans and parent q-values
  e363 = e02F[which(e02F$xloc %in% i),]
  e363 = e363[order(e363$cis_qvalue*e363$trans_qvalue*e363$parent_qvalue),]
  # pick the top effect
  div02F[i,] = as.data.frame(e363[1,])
}
div02M=as.data.frame(matrix(ncol=length(colnames(e02M)),nrow=0))
colnames(div02M) = colnames(e02M)
for (i in na.omit(unique(e02M$xloc))){
  # order F1s based on q values - multiply cis and trans q-values
  # modified 3/5/2017 --> multiply cis, trans and parent q-values
  e363 = e02M[which(e02M$xloc %in% i),]
  e363 = e363[order(e363$cis_qvalue*e363$trans_qvalue*e363$parent_qvalue),]
  # pick the top effect
  div02M[i,] = as.data.frame(e363[1,])
}

# combine
combinedEffects363 = gdata::combine(div01F,div01M,div02F,div02M)
effects363 = cbind(table(div01F$colClass),table(div01M$colClass),table(div02F$colClass),table(div02M$colClass))
colnames(effects363) = c('01F','01M','02F','02M')

# save
write.table(combinedEffects363,'AseReadCountsEffects363_bestQvalue_cisParentTrans_STAR_duprem_nosexFDR10.txt')
write.table(effects363,'effectTable363_bestQvalue_cisParentTrans_STAR_duprem_nosexFDR10.txt')


# plot
effects = gdata::combine(melt(effects363[rownames(effects363) %in% c('ambiguous')==F,]),names=c('Litc'))
colnames(effects) = c('class','ind','frequency','cross')
effects$class = factor(effects$class,levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved'))
#postscript('/fml/chones/projects/PJ035_JPV_CisTrans/ASE/boxplotAseReadCountsEffects363_bestQvalue_STAR_duprem.ps',width=7)
g = ggplot(effects)
g = g + geom_boxplot(aes(as.factor(class),frequency,fill=as.factor(class)))
g = g + geom_jitter(aes(as.factor(class),frequency),height=0,width=0.3)
g = g + facet_wrap(~cross)
g = g + theme(axis.text.x = element_text(angle=45, hjust=1,size=24,colour='black'),axis.title.x=element_blank(),panel.background=element_rect(fill='white',colour='grey78'),panel.grid.major=element_line(colour='grey78'),text=element_text(size=22),legend.position="none")
g = g + ylab('Frequency')
g = g + scale_fill_brewer(palette="Set1")
g
#dev.off()



#########
# PLOTS #
#########
library(RColorBrewer)
library(reshape)
library(ggplot2)

effects363 = as.matrix(read.table('effectTable363_bestQvalue_cisParentTrans_STAR_duprem_nosexFDR10.txt',header=T,stringsAsFactors=F))
effects212 = as.matrix(read.table('effectTable212_bestQvalue_cisParentTrans_STAR_duprem_nosexFDR10.txt',header=T,stringsAsFactors=F))
effects214 = as.matrix(read.table('effectTable214_bestQvalue_cisParentTrans_STAR_duprem_nosexFDR10.txt',header=T,stringsAsFactors=F))
effects172 = as.matrix(read.table('effectTable172_bestQvalue_cisParentTrans_STAR_duprem_nosexFDR10.txt',header=T,stringsAsFactors=F))

# as percentage
#x11(height=5,width=7)
combinedEffectsPCT = gdata::combine(melt(effects172/apply(effects172,2,sum)),melt(effects212/apply(effects212,2,sum)),melt(effects214/apply(effects214,2,sum)),melt(effects363/apply(effects363,2,sum)),names=c('Tyne','Forss','Shiel','Litc'))
colnames(combinedEffectsPCT) = c('class','ind','frequency','cross')
combinedEffectsPCT=combinedEffectsPCT[which(combinedEffectsPCT$class!='ambiguous'),]
combinedEffectsPCT$class = factor(combinedEffectsPCT$class,levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved'))
cols = brewer.pal(6,"Set1")
cols[6] = "#999999"
g = ggplot(combinedEffectsPCT)
#g = g + geom_boxplot(aes(class,frequency,fill=as.factor(class)))
g = g + geom_jitter(aes(class,frequency,colour=as.factor(class)),height=0,width=0.3,size=2,alpha=0.7)
g = g + facet_wrap(~cross,nrow=1,ncol=4)
g = g + theme(axis.text.x = element_text(angle=45, hjust=1,size=24,colour='black'),axis.title.x=element_blank(),panel.background=element_rect(fill='white',colour='grey78'),panel.grid.major=element_line(colour='grey78'),text=element_text(size=22),legend.position="none")
g = g + theme_classic()
g = g + theme(legend.position = "none")
g = g + ylab('Percentage')
g = g + scale_colour_manual(values=cols)
g
# save as boxplotAseReadCountsEffectsAllCrossesPercentage.pdf


##############
# Scatterpots
##############
rm(list=ls())

# load
effectTable = read.table('AseReadCountsEffects363_bestQvalue_cisParentTrans_STAR_duprem_nosexFDR10.txt',header=T,stringsAsFactors=F)

# fold change plot
effectTable$colClass = factor(effectTable$colClass,levels=c('cis','trans','cis+trans','cis-trans','compensatory','ambiguous','conserved'))

marker = list(color = brewer.pal(5, "Set1"))

Pplot = effectTable[effectTable[,'colClass']!='ambiguous' & effectTable[,'colClass']!='conserved',]
#Pplot = effectTable[effectTable[,'colClass']=='cis',]
#Pplot = Pplot[which(sign(Pplot$P.fm.m) == sign(Pplot$F1.fm.m)),]
#Pplot = effectTable[effectTable[,'colClass']=='trans',]
#Pplot = effectTable[effectTable[,'colClass']=='cis+trans',]
#Pplot = effectTable[effectTable[,'colClass']=='cis-trans',]
#Pplot = effectTable[effectTable[,'colClass']=='compensatory',]

g = ggplot(Pplot)
g = g + geom_point(aes(foldchange2logratio(P.fm.m),foldchange2logratio(F1.fm.m),colour=colClass),cex=2,alpha=0.7)
g = g + theme(panel.background=element_rect(fill='white',colour='black'),text=element_text(size=22),legend.position="none")
g = g + xlab('Log freshwater/marine (parents)')
g = g + ylab("Log freshwater/marine (F1's)")
#g = g + xlim(-10,10)
#g = g + ylim(-10,10)
g = g + scale_color_brewer(palette="Set1")
g = g + theme_classic()
g = g + theme(legend.position = "none")
g = g + theme(legend.position="none",axis.text.x = element_text(size=20,colour = "black"),axis.text.y = element_text(size=20,colour = "black"), text = element_text(size=20,colour = "black"),strip.background = element_blank())
#g = g + scale_color_manual(values=c(marker$color[4]))
g

# save as: scatterplotAseReadCountsEffects363.pdf

# load
effectTable = read.table('AseReadCountsEffects172_bestQvalue_cisParentTrans_STAR_duprem_nosexFDR10.txt',header=T,stringsAsFactors=F)

# fold change plot
effectTable$colClass = factor(effectTable$colClass,levels=c('cis','trans','cis+trans','cis-trans','compensatory','ambiguous','conserved'))

marker = list(color = brewer.pal(5, "Set1"))

Pplot = effectTable[effectTable[,'colClass']!='ambiguous' & effectTable[,'colClass']!='conserved',]
#Pplot = effectTable[effectTable[,'colClass']=='cis',]
#Pplot = Pplot[which(sign(Pplot$P.fm.m) == sign(Pplot$F1.fm.m)),]
#Pplot = effectTable[effectTable[,'colClass']=='trans',]
#Pplot = effectTable[effectTable[,'colClass']=='cis+trans',]
#Pplot = effectTable[effectTable[,'colClass']=='cis-trans',]
#Pplot = effectTable[effectTable[,'colClass']=='compensatory',]

g = ggplot(Pplot)
g = g + geom_point(aes(foldchange2logratio(P.fm.m),foldchange2logratio(F1.fm.m),colour=colClass),cex=2,alpha=0.7)
g = g + theme(panel.background=element_rect(fill='white',colour='black'),text=element_text(size=22),legend.position="none")
g = g + xlab('Log freshwater/marine (parents)')
g = g + ylab("Log freshwater/marine (F1's)")
#g = g + xlim(-10,10)
#g = g + ylim(-10,10)
g = g + scale_color_brewer(palette="Set1")
g = g + theme_classic()
g = g + theme(legend.position = "none")
g = g + theme(legend.position="none",axis.text.x = element_text(size=20,colour = "black"),axis.text.y = element_text(size=20,colour = "black"), text = element_text(size=20,colour = "black"),strip.background = element_blank())
#g = g + scale_color_manual(values=c(marker$color[4]))
g

# save as: scatterplotAseReadCountsEffects172.pdf


# load
effectTable = read.table('AseReadCountsEffects212_bestQvalue_cisParentTrans_STAR_duprem_nosexFDR10.txt',header=T,stringsAsFactors=F)

# fold change plot
effectTable$colClass = factor(effectTable$colClass,levels=c('cis','trans','cis+trans','cis-trans','compensatory','ambiguous','conserved'))

marker = list(color = brewer.pal(5, "Set1"))

Pplot = effectTable[effectTable[,'colClass']!='ambiguous' & effectTable[,'colClass']!='conserved',]
#Pplot = effectTable[effectTable[,'colClass']=='cis',]
#Pplot = Pplot[which(sign(Pplot$P.fm.m) == sign(Pplot$F1.fm.m)),]
#Pplot = effectTable[effectTable[,'colClass']=='trans',]
#Pplot = effectTable[effectTable[,'colClass']=='cis+trans',]
#Pplot = effectTable[effectTable[,'colClass']=='cis-trans',]
#Pplot = effectTable[effectTable[,'colClass']=='compensatory',]

g = ggplot(Pplot)
g = g + geom_point(aes(foldchange2logratio(P.fm.m),foldchange2logratio(F1.fm.m),colour=colClass),cex=2,alpha=0.7)
g = g + theme(panel.background=element_rect(fill='white',colour='black'),text=element_text(size=22),legend.position="none")
g = g + xlab('Log freshwater/marine (parents)')
g = g + ylab("Log freshwater/marine (F1's)")
#g = g + xlim(-10,10)
#g = g + ylim(-10,10)
g = g + scale_color_brewer(palette="Set1")
g = g + theme_classic()
g = g + theme(legend.position = "none")
g = g + theme(legend.position="none",axis.text.x = element_text(size=20,colour = "black"),axis.text.y = element_text(size=20,colour = "black"), text = element_text(size=20,colour = "black"),strip.background = element_blank())
#g = g + scale_color_manual(values=c(marker$color[4]))
g


# save as: scatterplotAseReadCountsEffects212.pdf


# load
effectTable = read.table('AseReadCountsEffects214_bestQvalue_cisParentTrans_STAR_duprem_nosexFDR10.txt',header=T,stringsAsFactors=F)

# fold change plot
effectTable$colClass = factor(effectTable$colClass,levels=c('cis','trans','cis+trans','cis-trans','compensatory','ambiguous','conserved'))

marker = list(color = brewer.pal(5, "Set1"))

Pplot = effectTable[effectTable[,'colClass']!='ambiguous' & effectTable[,'colClass']!='conserved',]
#Pplot = effectTable[effectTable[,'colClass']=='cis',]
#Pplot = Pplot[which(sign(Pplot$P.fm.m) == sign(Pplot$F1.fm.m)),]
#Pplot = effectTable[effectTable[,'colClass']=='trans',]
#Pplot = effectTable[effectTable[,'colClass']=='cis+trans',]
#Pplot = effectTable[effectTable[,'colClass']=='cis-trans',]
#Pplot = effectTable[effectTable[,'colClass']=='compensatory',]

g = ggplot(Pplot)
g = g + geom_point(aes(foldchange2logratio(P.fm.m),foldchange2logratio(F1.fm.m),colour=colClass),cex=2,alpha=0.7)
g = g + theme(panel.background=element_rect(fill='white',colour='black'),text=element_text(size=22),legend.position="none")
g = g + xlab('Log freshwater/marine (parents)')
g = g + ylab("Log freshwater/marine (F1's)")
#g = g + xlim(-10,10)
#g = g + ylim(-10,10)
g = g + scale_color_brewer(palette="Set1")
g = g + theme_classic()
g = g + theme(legend.position = "none")
g = g + theme(legend.position="none",axis.text.x = element_text(size=20,colour = "black"),axis.text.y = element_text(size=20,colour = "black"), text = element_text(size=20,colour = "black"),strip.background = element_blank())
#g = g + scale_color_manual(values=c(marker$color[4]))
g


# save as: scatterplotAseReadCountsEffects214.pdf

