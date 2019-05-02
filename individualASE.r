
#######################
# ASE for individuals #
#######################
library(qvalue)
library(gtools)

# normalized coverage - only 3.5 ppt individuals
normCounts = read.table('.../aseReadCounts_c172_F1_infSites_alleleCov_DESeqNormalized_STAR_duprem.txt',stringsAsFactors=F,header=T)

# select for SNPs that have at least 10 reads -- otherwise you're just testing for nothing and the FDR correction will be too strict -- this filtering makes a HUGE difference!!!!
normCounts = normCounts[apply(normCounts>=10,1,any),]

# ASE test -- MODIFY COVERAGE REQUIEMENT!!!
# requiring minimum coverage of 10 gives good results --  avoid unneccessary testing of low counts and too strickt FDR correction
altAseTest = function(x,prob){
    if (x[1]>10 | x[2]>10 & x[1] < c(x[1]+x[2])){
        return(unlist(binom.test(c(as.integer(x[1]),as.integer(x[2])),p=prob)['p.value']))
    }
    else {
        return(NA)
    }
}

# for individuals in 3.5 ppt
inds=c('c172_F1_20_F','c172_F1_20_M','c172_F1_04_F','c172_F1_04_M')
aseResInd3 = data.frame()
for (y in inds){
    aseResInd3[1:nrow(normCounts),y] = apply(X=data.frame(normCounts[,paste(y,'REFCOV',sep='_')],normCounts[,paste(y,'ALTCOV',sep='_')]),1,FUN=altAseTest,
    prob=sum(normCounts[,paste(y,'REFCOV',sep='_')])/(sum(normCounts[,paste(y,'REFCOV',sep='_')])+sum(normCounts[,paste(y,'ALTCOV',sep='_')])))
    
}

# qvalues of binomial test
aseResInd3[rownames(aseResInd3[is.na(aseResInd3[,1])==F,]),1] = qvalue(na.omit(aseResInd3[,1]))$qvalue
aseResInd3[rownames(aseResInd3[is.na(aseResInd3[,2])==F,]),2] = qvalue(na.omit(aseResInd3[,2]))$qvalue
aseResInd3[rownames(aseResInd3[is.na(aseResInd3[,3])==F,]),3] = qvalue(na.omit(aseResInd3[,3]))$qvalue
aseResInd3[rownames(aseResInd3[is.na(aseResInd3[,4])==F,]),4] = qvalue(na.omit(aseResInd3[,4]))$qvalue

aseResInd3Sig = aseResInd3<0.1

length(which(aseResInd3Sig[,1]))
length(which(aseResInd3Sig[,2]))
length(which(aseResInd3Sig[,3]))
length(which(aseResInd3Sig[,4]))

length(which(aseResInd3Sig[,1]))/length(na.omit(aseResInd3[,1]))
length(which(aseResInd3Sig[,2]))/length(na.omit(aseResInd3[,2]))
length(which(aseResInd3Sig[,3]))/length(na.omit(aseResInd3[,3]))
length(which(aseResInd3Sig[,4]))/length(na.omit(aseResInd3[,4]))

# verification plot
plot(log(normCounts[,'c172_F1_20_M_REFCOV']),log(normCounts[,'c172_F1_20_M_ALTCOV']),pch=20,cex=.2)
points(log(normCounts[which(aseResInd3Sig[,2]),'c172_F1_20_M_REFCOV']),log(normCounts[which(aseResInd3Sig[,2]),'c172_F1_20_M_ALTCOV']),pch=20,cex=.2,col='red')
abline(a=0,b=1,col='red')

# ven diagram of ASE sharing
library(gplots)
venn(list(which(aseResInd3Sig[,1]),which(aseResInd3Sig[,2]),which(aseResInd3Sig[,3]),which(aseResInd3Sig[,4])))
# random expectation
#venn(list(sample(1:nrow(aseResInd3Sig),length(which(aseResInd3Sig[,1]))),sample(1:nrow(aseResInd3Sig),length(which(aseResInd3Sig[,2]))),sample(1:nrow(aseResInd3Sig),length(which(aseResInd3Sig[,3]))),sample(1:nrow(aseResInd3Sig),length(which(aseResInd3Sig[,4])))))

# polarise ASE
genoInf = read.table('.../c172_P_genomic_recalibrated_filtered_SNPs_informative_minCov10_genotypes.txt',header=T,stringsAsFactors=F)
rownames(genoInf) = paste(genoInf[,1],genoInf[,2],sep=':')

#contruct a bedtrack with qvalue for each tested position --> with mean counts, all F1's need to have counts in both ref and alt positions (otherwise mean function will give NA), this is quite strickt
ind='c172_F1_20_M'

aseTableInd3 = data.frame(rownames(normCounts),aseResInd3[,ind])
aseTableInd3[,3] = foldchange(normCounts[,paste(ind,'REFCOV',sep='_')],normCounts[,paste(ind,'ALTCOV',sep='_')])
aseTableInd3[,4] = normCounts[,paste(ind,'REFCOV',sep='_')]
aseTableInd3[,5] = normCounts[,paste(ind,'ALTCOV',sep='_')]
colnames(aseTableInd3)[c(1,2,3,4,5)] = c('pos','qvalue','ref:alt','refCount','altCount')
rownames(aseTableInd3) = aseTableInd3[,1]
aseTableInd3 = na.omit(aseTableInd3)

# add polarized parent fold change to aseTable
# --> female:male ratio
# --> invert the ref:alt fold change in cases where the female is alt
aseTableInd3[,6] = aseTableInd3[,3]
fAlt = genoInf[rownames(aseTableInd3),'female_P'] == '1/1'
names(fAlt) = rownames(aseTableInd3)
aseTableInd3[names(which(fAlt)),6] = -1*(aseTableInd3[names(which(fAlt)),6])
colnames(aseTableInd3)[6] = c('fm:m')

#convert to bed format
f = function(x){c(unlist(strsplit(x[1],':')),as.numeric(unlist(strsplit(x[1],':'))[2])+1,as.numeric(x[2]),as.numeric(x[3]),as.numeric(x[4]),as.numeric(x[5]),as.numeric(x[6]))}
aseInd3 = apply(aseTableInd3,1,f)
aseBedInd3 = matrix(aseInd3,ncol=8,byrow=T)
head(aseBedInd3)
colnames(aseBedInd3) = c('chr','start','end','qvalue','ref:alt','refCount','altCount','fm:m')

write.table(aseBedInd3,paste('.../aseTable_',ind,'_FDR10_minCov10_aseReadCounts_refCorrected_F1normalized_genomicMask_STAR_duprem.txt',sep=''),row.names=F)








