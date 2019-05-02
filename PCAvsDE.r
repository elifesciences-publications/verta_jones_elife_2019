library(gtools)
library(gplots)
library(ggplot2)

#--------------------------------
geneFeatures = read.table('.../assembly/merged_assembly_STAR_UCSCgtf/mergedParsed.txt',stringsAsFactors=F,sep='\t',skip=1)
colnames(geneFeatures) = c('chrom','source','type','start','end','V6','strand','V8','V9','gene_id','tcons','exon','name','oID','ensGene')

#--------------------------------
qThreshold = 0.2

#--------------------------------
# DE test results - TYNE
tyneDEtest = read.table('.../expression/cuffdiff/STAR/fw_vs_mar_Tyne_UCSCgtf/gene_exp.diff',header=T,stringsAsFactors=F,sep='\t')
tyneDEsig = tyneDEtest[which(tyneDEtest$q_value<qThreshold),]
tyneDE = tyneDEsig$gene_id
tyneFC = tyneDEtest$log2.fold_change.
names(tyneFC) = tyneDEtest$gene_id
tyneLR = foldchange2logratio(tyneFC)

#--------------------------------
# DE test results - LITC
litcDEtest = read.table('.../expression/cuffdiff/STAR/fw_vs_mar_Litc_UCSCgtf/gene_exp.diff',header=T,stringsAsFactors=F,sep='\t')
litcDEsig = litcDEtest[which(litcDEtest$q_value<qThreshold),]
litcDE = litcDEsig$gene_id
litcFC = litcDEtest$log2.fold_change.
names(litcFC) = litcDEtest$gene_id
litcLR = foldchange2logratio(litcFC)

#--------------------------------
# DE test results - both together
parallelDEtest = read.table('.../expression/cuffdiff/STAR/fw_vs_mar_Litc_Tyne_UCSCgtf/gene_exp.diff',header=T,stringsAsFactors=F,sep='\t')
rownames(parallelDEtest) = parallelDEtest$gene_id
parallelDEsig = parallelDEtest[which(parallelDEtest$q_value<qThreshold),]
parallelDE = parallelDEsig$gene_id

#----------------------------------
# parallel
parallel = parallelDE[which(sign(litcFC[parallelDE]) == sign(tyneFC[parallelDE]))]

#----------------------------------
# composite PC
PC2 = read.table('.../Tyne_Litc_PC2_loadings.txt')
PC5 = read.table('.../Tyne_Litc_PC5_loadings.txt')

PC = data.frame(x=0.145*PC2$x+0.063*-1*PC5$x)
rownames(PC) = rownames(PC2)
extremePCglobal = rownames(PC)[which(PC$x<quantile(PC$x,probs=c(0.01)) | PC$x>quantile(PC$x,probs=c(0.99)))]
PC$zPC = scale(PC$x)

#----------------------------------
# composite PC loading versus parallel DE p-value
ggplot(data.frame(qvalue=parallelDEtest$q_value,PC=abs(PC[rownames(parallelDEtest),'zPC'])))+
geom_point(aes(qvalue,PC),fill='grey48')+
geom_point(data=data.frame(qvalue=parallelDEtest[parallel,'q_value'],PC=abs(PC[parallel,'zPC'])),aes(qvalue,PC),col='blue')+
geom_smooth(aes(qvalue,PC),method='lm')+
theme(axis.text = element_text(hjust=1,size=24,colour='black'),panel.background=element_rect(fill='white',colour='grey78'),text=element_text(size=22),legend.position="none")
ggsave('.../compositePCloadingVsQvalue.ps')

fit = lm(abs(PC[rownames(parallelDEtest),'zPC']) ~ parallelDEtest$p_value)
summary(fit)

# density plot
ggplot(data.frame(PC=PC$zPC[which(rownames(PC) %in% rownames(parallelDEtest))]))+
geom_density(aes(abs(PC)),size=1.5)+
theme(axis.text = element_text(hjust=1,size=24,colour='black'),panel.background=element_rect(fill='white',colour='grey78'),text=element_text(size=22),legend.position="none")+
coord_flip()+
geom_density(data=data.frame(PC=PC$zPC[which(rownames(PC) %in% parallel)]),aes(abs(PC)),colour='blue',size=1.5)
ggsave('.../compositePCloadingVsQvalueDensity.ps')

# overlap probability
cufflinksPCAoverlap = attr(venn(list(extremePCglobal,parallel)),"intersect")$'A:B'
total = length(which(parallelDEtest$status != 'NOTEST'))
deN = length(parallel)
pcN = length(extremePCglobal)
overlap = length(cufflinksPCAoverlap)
phyper(overlap - 1, deN, total - deN, pcN, lower.tail=F)






