library(gplots)
library(ggplot2)
library(gdata)
library(reshape)
library(gtools)
library(RColorBrewer)

# define xloc with the same sign of expression between freshwater and marine parents
normCounts = read.table('.../xlocCountTableVarTransformed.txt')

mar363 = c(
'c363_P_FC08_F_0',
'c358_P_FC09_M_0',
'c357_P_FC06_M_0',
'c353_P_FC05_F_0')
fw363 = c(
'c363_P_FC18_M_0',
'c358_P_FC12_F_0',
'c357_P_FC14_F_0',
'c353_P_FC15_M_0')

mar172 = c(
'c172_P_532_F_lane8_0',
'c169_P_342_M_0',
'c208_P_531_M_0',
'c209_P_341_F_0')
fw172 = c(
'c172_P_533_M_lane8_0',
'c169_P_432_F_0',
'c208_P_321_F_0',
'c209_P_422_M_0')

normLR363 = foldchange2logratio(foldchange(apply(normCounts[,mar363],1,mean,na.rm=T),apply(normCounts[,fw363],1,mean,na.rm=T)))
normLR172 = foldchange2logratio(foldchange(apply(normCounts[,mar172],1,mean,na.rm=T),apply(normCounts[,fw172],1,mean,na.rm=T)))
normLR = cbind(normLR363,normLR172)

# define xloc with the same sign of expression between freshwater and marine parents
signLR = sign(normLR)
sameSignLR = normLR[which(abs(apply(signLR[,c("normLR363","normLR172")],1,sum)) == 2), ]

################
# Composite PC #
################

PC2 = read.table('.../Tyne_Litc_PC2_loadings.txt')
PC5 = read.table('.../Tyne_Litc_PC5_loadings.txt')

PC = data.frame(x=0.145*PC2$x+0.063*-1*PC5$x)
rownames(PC) = rownames(PC2)
extremePCglobal = rownames(PC)[which(PC$x<quantile(PC$x,probs=c(0.01)) | PC$x>quantile(PC$x,probs=c(0.99)))]

#####################################################################
# load effect tables (combined F1's, best effect per xloc in each F1)
# reorder classes
# calculate summary freqs for effect classes
effects172 = read.table('.../AseReadCountsEffects172_bestQvalue_cisParentTrans_STAR_duprem_nosexFDR10.txt',header=T,stringsAsFactors=F)
effects172$colClass = factor(effects172$colClass,levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous'))
colList172=lapply(split(effects172,effects172$source),'[[','colClass')
freqList172 = lapply(colList172,table)
freqList172PCT = mapply("/",freqList172,lapply(freqList172,sum),SIMPLIFY = FALSE)
# assign compPC to xloc
effects172$compPC = PC[effects172$xloc,'x']

effects363 = read.table('.../AseReadCountsEffects363_bestQvalue_cisParentTrans_STAR_duprem_nosexFDR10.txt',header=T,stringsAsFactors=F)
effects363$colClass = factor(effects363$colClass,levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous'))
colList363=lapply(split(effects363,effects363$source),'[[','colClass')
freqList363 = lapply(colList363,table)
freqList363PCT = mapply("/",freqList363,lapply(freqList363,sum),SIMPLIFY = FALSE)
# assign compPC to xloc
effects363$compPC = PC[effects363$xloc,'x']


######################
# outliers in compPC #
######################

# c172
extreme172PC = effects172[which(effects172$xloc %in% extremePCglobal),]
extreme172PC = extreme172PC[which(extreme172PC$xloc %in% rownames(sameSignLR)),]
colListPCEx172=lapply(split(extreme172PC,extreme172PC$source),'[[','colClass')
freqListPCEx172 = lapply(colListPCEx172,table)
freqListPCEx172PCT = mapply("/",freqListPCEx172,lapply(freqListPCEx172,sum),SIMPLIFY = FALSE)

# c363
extreme363PC = effects363[which(effects363$xloc %in% extremePCglobal),]
extreme363PC = extreme363PC[which(extreme363PC$xloc %in% rownames(sameSignLR)),]
colListPCEx363=lapply(split(extreme363PC,extreme363PC$source),'[[','colClass')
freqListPCEx363 = lapply(colListPCEx363,table)
freqListPCEx363PCT = mapply("/",freqListPCEx363,lapply(freqListPCEx363,sum),SIMPLIFY = FALSE)

plot172 = data.frame(class=melt(freqListPCEx172PCT)$Var.1,value=melt(freqListPCEx172PCT)$value - melt(freqList172PCT)$value, F1=melt(freqListPCEx172PCT)$L1)
plot172$class = factor(plot172$class,levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous'))
plot172 = plot172[which(plot172$class %in% c('cis','trans','cis+trans','cis-trans','compensatory','conserved')),]

plot363 = data.frame(class=melt(freqListPCEx363PCT)$Var.1,value=melt(freqListPCEx363PCT)$value - melt(freqList363PCT)$value, F1=melt(freqListPCEx363PCT)$L1)
plot363$class = factor(plot363$class,levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous'))
plot363 = plot363[which(plot363$class %in% c('cis','trans','cis+trans','cis-trans','compensatory','conserved')),]

# random list of effects...
rand172=data.frame(class=melt(freqListPCEx172PCT)$Var.1,F1=melt(freqListPCEx172PCT)$L1)
rand363=data.frame(class=melt(freqListPCEx363PCT)$Var.1,F1=melt(freqListPCEx363PCT)$L1)

for (i in 3:1003){
    random363PC = effects363[sample(rownames(effects363),nrow(extreme363PC)),]
    colListPCRand363=lapply(split(random363PC,random363PC$source),'[[','colClass')
    freqListPCRand363 = lapply(colListPCRand363,table)
    freqListPCRand363PCT = mapply("/",freqListPCRand363,lapply(freqListPCRand363,sum),SIMPLIFY = FALSE)
    
    random172PC = effects172[sample(rownames(effects172),nrow(extreme172PC)),]
    colListPCRand172=lapply(split(random172PC,random172PC$source),'[[','colClass')
    freqListPCRand172 = lapply(colListPCRand172,table)
    freqListPCRand172PCT = mapply("/",freqListPCRand172,lapply(freqListPCRand172,sum),SIMPLIFY = FALSE)
    
    plot172 = data.frame(class=melt(freqListPCRand172PCT)$Var.1,value=melt(freqListPCRand172PCT)$value - melt(freqList172PCT)$value, F1=melt(freqListPCRand172PCT)$L1)
    plot172$class = factor(plot172$class,levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous'))
    plot363 = data.frame(class=melt(freqListPCRand363PCT)$Var.1,value=melt(freqListPCRand363PCT)$value - melt(freqList363PCT)$value, F1=melt(freqListPCRand363PCT)$L1)
    plot363$class = factor(plot363$class,levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous'))
    
    rand172[,i]=plot172$value
    rand363[,i]=plot363$value
    
}

plotRand = gdata::combine(melt(rand172),melt(rand363),names=c('Tyne','Litc'))
plotRand$class = factor(plotRand$class,levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous'))
plotRand = plotRand[which(plotRand$class != 'ambiguous'),]

# observed
plot172 = data.frame(class=melt(freqListPCEx172PCT)$Var.1,value=melt(freqListPCEx172PCT)$value - melt(freqList172PCT)$value, F1=melt(freqListPCEx172PCT)$L1)
plot172$class = factor(plot172$class,levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous'))
plot172 = plot172[which(plot172$class %in% c('cis','trans','cis+trans','cis-trans','compensatory','conserved')),]
plot363 = data.frame(class=melt(freqListPCEx363PCT)$Var.1,value=melt(freqListPCEx363PCT)$value - melt(freqList363PCT)$value, F1=melt(freqListPCEx363PCT)$L1)
plot363$class = factor(plot363$class,levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous'))
plot363 = plot363[which(plot363$class %in% c('cis','trans','cis+trans','cis-trans','compensatory','conserved')),]
plotAll = gdata::combine(plot172,plot363,names=c('Tyne','Litc'))

# plot
x11(height=3.5,width=3)
cols = brewer.pal(6,"Set1")
cols[6] = "#999999"
ggplot() + geom_boxplot(data=plotRand,aes(source,value),colour="grey48",outlier.alpha=0) + geom_jitter(data=plotAll,aes(source,value,colour=as.factor(class)),height=0,width=0.1,size=2) + scale_colour_manual(values=cols) + theme(axis.text.x = element_text(angle=45, hjust=1,size=24,colour='black'),axis.title.x=element_blank(),panel.background=element_rect(fill='white',colour='grey78'),panel.grid.major=element_line(colour='white'),text=element_text(size=22),legend.position="none")+ylim(-0.5,0.5) + facet_wrap(~class)
ggsave('.../bestEffectsTyneLitcGlobalPCsameSignExcessTyneLitc_withRandomDist.ps')



###############
# Parallel DE #
###############

#-----------------------------
# parallel DE and same sign LR
#-----------------------------
qThreshold = 0.2

#--------------------------------
# DE test results - TYNE
tyneDEtest = read.table('.../cuffdiff/STAR/fw_vs_mar_Tyne_UCSCgtf/gene_exp.diff',header=T,stringsAsFactors=F,sep='\t')
tyneDEsig = tyneDEtest[which(tyneDEtest$q_value<qThreshold),]
tyneDE = tyneDEsig$gene_id
tyneFC = tyneDEtest$log2.fold_change.
names(tyneFC) = tyneDEtest$gene_id
tyneLR = foldchange2logratio(tyneFC)

#--------------------------------
# DE test results - LITC
litcDEtest = read.table('.../cuffdiff/STAR/fw_vs_mar_Litc_UCSCgtf/gene_exp.diff',header=T,stringsAsFactors=F,sep='\t')
litcDEsig = litcDEtest[which(litcDEtest$q_value<qThreshold),]
litcDE = litcDEsig$gene_id
litcFC = litcDEtest$log2.fold_change.
names(litcFC) = litcDEtest$gene_id
litcLR = foldchange2logratio(litcFC)

#--------------------------------
# DE test results - both together
parallelDEtest = read.table('.../cuffdiff/STAR/fw_vs_mar_Litc_Tyne_UCSCgtf/gene_exp.diff',header=T,stringsAsFactors=F,sep='\t')
rownames(parallelDEtest) = parallelDEtest$gene_id
parallelDEsig = parallelDEtest[which(parallelDEtest$q_value<qThreshold),]
parallelDE = parallelDEsig$gene_id

#----------------------------------
# parallel
parallel = parallelDE[which(sign(litcFC[parallelDE]) == sign(tyneFC[parallelDE]))]


# classes assigned to parallel DE genes
DE172 = effects172[which(effects172$xloc %in% parallelDE),]
DE172$colClass = factor(DE172$colClass,levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous'))
colListDE172=lapply(split(DE172,DE172$source),'[[','colClass')
freqListDE172 = lapply(colListDE172,table)
freqListDE172PCT = mapply("/",freqListDE172,lapply(freqListDE172,sum),SIMPLIFY = FALSE)

DE363 = effects363[which(effects363$xloc %in% parallelDE),]
DE363$colClass = factor(DE363$colClass,levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous'))
colListDE363=lapply(split(DE363,DE363$source),'[[','colClass')
freqListDE363 = lapply(colListDE363,table)
freqListDE363PCT = mapply("/",freqListDE363,lapply(freqListDE363,sum),SIMPLIFY = FALSE)

plot172 = data.frame(class=melt(freqListDE172PCT)$Var.1,value=melt(freqListDE172PCT)$value - melt(freqList172PCT)$value, F1=melt(freqListDE172PCT)$L1)
plot172$class = factor(plot172$class,levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous'))
plot172 = plot172[which(plot172$class %in% c('cis','trans','cis+trans','cis-trans','compensatory','conserved')),]
plot363 = data.frame(class=melt(freqListDE363PCT)$Var.1,value=melt(freqListDE363PCT)$value - melt(freqList363PCT)$value, F1=melt(freqListDE363PCT)$L1)
plot363$class = factor(plot363$class,levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous'))
plot363 = plot363[which(plot363$class %in% c('cis','trans','cis+trans','cis-trans','compensatory','conserved')),]

# random list of effects...
rand172=data.frame(class=melt(freqListDE172PCT)$Var.1,F1=melt(freqListDE172PCT)$L1)
rand363=data.frame(class=melt(freqListDE363PCT)$Var.1,F1=melt(freqListDE363PCT)$L1)

for (i in 3:1003){
    random363DE = effects363[sample(rownames(effects363),nrow(DE363)),]
    colListDERand363=lapply(split(random363DE,random363DE$source),'[[','colClass')
    freqListDERand363 = lapply(colListDERand363,table)
    freqListDERand363PCT = mapply("/",freqListDERand363,lapply(freqListDERand363,sum),SIMPLIFY = FALSE)
    
    random172DE = effects172[sample(rownames(effects172),nrow(DE172)),]
    colListDERand172=lapply(split(random172DE,random172DE$source),'[[','colClass')
    freqListDERand172 = lapply(colListDERand172,table)
    freqListDERand172PCT = mapply("/",freqListDERand172,lapply(freqListDERand172,sum),SIMPLIFY = FALSE)
    
    plot172 = data.frame(class=melt(freqListDERand172PCT)$Var.1,value=melt(freqListDERand172PCT)$value - melt(freqList172PCT)$value, F1=melt(freqListDERand172PCT)$L1)
    plot172$class = factor(plot172$class,levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous'))
    plot363 = data.frame(class=melt(freqListDERand363PCT)$Var.1,value=melt(freqListDERand363PCT)$value - melt(freqList363PCT)$value, F1=melt(freqListDERand363PCT)$L1)
    plot363$class = factor(plot363$class,levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous'))
    
    rand172[,i]=plot172$value
    rand363[,i]=plot363$value
    
}

plotRand = gdata::combine(melt(rand172),melt(rand363),names=c('Tyne','Litc'))
plotRand$class = factor(plotRand$class,levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous'))
plotRand = plotRand[which(plotRand$class != 'ambiguous'),]
plotRand$source = factor(plotRand$source,levels=c('Tyne','Litc'))

# observed
plot172 = data.frame(class=melt(freqListDE172PCT)$Var.1,value=melt(freqListDE172PCT)$value - melt(freqList172PCT)$value, F1=melt(freqListDE172PCT)$L1)
plot172$class = factor(plot172$class,levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous'))
plot172 = plot172[which(plot172$class %in% c('cis','trans','cis+trans','cis-trans','compensatory','conserved')),]
plot363 = data.frame(class=melt(freqListDE363PCT)$Var.1,value=melt(freqListDE363PCT)$value - melt(freqList363PCT)$value, F1=melt(freqListDE363PCT)$L1)
plot363$class = factor(plot363$class,levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved','ambiguous'))
plot363 = plot363[which(plot363$class %in% c('cis','trans','cis+trans','cis-trans','compensatory','conserved')),]
plotAll = gdata::combine(plot172,plot363,names=c('Tyne','Litc'))

x11(height=7,width=6)
cols = brewer.pal(6,"Set1")
cols[6] = "#999999"
ggplot() + geom_boxplot(data=plotRand,aes(source,value),colour="grey48",outlier.alpha=0) + geom_jitter(data=plotAll,aes(source,value,colour=as.factor(class)),height=0,width=0.1,size=4) + scale_colour_manual(values=cols) + theme(axis.text.x = element_text(angle=45, hjust=1,size=24,colour='black'),axis.title.x=element_blank(),panel.background=element_rect(fill='white',colour='grey78'),panel.grid.major=element_line(colour='white'),text=element_text(size=22),legend.position="none")+ylim(-0.5,0.5) + facet_wrap(~class)
ggsave('.../bestEffectsCisParentTransTyneLitcDEexcessSameSignTyneLitc_withRandomDist.ps')





