library(reshape)
library(ggplot2)
library(zoo)
library(gtools)


########
# LITC #
########

# count table ==> expressed genes
normCounts = read.table('.../xlocCountTableVarTransformed.txt')
sel = c(
'c363_P_FC08_F_0',
'c363_P_FC18_M_0',
'c358_P_FC12_F_0',
'c358_P_FC09_M_0',
'c357_P_FC14_F_0',
'c357_P_FC06_M_0',
'c353_P_FC05_F_0',
'c353_P_FC15_M_0')
expressedXloc = rownames(normCounts)[which(apply(normCounts[,sel],1,mean)>9)] # expressed loci having at least as high mean expression as the target loci --> no effect of expression level possible

# stats for each xloc
xlocFst = read.table('.../xlocLitcFstTssmatrix500kb.txt')
xlocPiMar = read.table('.../xlocLitcPiMarTssmatrix500kb.txt')
xlocPiFw = read.table('.../xlocLitcPiFwTssmatrix500kb.txt')

########
# TYNE #
########

# count table ==> expressed genes
normCounts = read.table('.../xlocCountTableVarTransformed.txt')
sel = c(
'c172_P_532_F_lane8_0',
'c172_P_533_M_lane8_0')
expressedXloc = rownames(normCounts)[which(apply(normCounts[,sel],1,mean)>9)] # expressed loci having at least as high mean expression as the target loci --> no effect of expression level possible

# stats for each xloc
xlocFst = read.table('.../xlocTyneFstTssmatrix500kb.txt')
xlocPiMar = read.table('.../xlocTynePiMarTssmatrix500kb.txt')
xlocPiFw = read.table('.../xlocTynePiFwTssmatrix500kb.txt')


##############################
#-----------------------------
# composite PC2 + PC5
# 1% PC2 outliers specificly in Tyne
# outlier loci on composite PC
PC2 = read.table('.../Tyne_Litc_PC2_loadings.txt')
PC5 = read.table('.../Tyne_Litc_PC5_loadings.txt')

PC = data.frame(x=0.145*PC2$x+0.063*-1*PC5$x)
rownames(PC) = rownames(PC2)
extremePCglobal = rownames(PC)[which(PC$x<quantile(PC$x,probs=c(0.01)) | PC$x>quantile(PC$x,probs=c(0.99)))]

##############################
#-----------------------------
# subset to target and control
targetFst = xlocFst[which(rownames(xlocFst) %in% extremePCglobal),]
controlFst = xlocFst[which(rownames(xlocFst) %in% extremePCglobal == F),]

targetPiMar = xlocPiMar[which(rownames(xlocPiMar) %in% extremePCglobal),]
controlPiMar = xlocPiMar[which(rownames(xlocPiMar) %in% extremePCglobal == F),]

targetPiFw = xlocPiFw[which(rownames(xlocPiFw) %in% extremePCglobal),]
controlPiFw = xlocPiFw[which(rownames(xlocPiFw) %in% extremePCglobal == F),]

##############################
#-----------------------------
# parallel DE and same sign LR
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

# subset to target and control
targetFst = xlocFst[which(rownames(xlocFst) %in% parallel),]
controlFst = xlocFst[which(rownames(xlocFst) %in% parallel == F),]

targetPiMar = xlocPiMar[which(rownames(xlocPiMar) %in% parallel),]
controlPiMar = xlocPiMar[which(rownames(xlocPiMar) %in% parallel == F),]

targetPiFw = xlocPiFw[which(rownames(xlocPiFw) %in% parallel),]
controlPiFw = xlocPiFw[which(rownames(xlocPiFw) %in% parallel == F),]

############################
#---------------------------
# plotting function - SIMPLE
dataPlot = function(plotTarget,cand=FALSE){
    g = ggplot(plotTarget,aes(1:length(value),value))
    g = g + geom_point()
    g = g + geom_smooth(method='loess')
    if (cand==F){
        g = g + scale_x_continuous(breaks=c(50,150,250,350,450),labels=c("-200kb", "-100kb","TSS","100kb","200kb"))
    }
    else {
        g = g + scale_x_continuous(breaks=c(50,150,250,350,450)/5,labels=c("-200kb", "-100kb","TSS","100kb","200kb"))
    }
    g = g + theme(panel.background=element_rect(fill='white',colour='grey78'),panel.grid.major=element_line(colour='grey78'),legend.position="right",text=element_text(size=16))
    g
}

# plotting function - TARGET and CONTROL
dataPlot2 = function(plotTarget,plotControl,pi=FALSE){
    g = ggplot(plotTarget,aes(1:length(value),value))
    g = g + geom_point()
    g = g + geom_smooth(method='loess')
    g = g + geom_smooth(data=plotControl,aes(1:length(value),value),method='loess',color='grey',linetype='dotted')
    g = g + scale_x_continuous(breaks=c(50,150,250),labels=c("-100kb","TSS","100kb"))
    if (pi==T){
        g = g + ylim(values=c(0.0013,0.0022))
    }
    g = g + theme(panel.background=element_rect(fill='white',colour='grey78'),panel.grid.major=element_line(colour='grey78'),legend.position="right",text=element_text(size=16))
    g
}

# Fst plot TARGET
plotTarget = melt(apply(targetFst[,1:499],2,mean,na.rm=T))
dataPlot(plotTarget)
#ggsave('extremePC2LitcFstTssplot.pdf')

# Pi Mar plot TARGET
plotTarget = melt(apply(targetPiMar[,1:499],2,mean,na.rm=T))
dataPlot(plotTarget)
#ggsave('extremePC2LitcPiMarTssplot.pdf')

# Pi Fw plot TARGET
plotTarget = melt(apply(targetPiFw[,1:499],2,mean,na.rm=T))
dataPlot(plotTarget)
#ggsave('extremePC2LitcPiFwTssplot.pdf')

# Fst plot CANDS
#plotTarget = melt(rollapply(apply(targetFst[cands,],2,mean,na.rm=T),width=10,by=5,FUN=mean))
#dataPlot(plotTarget,cand=TRUE)

# Pi Mar plot CANDS
#plotTarget = melt(rollapply(apply(targetPiMar[cands,],2,mean,na.rm=T)/1000,width=10,by=5,FUN=mean))
#dataPlot(plotTarget,cand=T)

# Pi Fw plot CANDS
#plotTarget = melt(rollapply(apply(targetPiFw[cands,],2,mean,na.rm=T)/1000,width=10,by=5,FUN=mean))
#dataPlot(plotTarget,cand=T)

# Fst plot CONTROL
plotControl = melt(apply(controlFst[,1:499],2,mean,na.rm=T))
dataPlot(plotControl)

# Pi Mar plot CONTROL
plotControl = melt(apply(controlPiMar[,1:499],2,mean,na.rm=T)/1000)
dataPlot(plotControl)

# Fst plot TARGET + CONTROL
plotControl = melt(apply(controlFst[,100:400],2,mean,na.rm=T))
plotTarget = melt(apply(targetFst[,100:400],2,mean,na.rm=T))
dataPlot2(plotTarget,plotControl)
ggsave('.../parallelPCAtyneFstTssplotWithControl.ps')

# Pi Mar plot TARGET + CONTROL
plotControl = melt(apply(controlPiMar[,100:400]/100,2,mean,na.rm=T))
plotTarget = melt(apply(targetPiMar[,100:400]/100,2,mean,na.rm=T))
dataPlot2(plotTarget,plotControl,pi=TRUE)
ggsave('.../parallelPCAtynePiMarTssplotWithControl.ps')

# Pi Fw plot TARGET + CONTROL
plotControl = melt(apply(controlPiFw[,100:400],2,mean,na.rm=T)/100)
plotTarget = melt(apply(targetPiFw[100:400],2,mean,na.rm=T)/100)
dataPlot2(plotTarget,plotControl,pi=TRUE)
ggsave('.../parallelPCAtynePiFWTssplotWithControl.ps')



