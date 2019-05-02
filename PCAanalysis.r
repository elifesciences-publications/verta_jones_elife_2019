library(DESeq2)
library(extrafont)
library(ggplot2)
library(reshape)
library(ggbiplot)

setwd('.../expression/cuffnorm/STAR/UCSCgtf')
counts = read.table('genes.count_table',row.names=1,header=T,sep='\t')

# only pure ecotypes
sel = c(
'c172_P_532_F_lane8_0',
'c172_P_533_M_lane8_0',
'c172_P_532_F_lane3_0',
'c172_P_533_M_lane3_0',

'c169_P_432_F_0',
'c169_P_342_M_0',
'c208_P_321_F_0',
'c208_P_531_M_0',
'c209_P_341_F_0',
'c209_P_422_M_0',

'c212_P_454_F_0',
'c212_P_551_M_0',

'c214_P_524_F_0',
'c214_P_512_M_0',

'c363_P_FC08_F_0',
'c363_P_FC18_M_0',

'c358_P_FC12_F_0',
'c358_P_FC09_M_0',
'c357_P_FC14_F_0',
'c357_P_FC06_M_0',
'c353_P_FC05_F_0',
'c353_P_FC15_M_0')

river = c(rep('Tyne',10),rep('Forss',2),rep('Shiel',2),rep('Litc',8))

eco = c(

'M','F','M','F',

'F',
'M',
'F',
'M',
'M',
'F',

'M',
'F',

'M',
'F',

'M',
'F',

'F',
'M',
'F',
'M',
'M',
'F')

# normalize using DESeq2
subCounts = counts[,sel]
design = data.frame(row.names=sel,river = river,eco = eco)

iCounts = round(subCounts)
iCounts = DESeqDataSetFromMatrix(countData=iCounts, colData=design, design=~river+eco)
vCounts = varianceStabilizingTransformation(iCounts,blind=FALSE)
normCounts = assay(vCounts,normalized=T)

# save
write.table(normCounts,'.../xlocCountTableVarTransformed.txt')
# load
normCounts = read.table('.../xlocCountTableVarTransformed.txt')


#######
# PCA #
#######

# only Tyne and Litc pure strains
sel = c('c172_P_532_F_lane3_0',
'c172_P_533_M_lane3_0','c169_P_432_F_0',
'c169_P_342_M_0',
'c208_P_321_F_0',
'c208_P_531_M_0',
'c209_P_341_F_0',
'c209_P_422_M_0',
'c363_P_FC08_F_0',
'c363_P_FC18_M_0',
'c358_P_FC12_F_0',
'c358_P_FC09_M_0',
'c357_P_FC14_F_0',
'c357_P_FC06_M_0',
'c353_P_FC05_F_0',
'c353_P_FC15_M_0')

subNormCounts = normCounts[,sel]
subDesign = design[sel,]
subDesign$sex = c('Female','Male','Female','Male','Female','Male','Female','Male','Female','Male','Female','Male','Female','Male','Female','Male')

# PC2 versus PC5
#postscript('.../PC2PC5plotTyneLitc.ps')
pca = prcomp(t(subNormCounts),scale=F)
data <- cbind(pca$x,subDesign)
percentVar <- round(100 * pca$sdev^2/sum(pca$sdev^2))
ggplot(data, aes(PC2, PC5, color=eco, shape=river)) +
scale_color_manual(values=c('darkblue','darkred')) +
geom_point(size=8,stroke=2) +
xlab(paste0("PC2: ",percentVar[2],"% variance")) +
ylab(paste0("PC5: ",percentVar[5],"% variance")) +
scale_shape_manual(values=c(0,1))+
theme(axis.text=element_text(hjust=1,size=24,colour='black'),panel.background=element_rect(fill='white',colour='black'),panel.grid.major=element_line(colour='grey78'),text=element_text(size=22),legend.position="right")+
coord_fixed()
dev.off()

# project forss&shiel parents to PCA space --> PC 2 and 5
subCountsFS = normCounts[,which(colnames(normCounts) %in% c("c212_P_454_F_0","c212_P_551_M_0","c214_P_524_F_0","c214_P_512_M_0"))]
FSproject = data.frame(scale(t(subCountsFS), pca$center, pca$scale) %*% pca$rotation[,c(2,5)],
river = c('Forss','Forss','Shiel','Shiel'),
eco = c('M','F','M','F'))
# re=plot
#postscript('.../PC2PC5plotTyneLitcForssShielParentsProjected.ps')
percentVar <- round(100 * pca$sdev^2/sum(pca$sdev^2))
g = ggbiplot(pca, var.axes=F, choices=c(2,5), obs.scale = 1, var.scale = 1, groups = as.factor(paste(subDesign$eco,subDesign$river)), ellipse = FALSE, ellipse.prob=0.68, alpha=0)
g = g + scale_color_manual(values=c('darkblue','darkblue','darkblue','darkred','darkred','darkred'))
g = g + theme(axis.text=element_text(hjust=1,size=24,colour='black'),panel.background=element_rect(fill='white',colour='black'),panel.grid.major=element_line(colour='grey78'),text=element_text(size=22),legend.position="right")
g = g + geom_point(data=data, aes(PC2, PC5, color=eco, shape=river),size=8)
g = g + geom_point(data=FSproject,aes(PC2, PC5,color=eco, shape=river),size=8)
g = g + scale_shape_manual(values=c(0,16,1,18))
g
dev.off()

# add additional axis corresponding to composite PC2 + PC5 (invert PC5 sign - FW and MAR associated in opposite directions on PC2 and PC5)
data$compPC = 0.145*data$PC2+0.063*-1*data$PC5

# plot the projection of samples on composite PC2 and PC5
p25 = ggplot(data,aes(x=0,y=0.145*PC2+0.063*-1*PC5,color=eco, shape=river))+
scale_color_manual(values=c('darkblue','darkred')) +
geom_point(size=6,stroke=2)+
xlim(-5,5)+
ylim(-5,5)+
scale_shape_manual(values=c(16,18))+
theme(axis.text=element_text(hjust=1,size=16,colour='black'),panel.background=element_rect(fill='white',colour='black'),panel.grid.major=element_blank(),text=element_text(size=16),legend.position="right",axis.text.x=element_blank())+
coord_fixed(ratio = 1)
p25
ggsave('.../LitcTyneCompPCprojection.ps')

# plot the projection of samples on PC2
p2 = ggplot(data,aes(x=0,y=PC2,color=eco, shape=river))+
scale_color_manual(values=c('darkblue','darkred')) +
geom_point(size=6,stroke=2)+
xlim(-5,5)+
#ylim(-5,5)+
scale_shape_manual(values=c(16,18))+
theme(axis.text=element_text(hjust=1,size=16,colour='black'),panel.background=element_rect(fill='white',colour='black'),panel.grid.major=element_blank(),text=element_text(size=16),legend.position="right",axis.text.x=element_blank())+
coord_fixed(ratio = 1)
p2
ggsave('.../LitcTynePC2projection.ps')


# save loadings
PC2 = pca$rotation[,'PC2']
write.table(PC2,'.../Tyne_Litc_PC2_loadings.txt')

PC5 = pca$rotation[,'PC5']
write.table(PC5,'.../Tyne_Litc_PC5_loadings.txt')


#-----------------------
# expression log ratios
#-----------------------

library(gtools)
tyneFwMar = foldchange2logratio(foldchange(apply(normCounts[,which(design$river=='Tyne' & design$eco=='F')],1,mean),apply(normCounts[,which(design$river=='Tyne' & design$eco=='M')],1,mean)))

litcFwMar = foldchange2logratio(foldchange(apply(normCounts[,which(design$river=='Litc' & design$eco=='F')],1,mean),apply(normCounts[,which(design$river=='Litc' & design$eco=='M')],1,mean)))

forssFwMar = foldchange2logratio(foldchange(normCounts[,which(design$river=='Forss' & design$eco=='F')],normCounts[,which(design$river=='Forss' & design$eco=='M')]))

shielFwMar = foldchange2logratio(foldchange(normCounts[,which(design$river=='Shiel' & design$eco=='F')],normCounts[,which(design$river=='Shiel' & design$eco=='M')]))

# Litc FC versus PC2
fit = lm(litcFwMar~pca$rotation[,"PC2"])
g = ggplot(data.frame(LR=litcFwMar,PC=pca$rotation[,"PC2"]))
g = g + geom_point(aes(LR,PC),size=3,alpha=0.5)
g = g + theme(axis.text=element_text(hjust=1,size=24,colour='black'),panel.background=element_rect(fill='white',colour='black'),panel.grid.major=element_line(colour='grey78'),text=element_text(size=26),legend.position="right")
g = g + geom_smooth(aes(LR,PC),method="lm")
g = g + labs(subtitle=paste(paste('beta',round(fit$coefficients[2],digits=3),sep='='),paste('R2',round(signif(summary(fit)$adj.r.squared,digits=3), 5),sep='=')),size=12)
#g = g + ylim(-0.02,0.02)
#g = g + xlim(-0.9,0.9)
g
ggsave('.../litcFwMarPC2.pdf')

# Tyne FC versus PC5
fit = lm(tyneFwMar~pca$rotation[,"PC5"])
g = ggplot(data.frame(LR=tyneFwMar,PC=pca$rotation[,"PC5"]))
g = g + geom_point(aes(LR,PC),size=3,alpha=0.5)
g = g + theme(axis.text=element_text(hjust=1,size=24,colour='black'),panel.background=element_rect(fill='white',colour='black'),panel.grid.major=element_line(colour='grey78'),text=element_text(size=26),legend.position="right")
g = g + geom_smooth(aes(LR,PC),method="lm")
g = g + labs(subtitle=paste(paste('beta',round(fit$coefficients[2],digits=3),sep='='),paste('R2',round(signif(summary(fit)$adj.r.squared,digits=3), 5),sep='=')),size=12)
#g = g + ylim(-0.02,0.02)
#g = g + xlim(-0.9,0.9)
g
ggsave('.../tyneFwMarPC5.pdf')

# composite PC rotation
compRotation = 0.145*pca$rotation[,"PC2"]+0.063*-1*pca$rotation[,"PC5"]

# Litc FC versus compPC
fit = lm(litcFwMar~compRotation)
g = ggplot(data.frame(LR=litcFwMar,PC=compRotation))
g = g + geom_point(aes(LR,PC),size=3,alpha=0.5)
g = g + theme(axis.text=element_text(hjust=1,size=24,colour='black'),panel.background=element_rect(fill='white',colour='black'),panel.grid.major=element_line(colour='grey78'),text=element_text(size=26),legend.position="right")
g = g + geom_smooth(aes(LR,PC),method="lm")
g = g + labs(subtitle=paste(paste('beta',round(fit$coefficients[2],digits=3),sep='='),paste('R2',round(signif(summary(fit)$adj.r.squared,digits=3), 5),sep='=')),size=12)
g = g + ylim(-0.02,0.02)
g = g + xlim(-0.9,0.9)
g
ggsave('.../litcFwMarCompPC.pdf')

# Tyne FC versus compPC
fit = lm(tyneFwMar~compRotation)
g = ggplot(data.frame(LR=tyneFwMar,PC=compRotation))
g = g + geom_point(aes(LR,PC),size=3,alpha=0.5)
g = g + theme(axis.text=element_text(hjust=1,size=24,colour='black'),panel.background=element_rect(fill='white',colour='black'),panel.grid.major=element_line(colour='grey78'),text=element_text(size=26),legend.position="right")
g = g + geom_smooth(aes(LR,PC),method="lm")
g = g + labs(subtitle=paste(paste('beta',round(fit$coefficients[2],digits=3),sep='='),paste('R2',round(signif(summary(fit)$adj.r.squared,digits=3), 5),sep='=')),size=12)
g = g + ylim(-0.02,0.02)
g = g + xlim(-0.9,0.9)
g
ggsave('.../tyneFwMarCompPC.pdf')


