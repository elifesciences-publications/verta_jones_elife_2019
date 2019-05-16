library(ggplot2)
library(gdata)
library(reshape)
library(gtools)
library(RColorBrewer)


#####################################################################
# load effect tables (combined F1's, best effect per xloc in each F1)
# reorder classes
# calculate summary freqs for effect classes
effectTable172 = read.table('AseReadCountsEffects172_bestQvalue_cisParentTrans_STAR_duprem_nosexFDR10.txt',header=T,stringsAsFactors=F)
effects172 = table(effectTable172[,c("source","colClass")])

effectTable363 = read.table('AseReadCountsEffects363_bestQvalue_cisParentTrans_STAR_duprem_nosexFDR10.txt',header=T,stringsAsFactors=F)
effects363 = table(effectTable363[,c("source","colClass")])

effectTable212 = read.table('AseReadCountsEffects212_bestQvalue_cisParentTrans_STAR_duprem_nosexFDR10.txt',header=T,stringsAsFactors=F)
effects212 = table(effectTable212[,c("source","colClass")])

effectTable214 = read.table('AseReadCountsEffects214_bestQvalue_cisParentTrans_STAR_duprem_nosexFDR10.txt',header=T,stringsAsFactors=F)
effects214 = table(effectTable214[,c("source","colClass")])


# as percentage of all
combinedEffectsPCT = gdata::combine(melt(effects363/apply(effects363,1,sum)),melt(effects172/apply(effects172,1,sum)),melt(effects212/apply(effects212,1,sum)),melt(effects214/apply(effects214,1,sum)),names=c('Litc','Tyne','Forss','Shiel'))
colnames(combinedEffectsPCT) = c('ind','class','frequency','cross')
combinedEffectsPCT=combinedEffectsPCT[which(combinedEffectsPCT$class!='ambiguous'),]
combinedEffectsPCT$class = factor(combinedEffectsPCT$class,levels=c('cis','trans','cis+trans','cis-trans','compensatory','conserved'))
combinedEffectsPCT$cross = factor(combinedEffectsPCT$cross,levels=c('Tyne','Shiel','Forss','Litc'))
cols = brewer.pal(6,"Set1")
cols[6] = "#999999"
g = ggplot(combinedEffectsPCT)
g = g + geom_jitter(aes(class,frequency*100,colour=as.factor(class)),height=0,width=0.3,size=4,alpha=0.7)
g = g + facet_wrap(~cross,nrow=1)
g = g + theme_classic()
g = g + theme(axis.text.y = element_text(size=24,colour='black'), axis.text.x = element_text(angle=45, hjust=1,size=24,colour=cols),axis.title.x=element_blank(),strip.background = element_blank(),text=element_text(size=22),legend.position="none")
g = g + ylab('Frequency %')
g = g + scale_colour_manual(values=cols)
g
# save as boxplotAseReadCountsEffectsAllCrossesPercentage.pdf
