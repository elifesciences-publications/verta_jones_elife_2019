
# load bed files
# define a unique set of positions
# parse individual bed file into a dataframe of all unique positions
# select lines where the is no NA's

inds = c('c172_F1_20_F', 'c172_F1_20_M', 'c172_F1_04_F', 'c172_F1_04_M', 'c172_F1_10_F', 'c172_F1_10_M', 'c172_F1_13_F', 'c172_F1_13_M', 'c172_F1_01_F', 'c172_F1_01_M', 'c172_F1_05_F', 'c172_F1_05_M','c172_P_532_F_lane3', 'c172_P_533_M_lane3', 'c172_P_532_F_lane8', 'c172_P_533_M_lane8')

cov = list()
pos = c()
for (i in inds){
	print(i)
    cov[[i]] = read.table(paste('.../SNPs/',i,'/STAR/UCSCgtf/genomicMask/',i,'_aseReadCounter_infSites_duprem.txt',sep=''),header=T,sep='\t') # use with parents!
    pos = append(pos,paste(cov[[i]][,1],cov[[i]][,2],sep=':'))
}	
#unique positions - rows in data.frame
uPos = unique(pos)

library(data.table)

#merge all coverage files letting NA's fill unmatched positions - use DATA TABLE !!!
Cov = data.frame(pos=uPos)
rownames(Cov) = Cov[,'pos']
Cov = setDT(Cov)
for (i in names(cov)){
	bed = data.table(cov[[i]])
	bed[,'pos'] = paste(cov[[i]][,1],cov[[i]][,2],sep=':')
	colnames(bed)[6] = paste(i,'REFCOV',sep='_') 
	colnames(bed)[7] = paste(i,'ALTCOV',sep='_') 
	Cov = merge(Cov,bed[,c(paste(i,'REFCOV',sep='_'),paste(i,'ALTCOV',sep='_'),'pos'),with=F],all=T,by='pos')
}	

# check that parents are indeed fully informative
hist(abs(Cov$c172_P_532_F_lane8_REFCOV-Cov$c172_P_532_F_lane8_ALTCOV)/(Cov$c172_P_532_F_lane8_REFCOV+Cov$c172_P_532_F_lane8_ALTCOV))
hist(abs(Cov$c172_P_532_F_lane3_REFCOV-Cov$c172_P_532_F_lane3_ALTCOV)/(Cov$c172_P_532_F_lane3_REFCOV+Cov$c172_P_532_F_lane3_ALTCOV))

# filter out SNPs with heterozygous counts in the parents
hetSNPs = which(abs(Cov$c172_P_532_F_lane3_REFCOV-Cov$c172_P_532_F_lane3_ALTCOV)/(Cov$c172_P_532_F_lane3_REFCOV+Cov$c172_P_532_F_lane3_ALTCOV)<0.99 |
                abs(Cov$c172_P_532_F_lane8_REFCOV-Cov$c172_P_532_F_lane8_ALTCOV)/(Cov$c172_P_532_F_lane8_REFCOV+Cov$c172_P_532_F_lane8_ALTCOV)<0.99 |
                abs(Cov$c172_P_533_M_lane3_REFCOV-Cov$c172_P_533_M_lane3_ALTCOV)/(Cov$c172_P_533_M_lane3_REFCOV+Cov$c172_P_533_M_lane3_ALTCOV)<0.99 |
                abs(Cov$c172_P_533_M_lane8_REFCOV-Cov$c172_P_533_M_lane8_ALTCOV)/(Cov$c172_P_533_M_lane8_REFCOV+Cov$c172_P_533_M_lane8_ALTCOV)<0.99)

Cov = Cov[!hetSNPs]

# remove parent columns if necessary
#Cov = Cov[,1:25]

# save table
write.table(Cov,'aseReadCounts_c172_PallF1_infSites_STAR_duprem.txt')

#merge all coverage files letting NA's fill unmatched positions - use DATA TABLE !!!
Cov = data.frame(pos=uPos)
rownames(Cov) = Cov[,'pos']
Cov = setDT(Cov)
for (i in names(cov)){
	bed = data.table(cov[[i]])
	bed[,'pos'] = paste(cov[[i]][,1],cov[[i]][,2],sep=':')
	colnames(bed)[6] = paste(i,'REFCOV',sep='_') 
	colnames(bed)[7] = paste(i,'ALTCOV',sep='_') 
	colnames(bed)[8] = paste(i,'TOTCOV',sep='_') 
	Cov = merge(Cov,bed[,c(paste(i,'TOTCOV',sep='_'),'pos'),with=F],all=T,by='pos')
}	

Cov = Cov[!hetSNPs]

# remove parent columns if necessary
#Cov = Cov[,1:13]

write.table(Cov,'aseReadCounts_c172_PallF1_infSites_totalCov_STAR_duprem.txt')



