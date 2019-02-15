

devtools::load_all('../gws')
devtools::load_all('~/moiR')

#### Read the sweep output ####

load('../gws/dataint/gsweep.rda')
head(s)
tail(s)

s$chr = fn(s$chr)
s$pos = fn(s$pos)


load('../gws/dataint/gmetrics.rda')
gmetrics$chr = fn(gmetrics$chr)
gmetrics$pos = fn(gmetrics$pos)


gwares<-readRDS('../gws/dataint/rFitness.rda')
head(gwares)

####  plot for the genes
#
# ACD6: AT4G14400
# Chr4:8292201-8299001
#
# MELA1: AT1G62045 (the QTL on chromosome 1)
# Chr1: 22934601-22936101
#
# MELA2: new annotation (homolog of MELA1)
# Chr1: 3964801-3963001
#

acd6plus<-function(){

  d<-rbind(
    c('acd6', 4,8292201,8299001),
    c('mela1', 1,22934601,22936101),
    c('mela2', 4,3964801,3963001)
    )
  colnames(d)<- c('gene','chr','pos1','pos2')
  rownames(d)<- c('acd6','mela1','mela2')

  d<-data.frame(d)
  d$chr <- moiR::fn(d$chr)
  d$pos1 <- moiR::fn(d$pos1)
  d$pos2 <- moiR::fn(d$pos2)


  return(d)
}
acd6plus()
genes<-acd6plus()


### Plot sweeps


pdf('Sweep_screening_for_Wangsheng.pdf', height = 5, width = 12)
par(mfrow(c(3,1)))

mygene<-'acd6'

tmp<- dplyr::filter(s,chr== genes[mygene,'chr'] & location > genes[mygene,'pos1'] - 5e4 & location < genes[mygene,'pos2'] + 5e4)
plot(tmp$LR ~ tmp$location,type='line',xlab='position (bp)', ylab = 'Sweep likelihood ratio', main='Acd6')
rect(xleft = genes[mygene,'pos1'] , xright =  genes[mygene,'pos2'],ybottom = -10,ytop = 10 ,col=transparent('black',0.1))

mygene<-'mela1'

tmp<- dplyr::filter(s,chr== genes[mygene,'chr'] & location > genes[mygene,'pos1'] - 5e4 & location < genes[mygene,'pos2'] + 5e4)
plot(tmp$LR ~ tmp$location,type='line',xlab='position (bp)', ylab = 'Sweep likelihood ratio', main='Mela1')
rect(xleft = genes[mygene,'pos1'] , xright =  genes[mygene,'pos2'],ybottom = -10,ytop = 10 ,col=transparent('black',0.1))

mygene<-'mela2'

tmp<- dplyr::filter(s,chr== genes[mygene,'chr'] & location > genes[mygene,'pos1'] - 5e4 & location < genes[mygene,'pos2'] + 5e4)
plot(tmp$LR ~ tmp$location,type='line',xlab='position (bp)', ylab = 'Sweep likelihood ratio', main='Mela2')
rect(xleft = genes[mygene,'pos1'] , xright =  genes[mygene,'pos2'],ybottom = -10,ytop = 10 ,col=transparent('black',0.1))

dev.off()


### Plot Fst
head(gmetrics)

pdf('Fst_screening_for_Wangsheng.pdf', height = 5, width = 12)
par(mfrow(c(3,1)))

myvar<- 'fst'
# myvar<- 'taj'
mygene<-'acd6'

tmp<- dplyr::filter(gmetrics,chr== genes[mygene,'chr'] & pos > genes[mygene,'pos1'] - 5e4 & pos < genes[mygene,'pos2'] + 5e4)
plot(tmp[,myvar] ~ tmp$pos,type='points',xlab='position (bp)', ylab = 'Fst', main='Acd6')
rect(xleft = genes[mygene,'pos1'] , xright =  genes[mygene,'pos2'],ybottom = -10,ytop = 10 ,col=transparent('black',0.1))

mygene<-'mela1'

tmp<- dplyr::filter(gmetrics,chr== genes[mygene,'chr'] & pos > genes[mygene,'pos1'] - 5e4 & pos < genes[mygene,'pos2'] + 5e4)
plot(tmp[,myvar] ~ tmp$pos,type='points',xlab='position (bp)', ylab = 'Fst', main='Mela1')
rect(xleft = genes[mygene,'pos1'] , xright =  genes[mygene,'pos2'],ybottom = -10,ytop = 10 ,col=transparent('black',0.1))

mygene<-'mela2'

tmp<- dplyr::filter(gmetrics,chr== genes[mygene,'chr'] & pos > genes[mygene,'pos1'] - 5e4 & pos < genes[mygene,'pos2'] + 5e4)
plot(tmp[,myvar] ~ tmp$pos,type='points',xlab='position (bp)', ylab = 'Fst', main='Mela2')
rect(xleft = genes[mygene,'pos1'] , xright =  genes[mygene,'pos2'],ybottom = -10,ytop = 10 ,col=transparent('black',0.1))

dev.off()


### Plot gwa
# gwares<-readRDS('../gws/dataint/rSeeds_thi.rda')

gwares<-readRDS('../gws/dataint/rFitness.rda')
tops<-readRDS('../gws/dataint/rFitnesstopSNPs.rda')



pdf('Relativefitnes_screening_for_Wangsheng.pdf', height = 5, width = 12)
par(mfrow(c(3,1)))


myvar<- 'beta'

mygene<-'acd6'

tmp<- dplyr::filter(gwares,chr== genes[mygene,'chr'] & pos > genes[mygene,'pos1'] - 5e4 & pos < genes[mygene,'pos2'] + 5e4)
plot(tmp[,myvar] ~ tmp$pos,type='points',xlab='position (bp)', ylab = 'Relative fitness effect in Madrid & drought', main='Acd6')
rect(xleft = genes[mygene,'pos1'] , xright =  genes[mygene,'pos2'],ybottom = -10,ytop = 10 ,col=transparent('black',0.1))

mygene<-'mela1'

tmp<- dplyr::filter(gwares,chr== genes[mygene,'chr'] & pos > genes[mygene,'pos1'] - 5e4 & pos < genes[mygene,'pos2'] + 5e4)
plot(tmp[,myvar] ~ tmp$pos,type='points',xlab='position (bp)', ylab = 'Relative fitness effect in Madrid & drought', main='Mela1')
rect(xleft = genes[mygene,'pos1'] , xright =  genes[mygene,'pos2'],ybottom = -10,ytop = 10 ,col=transparent('black',0.1))

mygene<-'mela2'

tmp<- dplyr::filter(gwares,chr== genes[mygene,'chr'] & pos > genes[mygene,'pos1'] - 5e4 & pos < genes[mygene,'pos2'] + 5e4)
plot(tmp[,myvar] ~ tmp$pos,type='points',xlab='position (bp)', ylab = 'Relative fitness effect in Madrid & drought', main='Mela2')
rect(xleft = genes[mygene,'pos1'] , xright =  genes[mygene,'pos2'],ybottom = -10,ytop = 10 ,col=transparent('black',0.1))



myvar<- 'p_lrt'

mygene<-'acd6'

tmp<- dplyr::filter(gwares,chr== genes[mygene,'chr'] & pos > genes[mygene,'pos1'] - 5e4 & pos < genes[mygene,'pos2'] + 5e4)
plot(-log10(tmp[,myvar]) ~ tmp$pos,type='points',xlab='position (bp)', ylab = 'Significance (-log10 p value) of effect in Madrid & drought', main='Acd6')
rect(xleft = genes[mygene,'pos1'] , xright =  genes[mygene,'pos2'],ybottom = -10,ytop = 10 ,col=transparent('black',0.1))

mygene<-'mela1'

tmp<- dplyr::filter(gwares,chr== genes[mygene,'chr'] & pos > genes[mygene,'pos1'] - 5e4 & pos < genes[mygene,'pos2'] + 5e4)
plot(-log10(tmp[,myvar])~ tmp$pos,type='points',xlab='position (bp)', ylab = 'Significance (-log10 p value) of effect in Madrid & drought', main='Mela1')
rect(xleft = genes[mygene,'pos1'] , xright =  genes[mygene,'pos2'],ybottom = -10,ytop = 10 ,col=transparent('black',0.1))

mygene<-'mela2'

tmp<- dplyr::filter(gwares,chr== genes[mygene,'chr'] & pos > genes[mygene,'pos1'] - 5e4 & pos < genes[mygene,'pos2'] + 5e4)
plot(-log10(tmp[,myvar])~ tmp$pos,type='points',xlab='position (bp)', ylab = 'Significance (-log10 p value) of effect in Madrid & drought', main='Mela2')
rect(xleft = genes[mygene,'pos1'] , xright =  genes[mygene,'pos2'],ybottom = -10,ytop = 10 ,col=transparent('black',0.1))

dev.off()

## At the accession level

necrosis<-read.table('../acd6/accessionsnecrosis.tsv',fill=T,header=T)
necrosis
head(necrosis)

d<-merge(necrosis, by.x='Accession', dry, by.y='id')
head(d)

require(ggplot2)
require(latex2exp)
require(cowplot)

ggdotscolor(d$Necrosis_severity, d$Fitness_thi,
            xlab='Necrosis severity', ylab='Fitness in field Madrid (# seeds)') +
            geom_hline(yintercept=0, lty='dotted')+
            xlim(c(0,6))+
            stat_smooth(method='glm') 


dplyr::select(d, contains('Fitness') , Necrosis_severity) %>%
  moiR::cormatrix_r_p()  %>%
  write.tsv(.,file='Correlation_fitness_necrossis.tsv')

dplyr::select(d, contains('Fitness') , Necrosis_severity) %>%
  filter(Necrosis_severity < 6) %>%
  moiR::cormatrix_r_p( )%>%
    write.tsv(.,file='Correlation_fitness_necrossis_lowerthan6.tsv')

dplyr::select(d, contains('Survival_fruit') , Necrosis_severity) %>%
  # filter(Necrosis_severity < 6) %>%
  moiR::cormatrix_r_p() %>%
  write.tsv(.,file='Correlation_survival_necrossis.tsv')

dplyr::select(d, contains('Seeds') , Necrosis_severity) %>%
  # filter(Necrosis_severity < 6) %>%
  moiR::cormatrix_r_p() 
  write.tsv(.,file='Correlation_seeds_necrossis.tsv')

  
#### drought phenotype ####

d2<-merge(necrosis, by.x='Accession', drought, by.y='id')
head(d2)
dim(d2)  

plot(d2$polqua, d2$Necrosis_severity)
cor.test(d2$polqua, d2$Necrosis_severity)
cor.test(d2$polqua, d2$Necrosis_severity)

d2 %>%
  dplyr::select(.,contains(day), Necrosis_severity) %>%
  moiR::cormatrix_r_p() 
head(d2)
  