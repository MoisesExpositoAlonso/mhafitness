####************************************************************************####
devtools::load_all('.')
library(dplyr)

####************************************************************************####

####  load genomes
genomes<-readRDS('../gws/databig/genomenofilters.rda')
map<-genomes$map
head(map)
  map$SNP<-paste(map$chromosome, map$physical.pos,sep='_')
fam<-genomes$fam
head(fam)

####  whangsheng positions
# MHA1 Chr1: 22,935,037, p=10-11.9), and MHA2 (Chr4: 11,019,243
# mha<-c('1_22935037','4_11019243','4_8298247')
# mha<-c('1_22935037','4_11019243','4_8298244')
mha<-c('1_22935037','4_11019243','4_8298244')

map$SNP[grep('4_829824' , map$SNP)]


####  find what position is the marker and extract genome
# map$marker.ID [grep('1_2293', map$marker.ID) ]
source('R/readgenome.R')
# sg<-read_n_subset_genome(selectedsnps = mha,genomefile='databig/genomenofilters.rda')
# sg<-read_n_subset_genome(selectedsnps = mha,genomefile='databig/1001gbi.rda')
# sg<-read_n_subset_genome(selectedsnps = mha,genomefile='databig/1001gbi.rda')

genomefile='databig/mha.rda'
genomes<-readRDS(genomefile)
map<-genomes$map
fam<-genomes$fam
# map$SNP<-paste0(map$chromosome, "_",map$physical.pos)
Go<-bigmemory::attach.big.matrix(paste0(tools::file_path_sans_ext(genomefile),'.desc'))

sg<- data.frame(Go[,])
sg<-sg / 2
head(sg)
# colnames(sg)<-paste0('m',colnames(sg) )
colnames(sg)<-paste0('m',mha )
head(sg)

#####**********************************************************************#####
####  necrosis
dat<-read.delim('../acd6/data-raw/accessionsnecrosis.tsv')
dm<-merge(fam, by.x='sample.ID', dat, by.y='Accession',all.x=T)
dm<-cbind(dm,sg)

dim(dm)

table(dat$Accession %in%fam$sample.ID)


# lm(data = dm, Necrosis_severity ~ m1_22935037 + m4_11019243 ) %>%  summary
# lm(data = dm, Necrosis_severity ~ m1_22935037 * m4_8298244 ) %>%  summary
# lm(data = dm, Necrosis_severity ~ m4_11019243 * m4_8298244 ) %>%  summary
# lm(data = dm, Necrosis_severity ~ m1_22935037  ) %>%  summary
# lm(data = dm, Necrosis_severity ~ m4_11019243 ) %>%  summary
# lm(data = dm, Necrosis_severity ~ m4_8298244 ) %>%  summary
# 
# lm(data = dm, Necrosis_severity ~ m1_22935037 * m4_11019243  ) %>%  summary
# 
# dmnona<-na.omit(dm)
# randomForest(data = dmnona, Necrosis_severity ~ m1_22935037 * m4_11019243)

#####**********************************************************************#####
####  fitness associations
load('../gws/data/field.s.rda')
fit<-field.s

head(sg)
head(fam)
head(sg)

# merge
genomemat<-data.frame(fam$sample.ID,sg)
head(genomemat)
fitg<-merge(fit,by.x='id',genomemat, by.y='fam.sample.ID')
head(fitg)


# lm(data = fieldfilter(fitg,'mli'), Seeds ~ m1_22935037 * m4_11019243 *m4_8298244 ) %>%  summary
# lm(data = fieldfilter(fitg,'mlp'), Seeds ~ m1_22935037 * m4_11019243 *m4_8298244 ) %>%  summary
# lm(data = fieldfilter(fitg,'mhi'), Seeds ~ m1_22935037 * m4_11019243 *m4_8298244 ) %>%  summary
# lm(data = fieldfilter(fitg,'mhp'), Seeds ~ m1_22935037 * m4_11019243 *m4_8298244 ) %>%  summary
# lm(data = fieldfilter(fitg,'thi'), Seeds ~ m1_22935037 * m4_11019243 *m4_8298244 ) %>%  summary
# lm(data = fieldfilter(fitg,'thp'), Seeds ~ m1_22935037 * m4_11019243 *m4_8298244 ) %>%  summary
# lm(data = fieldfilter(fitg,'tlp'), Seeds ~ m1_22935037 * m4_11019243 *m4_8298244 ) %>%  summary
# lm(data = fieldfilter(fitg,'tli'), Seeds ~ m1_22935037 * m4_11019243 *m4_8298244 ) %>%  summary

fitg$env<-paste(fitg$site,fitg$water,fitg$indpop,sep="_")
fitg$m1_22935037[is.na(fitg$m1_22935037)]<-0.5
fitg$m4_11019243[is.na(fitg$m4_11019243)]<-0.5
fitg$m4_8298244[is.na(fitg$m4_8298244)]<-0.5

library(MCMCglmm)
m<-MCMCglmm(data = fitg, 
         fixed = Seeds ~ 
           m1_22935037:m4_11019243 +
           m1_22935037:m4_8298244 +
           m4_11019243:m4_8298244,
         random = ~env
         ) 
m %>%  summary
