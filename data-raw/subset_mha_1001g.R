library(devtools)
library(bigmemory)
library(Rcpp)
load_all('../gws')

################################################################################
### Read the plink files
################################################################################

# system("rm databig/g.*")
# g012file<-'data-raw/515g.012'
# myname<-"genome"
#
# system(paste0('rm databig/',myname,'*') )
#
# # Reading the bedfile and storing the data in temporary directory
# genomes <- readplink(g012file, backingfile = myname,backingpath = 'databig/')
#
# map<-genomes$map[,-c(2:3)]
# colnames(map)[1:2]<-c('chr','pos')
# map$SNP<-paste(map$chr, map$pos,sep="_")
#
# fam<-genomes$fam
# fam$id<-fam$sample.ID

####  save
# devtools::use_data(genomes,overwrite = TRUE)
# devtools::use_data(map,overwrite = F)
# devtools::use_data(fam,overwrite = F)


#####**********************************************************************#####
#### all biallelic positions of 1001 G ####
g012file<-'mha.012'
myname<-"mha"

# Reading the bedfile and storing the data in temporary directory
genomes <- readplink(g012file, backingfile = myname,backingpath = 'databig/')
saveRDS(file=paste0("databig/",myname,".rda" ), genomes)
