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
