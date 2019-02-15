
outplink='mha'

plink  --bfile "../plink/1001gbi" --extract mhagenes.txt --recode --make-bed  --out $outplink # the nofilters vcf needs to be biallelic for downstream analyses

plink --bfile $outplink --recode --tab --out $outplink
prettymap $outplink".map"
prettybim $outplink".bim"

plink --bfile $outplink  --recodeA --out $outplink # only recodeA to skip heterosygote call

cat $outplink".raw" | cut -d " "  -f 7-  > $outplink".012"


