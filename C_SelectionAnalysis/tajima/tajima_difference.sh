# Get Tajima's D Difference between populations

source params_base.sh

TAJIMADIR=${OUTDIR}/analyses/tajima
POP1=pop1
POP2=pop1
WIN=50000
OUTPREFIX=some_name

(echo -e "chromo\tposition\tTajima"; awk 'BEGIN {OFS="\t"} NR==FNR && FNR>1 {data[$2,$3] = $9; next} FNR>1 && ($2,$3) in data {print $2, $3, $9 - data[$2,$3]}' ${TAJIMADIR}/${POP1}/${POP1}.Tajima.${WIN}.Ztransformed.csv ${TAJIMADIR}/${POP2}/${POP2}.Tajima.${WIN}.Ztransformed.csv) \
| grep -v 'NW' | grep -v 'Z' >> ${TAJIMADIR}/${OUTPREFIX}.txt