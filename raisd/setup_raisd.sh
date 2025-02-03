# make directories for intermediate files-- will fail if these don't exist

mkdir -p ${OUTDIR}/analyses/{raisd,genelist}
mkdir -p ${OUTDIR}/datafiles/{safs,mls}
mkdir -p ${OUTDIR}/analyses/genelist/${POP1}_${POP2}/
mkdir -p ${OUTDIR}/analyses/genelist/${POP1}_${POP2}/${WIN}
mkdir -p ${OUTDIR}/analyses/raisd/${POP1}_${POP2}
mkdir -p ${OUTDIR}/analyses/raisd/${POP1}_${POP2}/${WIN}
mkdir -p ${OUTDIR}/analyses/raisd/${POP1}_${POP2}/snps