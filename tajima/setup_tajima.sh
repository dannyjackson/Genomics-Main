# make directories for intermediate files-- will fail if these don't exist

mkdir -p ${OUTDIR}/analyses/{tajima,genelist,thetas}
mkdir -p ${OUTDIR}/datafiles/{safs,mls}
mkdir -p ${OUTDIR}/analyses/genelist/${POP}/
mkdir -p ${OUTDIR}/analyses/genelist/${POP}/${WIN}
mkdir -p ${OUTDIR}/analyses/tajima/${POP}
mkdir -p ${OUTDIR}/analyses/tajima/${POP}/${WIN}
mkdir -p ${OUTDIR}/analyses/tajima/${POP}/snps
