# make directories for intermediate files-- will fail if these don't exist

mkdir -p ${OUTDIR}/analyses/{genelist,thetas}
mkdir -p ${OUTDIR}/datafiles/{safs,mls}
mkdir -p ${OUTDIR}/analyses/genelist/${POP}/
mkdir -p ${OUTDIR}/analyses/genelist/${POP}/${WIN}
mkdir -p ${OUTDIR}/analyses/thetas/${POP}
mkdir -p ${OUTDIR}/analyses/thetas/${POP}/${WIN}
mkdir -p ${OUTDIR}/analyses/thetas/${POP}/snps
mkdir -p ${OUTDIR}/analyses/Tajima/
mkdir -p ${OUTDIR}/analyses/Tajima/${POP}
mkdir -p ${OUTDIR}/analyses/Tajima/${POP}/${WIN}

