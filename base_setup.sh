# base_setup.sh

# make main directories
# specific to selection analyses (fst, dxy, Tajima's D, RAiSD)
# make directories for intermediate files-- will fail if these don't exist

mkdir -p ${OUTDIR}/analyses/
mkdir -p ${OUTDIR}/datafiles/
mkdir -p ${OUTDIR}/referencelists/


# make reference files

# Generate scaffold list
if [ -f "${OUTDIR}/referencelists/SCAFFOLDS.txt" ];
        then
            echo "SCAFFOLDS.txt already exists, moving on!"
        else
        awk '{print $1}' "${REF}.fai" > "${OUTDIR}/referencelists/SCAFFOLDS.all.txt"
        grep "$CHRLEAD" "${OUTDIR}/referencelists/SCAFFOLDS.all.txt" > "${OUTDIR}/referencelists/SCAFFOLDS.chroms.txt"
        grep -v "$SEXCHR" "${OUTDIR}/referencelists/SCAFFOLDS.chroms.txt" > "${OUTDIR}/referencelists/SCAFFOLDS.txt"
fi


# Make a file with chromosome name and length of chromosome
awk 'BEGIN {OFS = "\t"} {print $1,$2}' ${REF}.fai | grep ${CHRLEAD} | grep -v ${SEXCHR} > ${OUTDIR}/referencelists/autosomes_lengths.txt

while IFS=',' read -r first second; do
    sed -i "s/$second/$first/g" ${OUTDIR}/referencelists/autosomes_lengths.txt 
done <<< "$CHROM"

# Make a comma separated chromosome conversion file without a header where the first column is the name of the chromosome and the second is the name of the associated scaffold in the reference genome:

if [ -f "${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt" ]
        then
            echo "Chromosome conversion table already complete, moving on!"
        else
        echo '1,NC_044571.1' > ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '2,NC_044572.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '3,NC_044573.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '4,NC_044574.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '5,NC_044575.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '6,NC_044576.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '7,NC_044577.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '8,NC_044578.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '9,NC_044579.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '10,NC_044580.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '11,NC_044581.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '12,NC_044582.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '13,NC_044583.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '14,NC_044584.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '15,NC_044585.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '1A,NC_044586.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '17,NC_044587.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '18,NC_044588.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '19,NC_044589.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '20,NC_044590.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '21,NC_044591.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '22,NC_044592.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '23,NC_044593.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '24,NC_044594.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '25,NC_044595.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '26,NC_044596.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '27,NC_044597.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '28,NC_044598.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '29,NC_044599.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '4A,NC_044600.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo 'Z,NC_044601.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
fi


