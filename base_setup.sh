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
        grep -v "$CHRLEAD" "${OUTDIR}/referencelists/SCAFFOLDS.all.txt" > "${OUTDIR}/referencelists/SCAFFOLDS.unplaced.txt"
        grep -v "$SEXCHR" "${OUTDIR}/referencelists/SCAFFOLDS.chroms.txt" > "${OUTDIR}/referencelists/SCAFFOLDS.txt"
fi

# Make a comma separated chromosome conversion file without a header where the first column is the name of the chromosome and the second is the name of the associated scaffold in the reference genome:

if [ -f "${CHR_FILE}" ]
        then
            echo "Chromosome conversion table already complete, moving on!"
        else
            echo "Creating chromosome conversion table from NCBI RefGenome..."

            if micromamba env list | grep -Eq "\\bncbi_datasets\\b"
                then 
                    echo "ncbi_datasets environment is found"
                else
                    echo "ncbi_datasets environment is not found, creating it now..."
                    micromamba create -n ncbi_datasets -c conda-forge ncbi-datasets-cli
                fi

            micromamba activate ncbi_datasets  # Acivate environment with NCBI toolkit installed
            # Make chrom conversion file
            datasets summary genome accession ${REF_ACC} --report sequence --as-json-lines | dataformat tsv genome-seq --fields chr-name,refseq-seq-acc > ${OUTDIR}/referencelists/chrom_name_mapping.txt #Unfortunately I don't think ncbi_datasets toolkit has the ability to output as csv, so we'll have to do that ourselves in the python script below.
            python ${SCRIPTDIR}/Genomics-Main/general_scripts/make_chrom_conversion_file.py -i ${OUTDIR}/referencelists/chrom_name_mapping.txt -o ${CHR_FILE} -e ${SEXCHR},${SCAF_LEAD}
            
            # Make scaffold lengths file
            datasets summary genome accession ${REF_ACC} --report sequence --as-json-lines | dataformat tsv genome-seq --fields refseq-seq-acc,seq-length > ${OUTDIR}/referencelists/scaffold_lengths.txt # Get scaffold lengths
            grep ${CHRLEAD} scaffold_lengths.txt | grep -v ${SEXCHR} > autosomes_lengths.txt # Get just a file of autosome lengths
            
            micromamba deactivate # Deactivate NCBI datasets environment
fi


