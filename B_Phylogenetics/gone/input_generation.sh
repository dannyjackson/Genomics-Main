#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file>

This uses VCFtools and PLINK to generate input files for GONE2
I recommend running it as a slurm array to pass individuals to sbatch jobs for maximum efficiency

Required argument:
  -p  Path to the GONE parameter file (e.g., params_gone2.sh in the GitHub repository).
  -s  Population specific-parameter file.
  -r  Run name, required for providing a unique name to output files (especially when using in an array)."
    exit 1
fi

# Parse command-line arguments
while getopts p:s:r: option; do
    case "${option}" in
        p) PARAMS=${OPTARG};;
        s) POPPARAMS=${OPTARG};;
        r) RUNNAME=${OPTARG};;
        *) echo "Invalid option: -${OPTARG}" >&2; exit 1;;
    esac
done

if [ -z "${PARAMS}" ]; then
    echo "Error: No parameter file provided." >&2
    exit 1
fi

# Load parameters
source "${PARAMS}"
source "${POPPARAMS}"

printf "\n\n\n\n"
date
echo "Current script: input_generation.sh"


echo ARRAY NAME: $RUNNAME

echo VCF FILE PATH: $VCFFILE
echo SCAFFOLD SUBSET: $SCAFFOLD_LIST
echo POPULATION: $POPNAME
echo MAX MISSING DATA: $MAX_MISSING
echo MINOR ALLELE FREQ: $MAF

CHR_SUBSET_FLAGS=$(for name in $(cat $SCAFFOLD_LIST); do echo --chr $name; done)

RESULT_DIR=${OUTDIR}/datafiles/gone2_inputs/${RUNNAME}/${POPNAME}

if [ ! -d "$RESULT_DIR" ]; then
  echo "Directory for gone2 input for ${RUNNAME} with ${POPNAME} does not exist. Creating it now..."
  mkdir -p "$RESULT_DIR" # -p creates parent directories if they don't exist
else
  echo "Directory for gone2 input for ${RUNNAME} with ${POPNAME} already exists. WARNING: Existing files in this directory may be overwritten."
fi

# Create a filtered VCF to only include a specific subset of chromosomes
# --max-missing to set proportion of missing data you'll permit.
vcftools $CHR_SUBSET_FLAGS --gzvcf $VCFFILE --recode --recode-INFO-all --max-missing $MAX_MISSING --out ${RESULT_DIR}/${POPNAME}

# Convert VCF to plink formats
# --allow-extra-chr to deal with non-standard chromosome names
# --thin-count to randomly sample a subset of SNPs from the VCF (needed since GONE2 can't handle too many SNPs by default)
plink --vcf ${RESULT_DIR}/${POPNAME}.recode.vcf --maf $MAF --allow-extra-chr --recode --out ${RESULT_DIR}/${POPNAME}
# Note that if you don't have info on SNP position in a genetic map (cM), you'll probably need to set a fixed recombination rate when running GONE2 or do some extra analyses to find this info yourself. (In this case the .map file is not useful)

# Due to not having standard chromosome names and that outputted files aren't always consistent with what GONE2 wants, we use these python scripts to reformat the data. (No filtering or analyses done here, just moving the numbers around)
python map_clean.py --map ${RESULT_DIR}/${POPNAME}.map
python ped_clean.py --ped ${RESULT_DIR}/${POPNAME}.ped