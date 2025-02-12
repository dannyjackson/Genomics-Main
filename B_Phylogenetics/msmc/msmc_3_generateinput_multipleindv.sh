4.4.5_msmcinput_multiplegenomes.sh



# MSMC on Multiple Genomes 
# step 2 -- generate input

grep 'BFBO' /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_samplecodes.mindepth3x.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/BFBO_samplecodes.txt
grep 'PEBO' /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_samplecodes.mindepth3x.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/PEBO_samplecodes.txt
grep 'MABO' /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_samplecodes.mindepth3x.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/MABO_samplecodes.txt
grep 'NABO' /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_samplecodes.mindepth3x.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/NABO_samplecodes.txt
grep 'BRBO' /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_samplecodes.mindepth3x.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/BRBO_samplecodes.txt
grep 'RFBO' /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_samplecodes.mindepth3x.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/RFBO_samplecodes.txt

cp /xdisk/mcnew/dannyjackson/sulidae/referencelists/BFBO_samplecodes.txt BFBO_IND.txt
cp /xdisk/mcnew/dannyjackson/sulidae/referencelists/PEBO_samplecodes.txt PEBO_IND.txt
cp /xdisk/mcnew/dannyjackson/sulidae/referencelists/MABO_samplecodes.txt MABO_IND.txt
cp /xdisk/mcnew/dannyjackson/sulidae/referencelists/NABO_samplecodes.txt NABO_IND.txt
cp /xdisk/mcnew/dannyjackson/sulidae/referencelists/BRBO_samplecodes.txt BRBO_IND.txt
cp /xdisk/mcnew/dannyjackson/sulidae/referencelists/RFBO_samplecodes.txt RFBO_IND.txt

echo -e 'BFBO\nPEBO\nMABO\nNABO\nBRBO\nRFBO' > /xdisk/mcnew/dannyjackson/sulidae/referencelists/speciescodes.txt


# need to gzip all vcfs

ls vcf/ > vcf_filenames.txt


#!/bin/bash
module load parallel
parallel -j 12 -a vcf_filenames.txt 'C={}; bgzip -c vcf/$file > gzvcf/$file.gz'

sbatch --account=mcnew \
--job-name=bgzipvcfs \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/bgzipvcfs.%j \
--nodes=1 \
--ntasks-per-node=12 \
--time=6:00:00 \
./bgzipvcfs.sh

# Submitted batch job 3691397

# To combine multiple individuals, the steps generally follow what is detailed above. However, in step 2, you will instead use the script msmc_2_generateInput_multiInd.sh, which will combine all of the individuals from a given population. For this script, which is submitted using submit_2_multi.txt, you will need a file with the names of the individuals you want combined in your MSMC run (named POP_IND.txt, where "POP" is your population name).

# remove sex chromosomes from scaffolds file
# sex_chr=NC_087547.1,NC_087548.1 # name of sex chromosome to omit in analyses
sed -i 's/NC_087547.1//g' ${OUTDIR}/referencelists/SCAFFOLDS.txt
sed -i 's/NC_087548.1//g' ${OUTDIR}/referencelists/SCAFFOLDS.txt


for POP in `cat /xdisk/mcnew/dannyjackson/sulidae/referencelists/speciescodes.txt`;
	do echo $POP
	sbatch --account=mcnew \
	--job-name=msmc_run.${POP} \
        --partition=standard \
	--mail-type=ALL \
	--output=slurm_output/msmc_run.${POP}.%j \
	--nodes=1 \
	--ntasks-per-node=1 \
	--time=10:00:00 \
    ~/programs/msmc2_scripts/msmc_2_generateInput_multiInd.sh `echo ${POP}_IND.txt` `echo $POP`
done

BFBO
Submitted batch job 3692205
PEBO
Submitted batch job 3692206
MABO
Submitted batch job 3692207
NABO
Submitted batch job 3692208
BRBO
Submitted batch job 3692209
RFBO
Submitted batch job 3692210

sbatch --account=mcnew \
	--job-name=msmc_run.${POP} \
    --partition=standard \
	--mail-type=ALL \
	--output=slurm_output/msmc_run.${POP}.%j \
	--nodes=1 \
	--ntasks-per-node=1 \
	--time=10:00:00 \
    ~/programs/msmc2_scripts/msmc_2_generateInput_multiInd.sh BFBO_IND.txt BFBO









# test a bunch of options with different -p parameters
# parone
# P_PAR=1*2+25*1+1*2+1*3
# P_PAR=1*2+15*1+1*2
# store output in par1 or par2 directory


for POP in `cat /xdisk/mcnew/dannyjackson/sulidae/referencelists/speciescodes.txt`;
	do echo $POP
	sbatch --account=mcnew \
	--job-name=msmc_run.${POP} \
    --partition=standard \
	--mail-type=ALL \
	--output=slurm_output/msmc_run.${POP}.%j \
	--nodes=1 \
	--ntasks-per-node=1 \
	--time=10:00:00 \
    ~/programs/msmc2_scripts/msmc_3_runMSMC.sh `echo $POP` par_2
done


# run this as soon as 3692253 is finished
sbatch --account=mcnew \
--job-name=msmc_run.${POP} \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/msmc_run.${POP}.%j \
--nodes=1 \
--ntasks-per-node=1 \
--time=10:00:00 \
~/programs/msmc2_scripts/msmc_3_runMSMC.sh BFBO

# Submitted batch job 3694138

BFBO
Submitted batch job 3692246 # cancelled
PEBO
Submitted batch job 3692247
MABO
Submitted batch job 3692248
NABO
Submitted batch job 3692249
BRBO
Submitted batch job 3692250
RFBO
Submitted batch job 3692251




# In step 3, you will then need to change -I, which specifies which haplotypes to analyze. For one individual, this is 0,1, for two individuals 0,1,2,3, and so on. This will be done automatically in the “multiInd" script, which can simply be run as follows:

# for multiple individuals in a population
sbatch msmc_3_runMSMC_slurm.sh
# or without using slurm





# now make plots in R

for i in `cat /xdisk/mcnew/dannyjackson/sulidae/referencelists/speciescodes.txt`;
	do echo $i
    d=`date +%m%d%y`
	Rscript msmc_plot.multiind.r /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc $i $d
done


scp -r dannyjackson@filexfer.hpc.arizona.edu:/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/plots /Users/danjack/Documents/SulidaePaper/Summer2024/RawFigures


/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/output/par_1/msmc_output.msmc_BFBO_012825.final.txt
/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/output/par_1/msmc_output.msmc_PEBO_012825.final.txt
/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/output/par_1/msmc_output.msmc_NABO_012825.final.txt
/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/output/par_1/msmc_output.msmc_MABO_012825.final.txt
/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/output/par_1/msmc_output.msmc_RFBO_012825.final.txt
/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/output/par_1/msmc_output.msmc_BRBO_012825.final.txt

mu <- 1.25e-8
gen <- 10 # average GenLength from Supplementary Table 4 in Generation lengths of the world's birds and their implications for extinction risk https://conbio.onlinelibrary.wiley.com/doi/epdf/10.1111/cobi.13486

# bfboDat<-read.table("/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/output/par_1/msmc_output.msmc_BFBO_012825.final.txt", header=TRUE)
peboDat<-read.table("/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/output/par_1/msmc_output.msmc_PEBO_012825.final.txt", header=TRUE)

maboDat<-read.table("/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/output/par_1/msmc_output.msmc_MABO_012825.final.txt", header=TRUE)
naboDat<-read.table("/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/output/par_1/msmc_output.msmc_NABO_012825.final.txt", header=TRUE)

brboDat<-read.table("/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/output/par_1/msmc_output.msmc_BRBO_012825.final.txt", header=TRUE)
rfboDat<-read.table("/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/output/par_1/msmc_output.msmc_RFBO_012825.final.txt", header=TRUE)


png(paste0("/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc//plots/all.msmc.012825.png"), width=500, height=500)
plot(maboDat$left_time_boundary/mu*gen, (1/maboDat$lambda)/(2*mu), log="x",ylim=c(0,100000),
     type="n", xlab="Years ago", ylab="effective population size")
# lines(bfboDat$left_time_boundary/mu*gen, (1/bfboDat$lambda)/(2*mu), type="s", col="green")
lines(peboDat$left_time_boundary/mu*gen, (1/peboDat$lambda)/(2*mu), type="s", col="blue")
lines(maboDat$left_time_boundary/mu*gen, (1/maboDat$lambda)/(2*mu), type="s", col="yellow")
lines(naboDat$left_time_boundary/mu*gen, (1/naboDat$lambda)/(2*mu), type="s", col="orange")
lines(brboDat$left_time_boundary/mu*gen, (1/brboDat$lambda)/(2*mu), type="s", col="brown")
lines(rfboDat$left_time_boundary/mu*gen, (1/rfboDat$lambda)/(2*mu), type="s", col="red")
legend("top",legend=c("Peruvian", "Masked", "Nazca", "Brown", "Red-footed"), col=c("blue", "yellow", "orange", "brown", "red"), lty=c(1,1))
dev.off()




png(paste0("/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc//plots/all.msmc.012825.full.png"), width=500, height=500)
plot(peboDat$left_time_boundary/mu*gen, (1/peboDat$lambda)/(2*mu), log="x",
     type="n", xlab="Years ago", ylab="effective population size")
# lines(bfboDat$left_time_boundary/mu*gen, (1/bfboDat$lambda)/(2*mu), type="s", col="green")
lines(peboDat$left_time_boundary/mu*gen, (1/peboDat$lambda)/(2*mu), type="s", col="blue")
lines(maboDat$left_time_boundary/mu*gen, (1/maboDat$lambda)/(2*mu), type="s", col="yellow")
lines(naboDat$left_time_boundary/mu*gen, (1/naboDat$lambda)/(2*mu), type="s", col="orange")
lines(brboDat$left_time_boundary/mu*gen, (1/brboDat$lambda)/(2*mu), type="s", col="brown")
lines(rfboDat$left_time_boundary/mu*gen, (1/rfboDat$lambda)/(2*mu), type="s", col="red")
legend("top",legend=c("Peruvian", "Masked", "Nazca", "Brown", "Red-footed"), col=c("blue", "yellow", "orange", "brown", "red"), lty=c(1,1))
dev.off()



# generate combined input 
cat BFBO_IND.txt > BFBO_PEBO.txt
cat PEBO_IND.txt >> BFBO_PEBO.txt

cat MABO_IND.txt > MABO_NABO.txt
cat NABO_IND.txt >> MABO_NABO.txt

sbatch --account=mcnew \
--job-name=msmc_run.MABO_NABO \
    --partition=standard \
--mail-type=ALL \
--output=slurm_output/msmc_run.MABO_NABO.%j \
--nodes=1 \
--ntasks-per-node=1 \
--time=5:00:00 \
~/programs/msmc2_scripts/msmc_2_generateInput_multiInd.sh `echo MABO_NABO.txt` MABO_NABO

# Submitted batch job 3694141
# Run MABO NABO msmc
less ~/programs/msmc2_scripts/msmc_3_runMSMC.sh

for i in `cat MABO_NABO.txt`
do echo $i
IND=$i
    for s in `cat ${OUTDIR}/referencelists/SCAFFOLDS.txt`
        do echo $s
        ls /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/input/msmc_input.${IND}.${s}.txt >> /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/input/SCAFS_INPUT_MABO_NABO
    done
done



#!/bin/sh

MSMC_INPUT=`cat /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/input/SCAFS_INPUT_MABO_NABO`
~/programs/msmc_2.0.0_linux64bit  -t 4 -p 1*2+15*1+1*2 -i 100 -I 0,0,0,0,1,1,1,1,1 -o /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/output/MABO_NABO $MSMC_INPUT

/programs/msmc_2.0.0_linux64bit -t 16 -p 1*2+15*1+1*2 -i 100 -o $MSMC_OUTPUT -I 0,1 $MSMC_INPUT

sbatch --account=mcnew \
--job-name=msmc_run.MABONABO \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/msmc_run.MABONABO.%j \
--nodes=1 \
--ntasks-per-node=8 \
--time=3:00:00 \
msmc_mabonabo.sh

Submitted batch job 12068012

# estimate cross coalescence
# ~/programs/msmc_2.0.0_linux64bit -t 16 -p $P_PAR -i 100 -o $MSMC_OUTPUT -I 0,1 $MSMC_INPUT

msmc2 -I 0,1,2,3 -o within_pop1 `cat INPUT_LIST.txt`
msmc2 -I 4,5,6,7 -o within_pop2 `cat INPUT_LIST.txt`
msmc2 -P 0,0,0,0,1,1,1,1 -o between_pop1-2 `cat INPUT_LIST.txt`

python combineCrossCoal.py \
    /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/output/par_2/between_pop1-2.final.txt \
    /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/output/par_2/msmc_output.msmc_MABO_012825.final.txt \
    /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/output/par_2/msmc_output.msmc_NABO_012825.final.txt > combined_pop1-2.final.txt

~/programs/msmc-tools/combineCrossCoal.py MABO_NABO.msmc2.final.txt $DIR/MABO.msmc2.final.txt \
    NABO.msmc2.final.txt > MABO_NABO.combined.msmc2.final.txt



mu <- 1.25e-8
gen <- 10
crossPopDat<-read.table("results/EUR_AFR.combined.msmc2.final.txt", header=TRUE)
plot(crossPopDat$left_time_boundary/mu*gen, 2 * crossPopDat$lambda_01 / (crossPopDat$lambda_00 + crossPopDat$lambda_11),
     xlim=c(1000,500000),ylim=c(0,1), type="n", xlab="Years ago", ylab="relative cross-coalescence rate")
lines(crossPopDat$left_time_boundary/mu*gen, 2 * crossPopDat$lambda_01 / (crossPopDat$lambda_00 + crossPopDat$lambda_11), type="s")