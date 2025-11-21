#!/bin/bash
# This script is functionally the same as the one given in the main GONE repo. This is just edited to work in a SLURM job of the UofA HPC

echo "Running GONE"
echo "Running GONE" >> ../timefile
START=$(date +%s)


arg=("$@")
options_for_GONE=""
i=0;
while [ $i -lt $[$#-2] ]
do
    if [ "${arg[i]}" == "-lc" ] || [ "${arg[i]}" == "-hc" ] || [ "${arg[i]}" == "-ne" ] || [ "${arg[i]}" == "-sd" ] || [ "${arg[i]}" == "-bs" ] || [ "${arg[i]}" == "-bn" ]
    then
        i=$[$i+2]
    elif [ "${arg[i]}" == "-sr" ] || [ "${arg[i]}" == "-vb" ]
    then
        i=$[$i+1]
    else
        break
    fi
done

for ((j=0;j<i;j++))
do
options_for_GONE="${options_for_GONE} ${arg[j]}"
done

j=$(($#-$i))

if [ $j -eq 3 ]
then
    echo "eq 3"
    fichero=${arg[i]}
    repeats=${arg[i+1]}
    threads=$(getconf _NPROCESSORS_ONLN)
    FILE=${arg[i+2]}
elif [ $j -eq 4 ]
then
echo "eq 4"
fichero=${arg[i]}
repeats=${arg[i+1]}
threads=${arg[i+2]}
FILE=${arg[i+2]}
else
	echo "  Usage: ./GONEparallel.sh [options_for_GONE] <filename> <number_of_repeats> [number_of_threads]"
	exit 1
fi

dirtemp="${fichero}_TEMP"

if [ -d $dirtemp ]
then
	rm -r $dirtemp
fi
mkdir $dirtemp
for ((nr=1; nr<=$repeats; nr++)); do echo $nr; done | xargs -I % -P $threads bash -c "./PROGRAMMES/GONE $options_for_GONE $fichero ${dirtemp}/${fichero}_%"

./PROGRAMMES/GONEaverage $dirtemp $fichero $repeats
echo " END OF ALL PROCESSES"

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "GONE run took $DIFF seconds"
echo "GONE run took $DIFF seconds" >> ../timefile

echo "END OF ANALYSES"
echo "END OF ANALYSES" >> ../timefile

mv outfileLD_Ne_estimates Output_Ne_$FILE
mv outfileLD_d2_sample Output_d2_$FILE
rm outfileLD
rm data.ped
rm data.map
rm out
mv outfileLD_TEMP TEMPORARY_FILES/

###################################################################
