#!/bin/bash  --login
#PBS -N TEMPLATE_workfct
#PBS -l select=1:ncpus=4:mem=20GB
#PBS -l walltime=00:10:00
#PBS -j oe

cd $PBS_O_WORKDIR
echo "<qsub_standard_output>"
date
echo "<qstat -f $PBS_JOBID>"
qstat -f $PBS_JOBID
echo "</qstat -f $PBS_JOBID>"

# Make sure any symbolic links are resolved to absolute path
export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)

# Set the number of threads to 1
#   This prevents any system libraries from automatically
#   using threading.
# export OMP_NUM_THREADS=1
env
echo "</qsub_standard_output>"

#to sync nodes
cd $PBS_O_WORKDIR

module load anaconda3/personal
source activate matchem

CONTINUE=false

CASE_IN=TEMPLATE

########### Python Script ###########

#sed -i 's,input_path =.*,'"input_path = \'$(pwd)\/\'"',' ~/PhD/photoemission/function_scripts/calculate_photo_properties.py
file_ending=___
# ARGS: directory file_ending mod_odi(bool) centered_structure?(bool)
python ~/PhD/photoemission/function_scripts/calculate_photo_properties.py $(pwd) $file_ending *** ___

exit_code=$?
echo the_exit_code=$exit_code

if [ "$CONTINUE" -eq true ]; then
    if [ $exit_code -eq 0 ]; then
        sed -i '0,/.*STATE=.*/s//STATE=workfct_success/' ${CASE_IN}_submission.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/' ${CASE_IN}_workfunction.sh
        ./${CASE_IN}_submission.sh
        exit
    else
        sed -i '0,/.*STATE=.*/s//STATE=workfct_fail/' ${CASE_IN}_submission.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/' ${CASE_IN}_workfunction.sh
        exit
    fi
fi
