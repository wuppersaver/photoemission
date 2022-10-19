#!/bin/bash  --login
#PBS -N TEMPLATE_workfct
#PBS -l select=1:ncpus=4:mem=20GB
#PBS -l walltime=00:10:00

cd $PBS_O_WORKDIR

module load anaconda3/personal
source activate matchem

CONTINUE=false

CASE_IN=TEMPLATE

########### Python Script ###########

sed -i 's,input_path =.*,'"input_path = \'$(pwd)\/\'"',' ~/PhD/photoemission/function_scripts/calculate_photo_properties.py

python ~/PhD/photoemission/function_scripts/calculate_photo_properties.py

exit_code=$?
echo the_exit_code=$exit_code

if [ "$CONTINUE" == true ]; then
    if [[ $exit_code -eq 0 ]] ; then
        sed -i '0,/.*STATE=.*/s//STATE=workfct_success/' ${CASE_IN}_submission.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/' ${CASE_IN}_workfunction.sh
        ./${CASE_IN}_submission.sh
    else
        sed -i '0,/.*STATE=.*/s//STATE=workfct_fail/' ${CASE_IN}_submission.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/' ${CASE_IN}_workfunction.sh
    fi
fi
