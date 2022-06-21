#!/bin/bash  --login
#PBS -N template_wrkfct
#PBS -l select=1:ncpus=4:mem=20GB
#PBS -l walltime=00:30:00

cd $PBS_O_WORKDIR

module load anaconda3/personal

CONTINUE=false

CASE_IN=template

########### Python Script ###########

source activate matchem

sed -i "0,/.*directory =.*/s//directory = '${pwd}'/" ~/function_scripts/calculate_photo_properties.py

python ~/function_scripts/calculate_photo_properties.py

if [ "$CONTINUE" == true ]; then
    if [[ $? == 0 ]] ; then
        sed -i '0,/.*STATE=.*/s//STATE=wrkfct_success/' ${CASE_IN}_subm.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/' ${CASE_IN}_set_work_fct.sh
        ./${CASE_IN}_subm.sh
    else
        sed -i '0,/.*STATE=.*/s//STATE=wrkfct_fail/' ${CASE_IN}_subm.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/' ${CASE_IN}_set_work_fct.sh
    fi
fi
