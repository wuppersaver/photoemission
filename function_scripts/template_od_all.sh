#!/bin/bash  --login
#PBS -N template_od_all
#PBS -l select=1:ncpus=1:mem=100GB
#PBS -l walltime=01:00:00

cd $PBS_O_WORKDIR

module load mpi intel-suite

CONTINUE=false

CASE_IN=template

########### OptaDOS ALL ###########

PRGM=~/modules_codes/optados/optados.x

cp ${CASE_IN}_all.odi ${CASE_IN}.odi

CASE_OUT=${CASE_IN}_od.out

$PRGM $CASE_IN 2>&1 | tee -a $CASE_OUT

mv ${CASE_IN}.odo ${CASE_IN}_all.odo

if [ "$CONTINUE" == true ]; then
    if [[ $? == 0 ]] ; then
        sed -i '0,/.*STATE=.*/s//STATE=od_all_success/' ${CASE_IN}_subm.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/' ${CASE_IN}_od_all.sh
        ./${CASE_IN}_subm.sh
    else
        sed -i '0,/.*STATE=.*/s//STATE=od_all_fail/' ${CASE_IN}_subm.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/' ${CASE_IN}_od_all.sh
    fi
fi
