#!/bin/bash  --login
#PBS -N TEMPLATE_od_fermi
#PBS -l select=1:ncpus=1:mem=100GB
#PBS -l walltime=01:00:00

cd $PBS_O_WORKDIR

module load mpi intel-suite

CONTINUE=false

CASE_IN=TEMPLATE

########### OptaDOS ALL ###########

PRGM=~/modules_codes/optados_photo_dev/optados/optados.x

cp ${CASE_IN}_optados_fermi.odi ${CASE_IN}.odi

CASE_OUT=${CASE_IN}_od_fermi.out

$PRGM $CASE_IN 2>&1 | tee -a $CASE_OUT
exit_code=$?
mv ${CASE_IN}.odo ${CASE_IN}_fermi.odo

echo the_exit_code=$exit_code

if [ "$CONTINUE" -eq true ]; then
    if [[ $exit_code -eq 0 ]] ; then
        sed -i '0,/.*STATE=.*/s//STATE=od_fermi_success/' ${CASE_IN}_submission.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/' ${CASE_IN}_optados_fermi.sh
        ./${CASE_IN}_submission.sh
        exit
    else
        sed -i '0,/.*STATE=.*/s//STATE=od_fermi_fail/' ${CASE_IN}_submission.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/' ${CASE_IN}_optados_fermi.sh
        exit
    fi
fi
