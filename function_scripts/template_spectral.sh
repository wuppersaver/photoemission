#!/bin/bash  --login
#PBS -N template_spectral
#PBS -l select=1:ncpus=64:mem=200GB
#PBS -l walltime=02:00:00

cd $PBS_O_WORKDIR

module load mpi intel-suite

CONTINUE=false

CASE_IN=template

########### Spectral Task ###########

PRGM=~/modules_codes/CASTEP-18.1/obj/linux_x86_64_ifort17/castep.mpi

cp ${CASE_IN}_spec.param ${CASE_IN}.param
cp ${CASE_IN}_spec.cell ${CASE_IN}.cell

CASE_OUT=${CASE_IN}_spec.out

mpiexec $PRGM $CASE_IN 2>&1 | tee -a $CASE_OUT


########### Modified Spectral Task ###########

PRGM=~/modules_codes/CASTEP-18.1_mod/obj/linux_x86_64_ifort17/castep.mpi

cp ${CASE_IN}_spec.param ${CASE_IN}.param
cp ${CASE_IN}_spec.cell ${CASE_IN}.cell

CASE_OUT=${CASE_IN}_spec.out

mpiexec $PRGM $CASE_IN 2>&1 | tee -a $CASE_OUT

if [ "$CONTINUE" == true ]; then
    if [[ $? == 0 ]] ; then
        sed -i '0,/.*STATE=.*/s//STATE=spectral_success/' ${CASE_IN}_subm.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/' ${CASE_IN}_spectral.sh
        ./${CASE_IN}_subm.sh
    else
        sed -i '0,/.*STATE=.*/s//STATE=spectral_fail/' ${CASE_IN}_subm.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/' ${CASE_IN}_spectral.sh
    fi
fi
