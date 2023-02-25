#!/bin/bash  --login
#PBS -N TEMPLATE_spectral
#PBS -l select=1:ncpus=64:mem=200GB
#PBS -l walltime=02:00:00

cd $PBS_O_WORKDIR

module load mpi intel-suite

CONTINUE=false

CASE_IN=TEMPLATE

########### Spectral Task ###########

PRGM=~/modules_codes/CASTEP-18.1_orig/obj/linux_x86_64_ifort17--mpi/castep.mpi

cp ${CASE_IN}_spectral.param ${CASE_IN}.param
cp ${CASE_IN}_spectral.cell ${CASE_IN}.cell

CASE_OUT=${CASE_IN}_spectral.out

mpiexec $PRGM $CASE_IN 2>&1 | tee -a $CASE_OUT

cp ${CASE_IN}.bands ${CASE_IN}.bands.spec


########### Modified Spectral Task ###########

PRGM=~/modules_codes/CASTEP-18.1_ome_new/obj/linux_x86_64_ifort17/castep.mpi

cp ${CASE_IN}_spectral.param ${CASE_IN}.param
cp ${CASE_IN}_spectral.cell ${CASE_IN}.cell

CASE_OUT=${CASE_IN}_spectral.out

mpiexec $PRGM $CASE_IN 2>&1 | tee -a $CASE_OUT
exit_code=$?
cp ${CASE_IN}.bands ${CASE_IN}.bands.spec.mod
echo the_exit_code=$exit_code

if [ "$CONTINUE" == true ]; then
    if [[ $exit_code -eq 0 ]] ; then
        sed -i '0,/.*STATE=.*/s//STATE=spectral_success/' ${CASE_IN}_submission.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/' ${CASE_IN}_spectral.sh
        ./${CASE_IN}_submission.sh
    else
        sed -i '0,/.*STATE=.*/s//STATE=spectral_fail/' ${CASE_IN}_submission.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/' ${CASE_IN}_spectral.sh
    fi
fi
