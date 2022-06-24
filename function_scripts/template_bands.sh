#!/bin/bash  --login
#PBS -N template_bands
#PBS -l select=1:ncpus=64:mem=200GB
#PBS -l walltime=08:00:00

cd $PBS_O_WORKDIR

module load mpi intel-suite

CONTINUE=false

CASE_IN=template

########### Bandstructure ###########

PRGMO2B=~/modules_codes/CASTEP-19.11/obj/linux_x86_64_ifort17--serial/orbitals2bands

PRGM=~/modules_codes/CASTEP-18.1/obj/linux_x86_64_ifort17/castep.mpi

cp ${CASE_IN}_bands.param ${CASE_IN}.param
cp ${CASE_IN}_bands.cell ${CASE_IN}.cell

CASE_OUT=${CASE_IN}_bands.out

mpiexec $PRGM $CASE_IN 2>&1 | tee -a $CASE_OUT

cp ${CASE_IN}.bands ${CASE_IN}.bands.orig

$PRGMO2B $CASE_IN 2>&1 | tee -a $CASE_OUT



if [ "$CONTINUE" == true ]; then
    if [[ $? == 0 ]] ; then
        sed -i '0,/.*STATE=.*/s//STATE=bands_success/' ${CASE_IN}_subm.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/' ${CASE_IN}_bands.sh
        ./${CASE_IN}_subm.sh
    else
        sed -i '0,/.*STATE=.*/s//STATE=bands_fail/' ${CASE_IN}_subm.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/' ${CASE_IN}_bands.sh
    fi
fi

