#!/bin/bash  --login
#PBS -N TEMPLATE_geom
#PBS -l select=1:ncpus=64:mem=200GB
#PBS -l walltime=06:00:00

cd $PBS_O_WORKDIR

module load mpi intel-suite

CONTINUE=false

CASE_IN=TEMPLATE

########### Geometry Optimization ###########

PRGM=~/modules_codes/CASTEP-18.1/obj/linux_x86_64_ifort17/castep.mpi

cp ${CASE_IN}_geom.cell ${CASE_IN}.cell
cp ${CASE_IN}_geom.param ${CASE_IN}.param

CASE_OUT=${CASE_IN}.out

mpiexec $PRGM $CASE_IN 2>&1 | tee $CASE_OUT

if [ "$CONTINUE" == true ]; then
    if [[ $? == 0 ]] ; then
        sed -i '0,/.*STATE=.*/s//STATE=geometry_success/' ${CASE_IN}_submission.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/' ${CASE_IN}_od_all.sh
        ./${CASE_IN}_submission.sh
    else
        sed -i '0,/.*STATE=.*/s//STATE=geometry_fail/' ${CASE_IN}_submission.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/' ${CASE_IN}_od_all.sh
    fi
fi
