#!/bin/bash  --login
#PBS -N TEMPLATE_geom
#PBS -l select=1:ncpus=64:mem=200GB
#PBS -l walltime=06:00:00

cd $PBS_O_WORKDIR

module load mpi intel-suite

CONTINUE=false

CASE_IN=TEMPLATE

########### Geometry Optimization ###########

PRGM=~/modules_codes/CASTEP-18.1_orig/obj/linux_x86_64_ifort17--mpi/castep.mpi

cp ${CASE_IN}_geometry.cell ${CASE_IN}.cell
cp ${CASE_IN}_geometry.param ${CASE_IN}.param

CASE_OUT=${CASE_IN}.out

mpiexec $PRGM $CASE_IN 2>&1 | tee $CASE_OUT
exit_code=$?
echo the_exit_code=$exit_code
if [ "$CONTINUE" == true ]; then
    if [[ exit_code == 0 ]] ; then
        sed -i '0,/.*STATE=.*/s//STATE=geometry_success/' ${CASE_IN}_submission.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/' ${CASE_IN}_geometry.sh
        ./${CASE_IN}_submission.sh
    else
        sed -i '0,/.*STATE=.*/s//STATE=geometry_fail/' ${CASE_IN}_submission.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/' ${CASE_IN}_geometry.sh
    fi
fi
