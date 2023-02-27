#!/bin/bash  --login
#PBS -N TEMPLATE_geom
#PBS -l select=1:ncpus=64:mem=250GB
#PBS -l walltime=06:00:00

cd $PBS_O_WORKDIR

module load mpi intel-suite

CONTINUE= false
STATE= false
TIGHT_CONV=$STATE
CASE_IN=TEMPLATE

########### Geometry Optimization ###########

PRGM=~/modules_codes/CASTEP-18.1_orig/obj/linux_x86_64_ifort17--mpi/castep.mpi

cp ${CASE_IN}_geometry.cell ${CASE_IN}.cell
cp ${CASE_IN}_geometry.param ${CASE_IN}.param

CASE_OUT=${CASE_IN}.out

mpiexec $PRGM $CASE_IN 2>&1 | tee $CASE_OUT
exit_code=$?

cp ${CASE_IN}.pot_fmt ${CASE_IN}.potfmt.geometry
cp ${CASE_IN}.den_fmt ${CASE_IN}.denfmt.geometry

echo the_exit_code=$exit_code
if CONTINUE; then
    if [[ $exit_code -eq 0 ]] ; then
        if ! test -e ${CASE_IN}.castep; then
            sed -i '0,/.*STATE=.*/s//STATE=geometry_fail/' ${CASE_IN}_submission.sh
            sed -i '0,/.*CONTINUE=.*/s//CONTINUE= false/' ${CASE_IN}_geometry.sh
            exit
        fi
        if grep -Fq "LBFGS: Geometry optimization completed successfully." ${CASE_IN}.castep; then
            if TIGHT_CONV; then
                sed -i '0,/.*STATE=.*/s//STATE=geometry_success/' ${CASE_IN}_submission.sh
                sed -i '0,/.*CONTINUE=.*/s//CONTINUE= false/' ${CASE_IN}_geometry.sh
                sed -i '0,/.*STATE=.*/s//STATE= false/' ${CASE_IN}_geometry.sh
                ./${CASE_IN}_submission.sh
                exit
            else
                mv ${CASE_IN}.castep ${CASE_IN}.castep.loose_conv
                sed -i '0,/.*STATE=.*/s//STATE=geometry_unfinished/' ${CASE_IN}_submission.sh
                sed -i '0,/.*ELEC_ENERGY_TOL:.*/s//ELEC_ENERGY_TOL: 1e-08 eV/' ${CASE_IN}_geometry.param
                sed -i '0,/.*STATE=.*/s//STATE= true/' ${CASE_IN}_geometry.sh
                sed -i '0,/.*CONTINUE=.*/s//CONTINUE= false/' ${CASE_IN}_geometry.sh
                ./${CASE_IN}_submission.sh
                exit
            fi
        else
            sed -i '0,/.*STATE=.*/s//STATE=geometry_unfinished/' ${CASE_IN}_submission.sh
            sed -i '0,/.*CONTINUE=.*/s//CONTINUE= false/' ${CASE_IN}_geometry.sh
            ./${CASE_IN}_submission.sh
            exit
        fi
    else
        sed -i '0,/.*STATE=.*/s//STATE=geometry_fail/' ${CASE_IN}_submission.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE= false/' ${CASE_IN}_geometry.sh
        exit
    fi
fi
