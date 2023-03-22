#!/bin/bash  --login
#PBS -N TEMPLATE_od_photo_sweep
#PBS -l select=1:ncpus=1:mem=50GB
#PBS -l walltime=02:30:00

cd $PBS_O_WORKDIR

module load intel-suite

CONTINUE=false
instance=1
CASE_IN=TEMPLATE
models=(___)
energies=($(seq -f "%'.2f" ___ *** ---))
jdos_maxs=(25) #($(seq -f "%'.2f" 5 1 20))
jdos_spacings=(0.1) #($(seq -f "%'.2f" 5 1 20))
iprint=___
directory=___
filename=${CASE_IN}_${instance}.odi
########### OptaDOS Photoemission ###########

OPTADOS=~/modules_codes/optados_photo_dev/optados/optados.x

cp ${CASE_IN}_optados_photo.odi ${CASE_IN}.odi

if [ "${directory}" != './' ]; then
    if [ ! -d $directory ]; then
        mkdir $directory
    fi
fi

for energy in ${energies[@]}
do
    for model in ${models[@]}
    do
        for jdos_max in ${jdos_maxs[@]}
        do
            for jdos_space in ${jdos_spacings[@]}
            do
                sed -i "s/.*iprint.*/iprint : $iprint/" $filename
                sed -i "s/.*photo_model.*/photo_model : $model/" $filename
                sed -i "s/.*jdos_spacing.*/jdos_spacing : $jdos_space/" $filename
                sed -i "s/.*jdos_max_energy.*/jdos_max_energy : $jdos_max/" $filename
                sed -i "s/.*photo_photon_energy.*/photo_photon_energy : $energy/" $filename
                CASE_OUT=${directory}${CASE_IN}.out
                $OPTADOS -multi_out $instance $CASE_IN 2>&1 | tee -a $CASE_OUT
                #mv ${CASE_IN}.odo ${directory}${CASE_IN}_${energy}_${model}_${jdos_space}_${jdos_max}_${iprint}.odo
        done 
    done
done

exit_code=$?
if [ "$CONTINUE" == 'true' ]; then
    echo $exit_code
    if [ $exit_code -eq 0 ] ; then
        sed -i '0,/.*STATE=.*/s//STATE=od_photo_sweep_success/' ${CASE_IN}_submission.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/' ${CASE_IN}_optados_photo_sweep.sh
        ./${CASE_IN}_submission.sh
        exit
    else
        sed -i '0,/.*STATE=.*/s//STATE=od_photo_sweep_fail/' ${CASE_IN}_submission.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/' ${CASE_IN}_optados_photo_sweep.sh
        exit
    fi
fi
