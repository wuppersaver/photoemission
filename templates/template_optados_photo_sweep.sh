#!/bin/bash  --login
#PBS -N TEMPLATE_od_photo_sweep
#PBS -l select=1:ncpus=1:mem=50GB
#PBS -l walltime=02:30:00

cd $PBS_O_WORKDIR

module load intel-suite

CONTINUE= false
CASE_IN=TEMPLATE
models=(___)
energies=($(seq -f "%'.2f" ___ *** ---))
jdos_maxs=(25) #($(seq -f "%'.2f" 5 1 20))
jdos_spacings=(0.1) #($(seq -f "%'.2f" 5 1 20))
iprint=___
directory=___

########### OptaDOS Photoemission ###########

OPTADOS=~/modules_codes/optados_photo_dev/optados/volume/optados.x

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
                sed -i "s/.*iprint.*/iprint : $iprint/" ${CASE_IN}.odi
                sed -i "s/.*photo_model.*/photo_model : $model/" ${CASE_IN}.odi
                sed -i "s/.*JDOS_SPACING.*/JDOS_SPACING : $jdos_space/" ${CASE_IN}.odi
                sed -i "s/.*JDOS_MAX_ENERGY.*/JDOS_MAX_ENERGY : $jdos_max/" ${CASE_IN}.odi
                sed -i "s/.*photo_photon_energy.*/photo_photon_energy : $energy/" ${CASE_IN}.odi
                CASE_OUT=${directory}${CASE_IN}.out
                $OPTADOS $CASE_IN 2>&1 | tee -a $CASE_OUT
                mv ${CASE_IN}.odo ${directory}${CASE_IN}_${energy}_${model}_${jdos_space}_${jdos_max}_${iprint}.odo
        done 
    done
done

exit_code=$?
if CONTINUE; then
    echo $exit_code
    if [[ $exit_code -eq 0 ]] ; then
        sed -i '0,/.*STATE=.*/s//STATE=od_photo_sweep_success/' ${CASE_IN}_submission.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE= false/' ${CASE_IN}_optados_photo_sweep.sh
        ./${CASE_IN}_submission.sh
        exit
    else
        sed -i '0,/.*STATE=.*/s//STATE=od_photo_sweep_fail/' ${CASE_IN}_submission.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE= false/' ${CASE_IN}_optados_photo_sweep.sh
        exit
    fi
fi
