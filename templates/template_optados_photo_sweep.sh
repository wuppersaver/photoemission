#!/bin/bash  --login
#PBS -N TEMPLATE_od_photo_sweep
#PBS -l select=1:ncpus=1:mem=50GB
#PBS -l walltime=02:30:00
#PBS -j oe

cd $PBS_O_WORKDIR

cleanup() {
    cd $PBS_O_WORKDIR
    mkdir ./tmp_$PBS_JOBID/
    cp $TMPDIR/* ./tmp_$PBS_JOBID/
}

exithandler() {
    echo "Job was killed on $date" | tee -a $TMPDIR/$CASE_IN.err
    cleanup
    exit
}

trap exithandler SIGTERM

echo "<qsub_standard_output>"
echo Start Date and Time
date

# Make sure any symbolic links are resolved to absolute path
export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)

# Set the number of threads to 1
#   This prevents any system libraries from automatically
#   using threading.
# export OMP_NUM_THREADS=1
env
echo "</qsub_standard_output>"

#to sync nodes
cd $PBS_O_WORKDIR

module load intel-suite

CONTINUE=false
CASE_IN=TEMPLATE
models=(___)
energies=($(seq -f "%'.2f" ___ *** ---))
jdos_maxs=(25) #($(seq -f "%'.2f" 5 1 20))
jdos_spacings=(0.1) #($(seq -f "%'.2f" 5 1 20))
iprint=___
directory=___

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
                sed -i "s/.*iprint.*/iprint : $iprint/" ${CASE_IN}.odi
                sed -i "s/.*photo_model.*/photo_model : $model/" ${CASE_IN}.odi
                sed -i "s/.*jdos_spacing.*/jdos_spacing : $jdos_space/" ${CASE_IN}.odi
                sed -i "s/.*jdos_max_energy.*/jdos_max_energy : $jdos_max/" ${CASE_IN}.odi
                sed -i "s/.*photo_photon_energy.*/photo_photon_energy : $energy/" ${CASE_IN}.odi
                CASE_OUT=${directory}${CASE_IN}.out
                $OPTADOS $CASE_IN 2>&1 | tee -a $CASE_OUT
                mv ${CASE_IN}.odo ${directory}${CASE_IN}_${energy}_${model}_${jdos_space}_${jdos_max}_${iprint}.odo
            done
        done 
    done
done

exit_code=$?
if [ "$CONTINUE" -eq true ]; then
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
