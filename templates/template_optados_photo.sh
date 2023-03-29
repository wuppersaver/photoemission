#!/bin/bash  --login
#PBS -N TEMPLATE_od_photo
#PBS -l select=1:ncpus=1:mem=50GB
#PBS -l walltime=02:30:00
#PBS -j oe

cd $PBS_O_WORKDIR
echo "<qsub_standard_output>"
date
echo "<qstat -f $PBS_JOBID>"
qstat -f $PBS_JOBID
echo "</qstat -f $PBS_JOBID>"

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
models=___
energy=___
CASE_IN=TEMPLATE

########### OptaDOS Photoemission ###########

OPTADOS=~/modules_codes/optados_photo_dev/optados/optados.x

cp ${CASE_IN}_optados_photo.odi ${CASE_IN}.odi

CASE_OUT=${CASE_IN}_${energy}_${model}.out

for model in ${models[@]}
do
    sed -i "s/.*photo_model.*/photo_model : $model/" ${CASE_IN}.odi
    $OPTADOS $CASE_IN 2>&1 | tee -a $CASE_OUT
    exit_code=$?

    mv ${CASE_IN}.odo ${CASE_IN}_${energy}_${model}.odo
done

if [ "$CONTINUE" -eq true ]; then
    echo $exit_code
    if [ $exit_code -eq 0 ]; then
        sed -i '0,/.*STATE=.*/s//STATE=od_photo_success/' ${CASE_IN}_submission.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/' ${CASE_IN}_optados_photo.sh
        ./${CASE_IN}_submission.sh
        exit
    else
        sed -i '0,/.*STATE=.*/s//STATE=od_photo_fail/' ${CASE_IN}_submission.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/' ${CASE_IN}_optados_photo.sh
        exit
    fi
fi
