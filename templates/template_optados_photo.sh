#!/bin/bash  --login
#PBS -N TEMPLATE_od_photo
#PBS -l select=1:ncpus=1:mem=50GB
#PBS -l walltime=02:30:00

cd $PBS_O_WORKDIR

module load intel-suite

CONTINUE=false
models=___
energy=___
CASE_IN=TEMPLATE

########### OptaDOS Photoemission ###########

OPTADOS=/rds/general/user/fcm19/home/modules_codes/optados/optados.x

cp ${CASE_IN}_optados_photo.odi ${CASE_IN}.odi

CASE_OUT=${CASE_IN}_${energy}_${model}.out

for model in ${models[@]}
do
    sed -i "s/.*photo_model.*/photo_model : $model/" ${CASE_IN}.odi
    $OPTADOS $CASE_IN 2>&1 | tee -a $CASE_OUT
    exit_code=$?

    mv ${CASE_IN}.odo ${CASE_IN}_${energy}_${model}.odo
done

if [ "$CONTINUE" == true ]; then
    echo $exit_code
    if [[ $exit_code -eq 0 ]] ; then
        sed -i '0,/.*STATE=.*/s//STATE=od_photo_success/' ${CASE_IN}_submission.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/' ${CASE_IN}_optados_photo.sh
        ./${CASE_IN}_submission.sh
    else
        sed -i '0,/.*STATE=.*/s//STATE=od_photo_fail/' ${CASE_IN}_submission.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/' ${CASE_IN}_optados_photo.sh
    fi
fi
