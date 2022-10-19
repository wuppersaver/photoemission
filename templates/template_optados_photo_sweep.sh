#!/bin/bash  --login
#PBS -N TEMPLATE_od_photo_sweep
#PBS -l select=1:ncpus=1:mem=50GB
#PBS -l walltime=02:30:00

cd $PBS_O_WORKDIR

module load intel-suite

CONTINUE=false

CASE_IN=TEMPLATE

########### OptaDOS Photoemission ###########

OPTADOS=/rds/general/user/fcm19/home/modules_codes/optados/optados.x

sweep_values=()

cp ${CASE_IN}_optados_photo_sweep.odi ${CASE_IN}.odi

sed -i 's/.*photo_model.*/photo_model : 1step/' ${CASE_IN}.odi
CASE_OUT=${CASE_IN}_1step.out

for i in ${sweep_values[@]}
do
	sed -i "s/.*photon_energy.*/photon_energy : $i/" ${CASE_IN}.odi
	$OPTADOS $CASE_IN 2>&1 | tee -a $CASE_OUT
	mv ${CASE_IN}.odo ${CASE_IN}_${i}_1step.odo
done

sed -i 's/.*photo_model.*/photo_model : 3step/' ${CASE_IN}.odi
CASE_OUT=${CASE_IN}_3step.out

for i in ${sweep_values[@]}
do
	sed -i "s/.*photon_energy.*/photon_energy : $i/" ${CASE_IN}.odi
	$OPTADOS $CASE_IN 2>&1 | tee -a $CASE_OUT
	mv ${CASE_IN}.odo ${CASE_IN}_${i}_3step.odo
done

exit_code=$?
if [ "$CONTINUE" == true ]; then
    echo $exit_code
    if [[ $exit_code -eq 0 ]] ; then
        sed -i '0,/.*STATE=.*/s//STATE=od_photo_sweep_success/' ${CASE_IN}_submission.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/' ${CASE_IN}_optados_photo_sweep.sh
        ./${CASE_IN}_submission.sh
    else
        sed -i '0,/.*STATE=.*/s//STATE=od_photo_sweep_fail/' ${CASE_IN}_submission.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/' ${CASE_IN}_optados_photo_sweep.sh
    fi
fi
