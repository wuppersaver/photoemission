#!/bin/bash  --login
#PBS -N template_od_photo
#PBS -l select=1:ncpus=1:mem=50GB
#PBS -l walltime=02:30:00

cd $PBS_O_WORKDIR

module load intel-suite

CONTINUE=false

CASE_IN=template

########### OptaDOS Photoemission ###########

OPTADOS=/rds/general/user/fcm19/home/modules_codes/optados/optados.x

energies=(4.0 4.07895 4.15789 4.23684 4.31579 4.39474 4.47368 4.55263 4.63158 4.71053 4.78947 4.86842 4.94737 5.02632 5.10526 5.18421 5.26316 5.34211 5.42105 5.5)

cp ${CASE_IN}_photo.odi ${CASE_IN}.odi

sed -i 's/.*photo_model.*/photo_model : 1step/' ${CASE_IN}.odi
CASE_OUT=${CASE_IN}_1step.out

for i in ${energies[@]}
do
	sed -i "s/.*photon_energy.*/photon_energy : $i/" ${CASE_IN}.odi
	$OPTADOS $CASE_IN 2>&1 | tee -a $CASE_OUT
	mv ${CASE_IN}.odo ${CASE_IN}_${i}_1step.odo
done

sed -i 's/.*photo_model.*/photo_model : 3step/' ${CASE_IN}.odi
CASE_OUT=${CASE_IN}_3step.out

for i in ${energies[@]}
do
	sed -i "s/.*photon_energy.*/photon_energy : $i/" ${CASE_IN}.odi
	$OPTADOS $CASE_IN 2>&1 | tee -a $CASE_OUT
	mv ${CASE_IN}.odo ${CASE_IN}_${i}_3step.odo
done

exit_code=$?
if [ "$CONTINUE" == true ]; then
    echo $exit_code
    if [[ exit_code == 0 ]] ; then
        sed -i '0,/.*STATE=.*/s//STATE=od_photo_success/' ${CASE_IN}_subm.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/' ${CASE_IN}_od_sweep_4-5.5.sh
        ./${CASE_IN}_subm.sh
    else
        sed -i '0,/.*STATE=.*/s//STATE=od_photo_fail/' ${CASE_IN}_subm.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/' ${CASE_IN}_od_sweep_4-5.5.sh
    fi
fi
