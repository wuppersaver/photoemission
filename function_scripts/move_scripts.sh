NEW=Cu100
DIRECTORY=../structures/${NEW}_victor_60A

files=("subm" "geometry" "spectral" "od_all" "set_work_fct" "od_sweep_4-5.5" "bands")

if [ ! -d ${DIRECTORY} ]; then
    mkdir ${DIRECTORY}
fi

for i in ${files[@]}; do
	cp ./template_${i}.sh ${NEW}_${i}.sh
    sed -i "s/CASE_IN=.*/CASE_IN=${NEW}/" ${NEW}_${i}.sh
    sed -i "s/#PBS -N.*/#PBS -N ${NEW}_${i}/" ${NEW}_${i}.sh
    mv ./${NEW}_${i}.sh ${DIRECTORY}/${NEW}_${i}.sh
done

cwd=$(pwd)
cd ${DIRECTORY}
chmod 744 *.sh
cd ${cwd}