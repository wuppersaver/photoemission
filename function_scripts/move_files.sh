pwd

OLD=Cu_surf_111

NEW=Cu_surf_100
SHORT=50A_Cu100
DIRECTORY=${NEW}_victor_50A

cp ./${OLD}_geom.cell ./${NEW}_geom.cell
sed -i "s/CASE_IN=.*/CASE_IN=${NEW}/" ${NEW}_geom.cell
mv ./${NEW}_geom.cell ../${DIRECTORY}/${NEW}_geom.cell

cp ./${OLD}_geom.param ./${NEW}_geom.param
sed -i "s/CASE_IN=.*/CASE_IN=${NEW}/" ${NEW}_geom.param
mv ${NEW}_geom.param ../${DIRECTORY}/${NEW}_geom.param

cp ./${OLD}_spec.cell ./${NEW}_spec.cell
sed -i "s/CASE_IN=.*/CASE_IN=${NEW}/" ${NEW}_spec.cell
mv ./${NEW}_spec.cell ../${DIRECTORY}/${NEW}_spec.cell

cp ./${OLD}_spec.param ./${NEW}_spec.param
sed -i "s/CASE_IN=.*/CASE_IN=${NEW}/" ${NEW}_spec.param
mv ./${NEW}_spec.param ../${DIRECTORY}/${NEW}_spec.param

cp ./${OLD}_bands.cell ./${NEW}_bands.cell
sed -i "s/CASE_IN=.*/CASE_IN=${NEW}/" ${NEW}_bands.cell
mv ./${NEW}_bands.cell ../${DIRECTORY}/${NEW}_bands.cell

cp ./${OLD}_bands.param ./${NEW}_bands.param
sed -i "s/CASE_IN=.*/CASE_IN=${NEW}/" ${NEW}_bands.param
mv ./${NEW}_bands.param ../${DIRECTORY}/${NEW}_bands.param

cp ./${OLD}_fermi.odi ../${DIRECTORY}/${NEW}_fermi.odi
cp ./${OLD}_photo.odi ../${DIRECTORY}/${NEW}_photo.odi

cp ./Castep_geom_opti.qsub ./${NEW}_geom_opti.qsub
sed -i "s/CASE_IN=.*/CASE_IN=${NEW}/" ${NEW}_geom_opti.qsub
sed -i "s/#PBS -N.*/#PBS -N ${SHORT}_geom/" ${NEW}_geom_opti.qsub
mv ./${NEW}_geom_opti.qsub ../${DIRECTORY}/${NEW}_geom_opti.qsub

cp ./Castep_spectral.qsub ./${NEW}_spectral.qsub
sed -i "s/CASE_IN=.*/CASE_IN=${NEW}/" ${NEW}_spectral.qsub
sed -i "s/#PBS -N.*/#PBS -N ${SHORT}_spectral/" ${NEW}_spectral.qsub
mv ./${NEW}_spectral.qsub ../${DIRECTORY}/${NEW}_spectral.qsub 

cp ./OptaDOS_sweep_4-5.5.qsub ../${DIRECTORY}/OptaDOS_sweep_4-5.5.qsub
sed -i "s/#PBS -N.*/#PBS -N ${SHORT}_od_sweep/" ../${DIRECTORY}/OptaDOS_sweep_4-5.5.qsub
sed -i "s/CASE_IN=.*/CASE_IN=${NEW}/" ../${DIRECTORY}/OptaDOS_sweep_4-5.5.qsub