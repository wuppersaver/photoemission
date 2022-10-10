CASE_IN=TEMPLATE

STATE=new
INTERNAL=$STATE
FILE_CASTEP=${CASE_IN}.castep

declare -A calculation
calculation[geometry]=${CASE_IN}_geometry.sh
calculation[spectral]=${CASE_IN}_spectral.sh
calculation[bands]=${CASE_IN}_bands.sh
calculation[optados_all]=${CASE_IN}_optados_fermi.sh
calculation[work_fct]=${CASE_IN}_workfunction.sh
calculation[optados_photo]=${CASE_IN}_od_sweep_4-5.5.sh
calculation[submission]=${CASE_IN}_submission.sh

if [[ $# -eq 0 ]] ; then #checking, if arguments are present in the bash call
    if [[ $INTERNAL == *stop* ]]; then
        exit
    fi
    if [[ $INTERNAL == new ]]; then
