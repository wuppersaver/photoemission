CASE_IN=TEMPLATE

STATE=new
INTERNAL=$STATE

declare -A calculation
calculation[geometry]=${CASE_IN}_geometry.sh
calculation[spectral]=${CASE_IN}_spectral.sh
calculation[bands]=${CASE_IN}_bands.sh
calculation[optados_all]=${CASE_IN}_optados_fermi.sh
calculation[work_fct]=${CASE_IN}_workfunction.sh
calculation[optados_photo_sweep]=${CASE_IN}_optados_photo_sweep.sh
calculation[optados_photo]=${CASE_IN}_optados_photo.sh
calculation[submission]=${CASE_IN}_submission.sh

if [[ $# -eq 0 ]] ; then #checking, if arguments are present in the bash call
    if [[ $INTERNAL -eq *stop* ]]; then
        exit
    fi
    if [[ $INTERNAL -eq new ]]; then
