CASE_IN=Cu_surf_111

FILE_CASTEP=${CASE_IN}.castep

declare -A single_calc
single_calc[geometry]=${CASE_IN}_geom_opti.qsub
single_calc[spectral]=${CASE_IN}_spectral.qsub
single_calc[optados]=${CASE_IN}_sweep_4-5.5.qsub



if [[ $# -eq 0 ]] ; then
    if [ -f ./${FILE_CASTEP} ]; then
        TEMP=$(tail ${FILE_CASTEP})
        if [[ -f ./${CASE_IN}.ome_bin && -f ./${CASE_IN}.fem_bin ]]; then
            qsub OptaDOS_sweep_4-5.5.qsub
            exit
        fi
        if [[ $TEMP == *"Total time"* ]]; then
            qsub ${CASE_IN}_spectral.qsub
            exit
        else
            echo "CONTINUATION" | tee -a ${CASE_IN}_geom.param
            qsub ${CASE_IN}_geom_opti.qsub
            exit
        fi
    else
        qsub ${CASE_IN}_geom_opti.qsub
        exit
    fi
else
    if [[ $# -gt 0 && $# -lt 2 ]] ; then
        if [ ! -v single_calc[$1] ]; then
            echo "Wrong argument was given!!"
            exit 1
        else
            qsub ${single_calc[$1]}
            exit
        fi
    else
        echo "Too many arguments given!!"
        exit 1
    fi
fi



