CASE_IN=Cu_surf_111

FILE_CASTEP=${CASE_IN}.castep

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

