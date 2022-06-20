CASE_IN=Cu_surf_100

STATE=geometry_run
INTERNAL=$STATE
FILE_CASTEP=${CASE_IN}.castep

declare -A calculation
calculation[geometry]=${CASE_IN}_geometry.qsub
calculation[spectral]=${CASE_IN}_spectral.qsub
calculation[bands]=${CASE_IN}_bands.qsub
calculation[optados_all]=${CASE_IN}_od_all.qsub
calculation[optados_photo]=${CASE_IN}_od_sweep_4-5.5.qsub

if [[ $# -eq 0 ]] ; then #checking, if arguments are present in the bash call
    if [[ $INTERNAL == *stop* ]]; then
        exit
    fi
    if [[ $INTERNAL == new ]]; then    
        sed -i '0,/.*STATE=.*/s//STATE=geometry_run/' ${CASE_IN}_subm.sh
        echo "GeometryOptimization"
        qsub ${calculation[geometry]} 
        exit
    fi
    if [[ $INTERNAL == geometry_run ]]; then    
        sed -i '0,/.*STATE=.*/s//STATE=geometry_cont/' ${CASE_IN}_subm.sh
        echo "CONTINUATION" | tee -a ${CASE_IN}_geom.param
        echo "GeometryOptimization Continued"
        qsub ${calculation[geometry]}
        exit
    fi
    if [[ $INTERNAL == geometry_success ]]; then 
        TEMP=$(tail ${FILE_CASTEP})
        if [[ $TEMP == *"Total time"* ]]; then #checking, that the geometry optimisation finished successfully
            sed -i '0,/.*STATE=.*/s//STATE=spectral_run/' ${CASE_IN}_subm.sh
            echo "Spectral"
            qsub ${calculation[spectral]}
            exit
        fi
    fi
    if [[ $INTERNAL == spectral_success ]]; then
        sed -i '0,/.*STATE=.*/s//STATE=od_all_run/' ${CASE_IN}_subm.sh
        echo "OptaDOS misc"
        qsub ${calculation[optados_all]}
        exit
    fi
    if [[ $INTERNAL == od_all_success ]]; then
        sed -i '0,/.*STATE=.*/s//STATE=bands_run/' ${CASE_IN}_subm.sh
        echo "Bandstructure"
        qsub ${calculation[bands]}
        exit
    fi
    if [[ $INTERNAL == bands_success ]]; then   #checking, that both unmodified and modified spectral tasks were run
        sed -i '0,/.*STATE=.*/s//STATE=od_photo_run/' ${CASE_IN}_subm.sh
        echo "OptaDOS Photoemission"
        qsub ${calculation[optados_photo]}
        exit
fi
else
    if [[ $# -gt 0 && $# -lt 2 ]] ; then
        if [ ! -v calculation[$1] ]; then
            echo "Wrong argument was given!!"
            exit 1
        else
             sed -i "0,/.*STATE=.*/s//STATE=${$1}_stop/" ${CASE_IN}_subm.sh
            qsub ${calculation[$1]}
            exit
        fi
    else
        echo "Too many arguments given!!"
        exit 1
    fi
fi

