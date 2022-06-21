CASE_IN=template

STATE=new
INTERNAL=$STATE
FILE_CASTEP=${CASE_IN}.castep

declare -A calculation
calculation[geometry]=${CASE_IN}_geometry.sh
calculation[spectral]=${CASE_IN}_spectral.sh
calculation[bands]=${CASE_IN}_bands.sh
calculation[optados_all]=${CASE_IN}_od_all.sh
calculation[work_fct]=${CASE_IN}_set_work_fct.sh
calculation[optados_photo]=${CASE_IN}_od_sweep_4-5.5.sh

if [[ $# -eq 0 ]] ; then #checking, if arguments are present in the bash call
    if [[ $INTERNAL == *stop* ]]; then
        exit
    fi
    if [[ $INTERNAL == new ]]; then    
        sed -i '0,/.*STATE=.*/s//STATE=geometry_run/' ${CASE_IN}_subm.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=true/'  ${calculation[geometry]} 
        echo "GeometryOptimization"
        qsub ${calculation[geometry]} 
        exit
    fi
    if [[ $INTERNAL == geometry_run ]]; then    
        sed -i '0,/.*STATE=.*/s//STATE=geometry_cont/' ${CASE_IN}_subm.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=true/'  ${calculation[geometry]} 
        echo "CONTINUATION" | tee -a ${CASE_IN}_geom.param
        echo "GeometryOptimization Continued"
        qsub ${calculation[geometry]}
        exit
    fi
    if [[ $INTERNAL == geometry_success ]]; then
        TEMP=$(tail ${FILE_CASTEP})
        if [[ $TEMP == *"Total time"* ]]; then
            sed -i '0,/.*STATE=.*/s//STATE=spectral_run/' ${CASE_IN}_subm.sh
            sed -i '0,/.*CONTINUE=.*/s//CONTINUE=true/'  ${calculation[spectral]} 
            echo "Spectral"
            qsub ${calculation[spectral]}
            exit
        fi
    fi
    if [[ $INTERNAL == spectral_success ]]; then
        sed -i '0,/.*STATE=.*/s//STATE=od_all_run/' ${CASE_IN}_subm.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=true/'  ${calculation[optados_all]} 
        echo "OptaDOS misc"
        qsub ${calculation[optados_all]}
        exit
    fi
    if [[ $INTERNAL == od_all_success ]]; then
        sed -i '0,/.*STATE=.*/s//STATE=bands_run/' ${CASE_IN}_subm.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=true/'  ${calculation[bands]} 
        echo "Bandstructure"
        qsub ${calculation[bands]}
        exit
    fi
    if [[ $INTERNAL == bands_success ]]; then
        sed -i '0,/.*STATE=.*/s//STATE=wrkfct_run/' ${CASE_IN}_subm.sh
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=true/'  ${calculation[work_fct]}
        echo "Setting Workfct and Volume/Area"
        qsub ${calculation[work_fct]]}
        exit
    fi
    if [[ $INTERNAL == wrkfct_success ]]; then
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
            sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/'  ${calculation[$1]} 
            qsub ${calculation[$1]}
            exit
        fi
    else
        echo "Too many arguments given!!"
        exit 1
    fi
fi


