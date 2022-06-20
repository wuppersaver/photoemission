CASE_IN=Cu_surf_111

FILE_CASTEP=${CASE_IN}.castep

declare -A calculation
calculation[geometry]=${CASE_IN}_geom_opti.qsub
calculation[spectral]=${CASE_IN}_spectral.qsub
calculation[optados_photo]=${CASE_IN}_sweep_4-5.5.qsub

if [[ $# -eq 0 ]] ; then #checking, if arguments are present in the bash call
    if [ -f ./${FILE_CASTEP} ]; then #checking, that a calculation has been run before
        TEMP=$(tail ${FILE_CASTEP})
        if [[ -f ./${CASE_IN}.ome_bin && -f ./${CASE_IN}.fem_bin ]]; then   #checking, that both unmodified and modified spectral tasks were run
            echo "OptaDOS Photoemission"
			qsub ${calculation[optados_photo]}
            exit
        fi
        if [[ $TEMP == *"Total time"* ]]; then #checking, that the geometry optimisation finished successfully
			echo "Spectral"
            qsub ${calculation[spectral]}
            exit
        else
            echo "CONTINUATION" | tee -a ${CASE_IN}_geom.param
            echo "GeometryOptimization Continued"
			qsub ${calculation[geometry]}
            exit
        fi
    else
		echo "GeometryOptimization"
        qsub ${calculation[geometry]}
        exit
    fi
else
    if [[ $# -gt 0 && $# -lt 2 ]] ; then
        if [ ! -v calculation[$1] ]; then
            echo "Wrong argument was given!!"
            exit 1
        else
            qsub ${calculation[$1]}
            exit
        fi
    else
        echo "Too many arguments given!!"
        exit 1
    fi
fi

