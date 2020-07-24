#!/bin/bash
#ALTERAR O PATH SE MUDAR DE DIRETÒRIO
cd "${0%/*}"

NExec=${1}
sample=${2}

#echo Calculating Mean and Std of execution ${sample}.
cd statistic/
./ConfluenceMean ${NExec} ${sample}
cd ..