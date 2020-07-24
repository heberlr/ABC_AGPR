#!/bin/bash
rm confluece/*.dat *~ *# &> /dev/null
#ALTERAR O PATH SE MUDAR DE DIRETÒRIO
cd "${0%/*}"
export OMP_NUM_THREADS=8

NExec=46
sample=0
# p1=0.008150
# p2=0.018370
# p3=0.130500
# p4=14.931850
# p5=0.024610

#make

cd statistic
g++ statistical.cpp -o ConfluenceMean
cd ..

for rep in $(seq 1 ${NExec})
do
    ./build/main.exe ${rep} ${sample} ${1} ${2}
done

#echo Calculating Mean and Std of execution ${sample}.
cd statistic/
./ConfluenceMean ${NExec} ${sample}
cd ..
