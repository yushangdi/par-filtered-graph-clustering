#!/bin/bash

# datasets=("Mallat"
# "UWaveGestureLibraryAll"
# "NonInvasiveFetalECGThorax2"
# "MixedShapesRegularTrain"
# "MixedShapesSmallTrain"
# "ECG5000"
# "NonInvasiveFetalECGThorax1"
# "MoteStrain"
# "HandOutlines"
# "UWaveGestureLibraryX"
# "CBF"
# "InsectWingbeatSound"
# "UWaveGestureLibraryY"
# "ShapesAll"
# "SonyAIBORobotSurface2"
# "FreezerSmallTrain"
# )

# sizes=(
#     2400
#     4478
#     3765
#     2925
#     2525
#     5000
#     3765
#     1272
#     1370
#     4478
#     930
#     2200
#     4478
#     1200
#     980
#     2878
# )

# workers=(96 48 36 24 12 4 1) #  
# round=3

round=1

datasets=("iris-D")

sizes=(150)

workers=(1) #  

methods=(
    "comp"
    "avg"
)


INPUTDIR="../../datasets"
OUTPUTDIR="outputs"
[ -d ${OUTPUTDIR} ] || mkdir ${OUTPUTDIR}

make

for method in  "${methods[@]}"; do
for wk in "${workers[@]}"; do
    ind=0
for dataset in "${datasets[@]}"; do
    if [[ "${wk}" -eq 1 ]];then
        command="PARLAY_NUM_THREADS=${wk} ./linkage ${INPUTDIR}/${dataset}.csv ${sizes[$ind]} ${OUTPUTDIR}/${dataset}_${method}_dendro ${method} ${round} > ${OUTPUTDIR}/${dataset}_${method}_${wk}th_timing.txt"
    else
        command="PARLAY_NUM_THREADS=${wk} numactl -i all ./linkage ${INPUTDIR}/${dataset}.csv ${sizes[$ind]} ${OUTPUTDIR}/${dataset}_${method}_dendro ${method} ${round} > ${OUTPUTDIR}/${dataset}_${method}_${wk}th_timing.txt"
    fi
    echo "$command"
    eval "$command"
    let ind++
done
done
done