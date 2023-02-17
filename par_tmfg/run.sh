#!/bin/bash

datasets=(
"Mallat"
"UWaveGestureLibraryAll"
"NonInvasiveFetalECGThorax2"
"MixedShapesRegularTrain"
"MixedShapesSmallTrain"
"ECG5000"
"NonInvasiveFetalECGThorax1"
"StarLightCurves"
"HandOutlines"
"UWaveGestureLibraryX"
"CBF"
"InsectWingbeatSound"
"UWaveGestureLibraryY"
"ShapesAll"
"SonyAIBORobotSurface2"
"FreezerSmallTrain"
"Crop"
"ElectricDevices"
)


sizes=(
    2400
    4478
    3765
    2925
    2525
    5000
    3765
    9236
    1370
    4478
    930
    2200
    4478
    1200
    980
    2878
    19412
    16160
)

OUTPUTDIR="outputs/cdbht"
INPUTDIR="~/datasets/UCR"


workers=(96 48 36 24 12 4 1) 
round=3
prefixs=(2 5 10 30 50 200)

# OUTPUTDIR="outputs"
# INPUTDIR="../datasets"

[ -d ${OUTPUTDIR} ] || mkdir ${OUTPUTDIR}
[ -d "outputs/Ps" ] || mkdir "outputs/Ps"
[ -d "outputs/Zs" ] || mkdir "outputs/Zs"

for wk in "${workers[@]}"; do
    ind=0
    for dataset in "${datasets[@]}"; do
    for prefix in "${prefixs[@]}"; do
	    if [[ "${wk}" -eq 1 ]];then
            command="PARLAY_NUM_THREADS=${wk} ./tmfg ${INPUTDIR}/${dataset}.dat ${dataset} ${sizes[$ind]} 0 prefix ${prefix} ${round} > ${OUTPUTDIR}/${dataset}_prefix_${prefix}_${wk}th.txt"
        else
            command="PARLAY_NUM_THREADS=${wk} numactl -i all ./tmfg ${INPUTDIR}/${dataset}.dat ${dataset} ${sizes[$ind]} 0 prefix ${prefix} ${round} > ${OUTPUTDIR}/${dataset}_prefix_${prefix}_${wk}th.txt"
        fi
        echo "$command"
        eval "$command"   
    done
    let ind++
    done
done


for wk in "${workers[@]}"; do
    ind=0
    for dataset in "${datasets[@]}"; do
	    if [[ "${wk}" -eq 1 ]];then
            command="PARLAY_NUM_THREADS=${wk} ./tmfg ${INPUTDIR}/${dataset}.dat ${dataset} ${sizes[$ind]} 0 exact 0 ${round} > ${OUTPUTDIR}/${dataset}_exact_${wk}th.txt"
        else
            command="PARLAY_NUM_THREADS=${wk} numactl -i all ./tmfg ${INPUTDIR}/${dataset}.dat ${dataset} ${sizes[$ind]} 0 exact 0 ${round} > ${OUTPUTDIR}/${dataset}_exact_${wk}th.txt"
        fi
        echo "$command"
        eval "$command"
    let ind++
    done
done
echo "done"