#!/bin/bash



# datasets=("Mallat"
# "MixedShapesRegularTrain"
# "MixedShapesSmallTrain"
# "MoteStrain"
# "HandOutlines"
# "CBF"
# "InsectWingbeatSound"
# "ShapesAll"
# "SonyAIBORobotSurface2"
# "FreezerSmallTrain"
# "UWaveGestureLibraryAll"
# "NonInvasiveFetalECGThorax2"
# "ECG5000"
# "NonInvasiveFetalECGThorax1"
# "UWaveGestureLibraryX"
# "UWaveGestureLibraryY"
# )

datasets=("iris")
OUTPUTDIR="outputs"
INPUTDIR="../datasets/"
MATLAB="matlab"

[ -d ${OUTPUTDIR} ] || mkdir ${OUTPUTDIR}
for dataset in "${datasets[@]}"; do
    command="${MATLAB} -nojvm -nosplash -nodesktop  -nodisplay  -r  'UCR_PMFG(\"${dataset}\", \"${INPUTDIR}\", \"${OUTPUTDIR}\"); exit'  -logfile ${OUTPUTDIR}/${dataset}_pmfg_timing.txt"
    echo "$command"
    eval "$command"
done

for dataset in "${datasets[@]}"; do
    command="${MATLAB} -nojvm -nosplash -nodesktop  -nodisplay  -r  'UCR_TMFG(\"${dataset}\", \"${INPUTDIR}\", \"${OUTPUTDIR}\"); exit'  -logfile ${OUTPUTDIR}/${dataset}_tmfg_timing.txt"
    echo "$command"
    eval "$command"
done