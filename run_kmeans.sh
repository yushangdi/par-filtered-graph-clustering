#!/bin/bash

datasets=("Mallat"
"UWaveGestureLibraryAll"
"NonInvasiveFetalECGThorax2"
"MixedShapesRegularTrain"
"MixedShapesSmallTrain"
"ECG5000"
"NonInvasiveFetalECGThorax1"
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

# workers=(24) 

# method="kmeans_spectral_quality"

# for wk in "${workers[@]}"; do
# for dataset in "${datasets[@]}"; do
#     command="OMP_NUM_THREADS=${wk} MKL_NUM_THREADS=${wk} OPENBLAS_NUM_THREADS=${wk} BLIS_NUM_THREADS=${wk} python3 kmeans_timing.py ${wk} ${dataset} ${method}"
#     echo "$command"
#     eval "$command"
# done
# done

workers=(48 36 24 12 4 1) 

method="kmeans_spectral"

for wk in "${workers[@]}"; do
for dataset in "${datasets[@]}"; do
    command="OMP_NUM_THREADS=${wk} MKL_NUM_THREADS=${wk} OPENBLAS_NUM_THREADS=${wk} BLIS_NUM_THREADS=${wk} python3 kmeans_timing.py ${wk} ${dataset} ${method}"
    echo "$command"
    eval "$command"
done
done

workers=(96) # for workers=96, only use 24 threads for openblas because using too many threads give "BLAS : Program is Terminated. Because you tried to allocate too many memory regions."

method="kmeans_spectral"

for wk in "${workers[@]}"; do
for dataset in "${datasets[@]}"; do
    command="OMP_NUM_THREADS=${wk} MKL_NUM_THREADS=${wk} OPENBLAS_NUM_THREADS=24 BLIS_NUM_THREADS=${wk} python3 kmeans_timing.py ${wk} ${dataset} ${method}"
    echo "$command"
    eval "$command"
done
done