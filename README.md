# par-filtered-graph-clustering

To get the submodules:
```bash
git pull
git submodule update --init
```


Implementations
--------
Each folder contains an implementation that we tested.

`par_tmfg` our parallel TMFG and DBHT implementation

`hac` the hierarchical agglomerative clustering (HAC) algorithm by [Yu et al.](https://arxiv.org/abs/2106.04727)

`Aste` Aste's MATLAB TMFG+DBHT implementation. This is modified from [DBHT](https://www.mathworks.com/matlabcentral/fileexchange/46750-dbht)
and [PMFG](https://www.mathworks.com/matlabcentral/fileexchange/38689-pmfg). The modifications include adding timers for benchmarking and substitute some subroutines for better performance. Speficically, we changed Aste's TMFG+DBHT implementation (DBHTs.m file) to use boost library's all pair's shortest path and breadth-first search implementation, because this gives significant speedup. Aste's MATLAB PMFG+DBHT implementation also uses boost's implementation.


## Installation Requirement

* g++ = 7.5.0 
* make
* [C++ boost library](https://www.boost.org/)
* MATLAB
* [MATLAB BGL](https://www.mathworks.com/matlabcentral/fileexchange/10922-matlabbgl)

After boost is installed, set the BOOST_ROOT variable in par_tmfg/Makefile to the address of boost folder

## Input

The input to both implementations is a symmetric matrix. 
The format of the file is space separated numbers.
An example file is in the `dataset` folder.

## Datasets

The UCR data sets can be downloaded from [here](https://www.cs.ucr.edu/~eamonn/time_series_data_2018/).
The stock data can be obtained using the [Yahoo Finance API](https://pypi.org/project/yfinance/). Our data is obtained in Nov. 2021.

You can also download our data [here](https://console.cloud.google.com/storage/browser/par-filtered-graph-clustering).
There is a readme.md in the data repository linked above that explains how to use the datasets.

# Running Tests
For running time tests, we use `numactl`. It can be installed using `apt install numactl`. 

## To run HAC:

run `make` in `hac/general_hac`

PARLAY_NUM_THREADS=`wk` numactl -i all ./linkage `dataset` `n` `outpout` `method` `round`

* `wk` is the number of workers to use
* `numactl -i all` is optional 
* `dataset` is the file name of the input distance matrix
* `n` is the number of data points
* `output` is the file name of the output file for the resulting dendrogram
* `method` can be "comp" or "avg" for complete linkage and average linkage respectively
* `round` is the number of times to run the program


```bash
cd hac/general_hac
make
PARLAY_NUM_THREADS=${wk} numactl -i all ./linkage ../../datasets/iris-D.csv 150 outputs/iris-D_comp_dendro comp 1
```


##  To run parallel TMFG + DBHT:

run `make` in `par_tmfg`

PARLAY_NUM_THREADS=`wk` numactl -i all ./tmfg `S` `output` `n` `D` `method` `prefix` `round`

* `wk` is the number of workers to use
* `numactl -i all` is optional 
* `S` is the file name of the input similarity matrix
* `output` is the file name prefix of the output file for the resulting dendrogram (-Z) and the resulting TMFG (-P)
* `n` is the number of data points
* `D` is the file name of the input dissimilarity matrix. If D=0, will use D = sqrt(2(1-s))
* `method` can be "exact" or "prefix".
* `prefix` is the prefix size to insert in each round. it is ignored when method is exact
* `round` is the number of times to run the program

```bash
cd par_tmfg
make
PARLAY_NUM_THREADS=${wk} numactl -i all ./tmfg ../datasets/iris-R.csv outputs/iris-R 150 0 prefix 2 1
PARLAY_NUM_THREADS=1 ./tmfg ../datasets/iris-R.csv outputs/iris-R 150 0 exact 0 1
```

## To run Aste's implementation of TMFG/PMFG+DBHT:

UCR_PMFG(`dataset`, `inputdir`, `outputdir`)

UCR_TMFG(`dataset`, `inputdir`, `outputdir`)

* `dataset` is the name of the dataset.
* `inputdir` is the directory of the input dataset
* `outputdir` is the output directory

```bash
cd Aste
matlab -nojvm -nosplash -nodesktop  -nodisplay  -r  'UCR_PMFG("iris", "../datasets/", "outputs"); exit'  -logfile outputs/iris_pmfg_timing.txt
matlab -nojvm -nosplash -nodesktop  -nodisplay  -r  'UCR_TMFG("iris", "../datasets/", "outputs"); exit'  -logfile outputs/iris_tmfg_timing.txt
```

## To run KMeans:

```python
from sklearn.cluster import SpectralClustering
from sklearn.cluster import KMeans

KMeans(n_clusters=k, random_state=0).fit(X)

SpectralClustering(n_clusters=k, affinity="nearest_neighbors",
                                n_neighbors=n_neighbor,
                                assign_labels='discretize',
                                random_state=1, n_jobs=worker).fit(X)
```