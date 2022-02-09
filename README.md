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


## Installation Requirement

* g++ = 7.5.0 
* make
* [boost library](https://www.boost.org/)

After boost is installed, set the BOOST_ROOT variable in par_tmfg/Makefile to the address of boost folder

## Input

The input to both implementations is a symmetric matrix. 
The format of the file is space separated numbers.
An example file is in the `dataset` folder.


## Running Tests
For running time tests, we use `numactl`. It can be installed using `apt install numactl`. 

To run HAC:

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
PARLAY_NUM_THREADS=${wk} numactl -i all ./linkage iris-D 150 iris_dendro comp 1
```


To run parallel TMFG + DBHT:

PARLAY_NUM_THREADS=`wk` numactl -i all ./tmfg `S` `output` `n` `D` `method` `prefix` `round`

* `wk` is the number of workers to use
* `numactl -i all` is optional 
* `S` is the file name of the input similarity matrix
* `output` is the file name of the output file for the resulting dendrogram (-Z) and the resulting TMFG (-P)
* `n` is the number of data points
* `D` is the file name of the input dissimilarity matrix. If D=0, will use D = sqrt(2(1-s))
* `method` can be "exact" or "prefix".
* `prefix` is the prefix size to insert in each round. it is ignored when method is exact
* `round` is the number of times to run the program

```bash
cd par_tmfg
make
PARLAY_NUM_THREADS=${wk} numactl -i all ./tmfg iris-R iris 150 0 exact 0 1
```

## datasets

The UCR data sets can be downloaded from [here](https://www.cs.ucr.edu/~eamonn/time_series_data_2018/)
The stock data can be obtained using the [Yahoo Finance API](https://pypi.org/project/yfinance/). Our data is obtained in Nov. 2021.

