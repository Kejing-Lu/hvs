# Project Title

HVS: Hierarchical Graph Structure Based on Voronoi Diagrams for Solving Approximate Nearest Neighbor Search

## Getting Started

### Prerequisites

* OpenCV 3.30+
* GCC 4.9+ with OpenMP
* CMake 2.8+
* Boost 1.55+
* TCMalloc

### Datasets and query sets

* Deep (https://yadi.sk/d/11eDCm7Dsn9GA)
* Other real datasets (https://www.cse.cuhk.edu.hk/systems/hash/gqr/datasets.html)

### Compile On Ubuntu

Complie HNSW

```shell
$ cd hnsw/
$ mkdir build/ && cd build/
cmake ..
make 
```

Compile HVS (based on NSG)

```shell
$ cd nsg-master/
$ mkdir build/ && cd build/
$ cmake ..
$ make
```

## Commands

* `K` is the value of top-K
* NSG parameters `L`, `R` and `C` were set to 50, 60 ,500 in our experiments
* `T` and `delta` are user-specified parameters of HVS

```shell
$ ./hnsw/build/main ${dataset} ${datasize} ${dimension} ${T} ${delta}
$ ./nsg-master/build/tests/test_nsg_index ${dataset} data_graph.fvecs ${L} ${R} ${C} data_index.fvecs ${T}
$ ./nsg-master/build/tests/test_nsg_optimized ${dataset} ${queryset} data_index.fvecs ${T} ${K} ${groundtruthset}
```
