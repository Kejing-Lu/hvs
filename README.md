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

* Tiny80M (https://drive.google.com/file/d/1PzW9cqi8VbzH9Bu_7UDWAXjA2DlBorG3/view?usp=sharing)
* Other real datasets (https://www.cse.cuhk.edu.hk/systems/hash/gqr/datasets.html)

### Compile On Ubuntu

Complie HVS (based on NSW)

```shell
$ cd hnsw/
$ mkdir build/ && cd build/
cmake ..
make 
```
## Commands

* `K` is the value of top-K, `L` is the value of efsearch and `qn` is the size of query set
* `T` and `delta` are user-specified parameters of HVS
* The data set, query set and the ground_truth set are stored in ${dPath}.ds, ${qPath}.q and truth.gt

Build HVS index
```shell
./hnsw/build/main ${dPath}.ds ${qPath}.q ${n} ${d} ${T} ${qn} ${K} ${L}
```
Search in HVS
```shell
./hnsw/build/main ${dPath}.ds nullptr ${n} ${d} ${T} -1 ${delta} -1
```

## A running example (ImageNet)
* Donwload the dataset, query set and the ground truth set of ImageNet from the following link
https://drive.google.com/file/d/1WV78sZT1j1oz-GoZTQdyKaNQgulRLe2O/view?usp=sharing
* Put all three files in the main folder
* Run the script
```shell
bash run_image.sh
```
