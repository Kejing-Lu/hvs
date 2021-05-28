#!/bin/bash
#make
#rm *.o

# ------------------------------------------------------------------------------
#  Parameters
# ------------------------------------------------------------------------------
dname=image
n=2340373
d=150
qn=200
T=1
delta=0.5
K=100
L=1000
dPath=./${dname}
qPath=./${dname}


./hnsw/build/main ${dPath}.ds nullptr ${n} ${d} ${T} -1 ${delta} -1

./hnsw/build/main ${dPath}.ds ${qPath}.q ${n} ${d} ${T} ${qn} ${K} ${L}
