export OMP_NUM_THREADS=4
icpc -O3 -xHost -fopenmp step-4.cpp
time bash ./4-part.sh
icpc -O3 -xHost step-3.cpp
time bash ./4-part.sh
