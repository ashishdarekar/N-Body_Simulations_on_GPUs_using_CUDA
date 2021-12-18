CC=g++

CUDA_INSTALL_DIR=/usr/local/cuda

NVCC=nvcc
CFLAGS=-O3 

STD = -std=c++11

ARCHFLAG = -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35 -gencode arch=compute_37,code=sm_37 -gencode arch=compute_50,code=sm_50 -gencode arch=compute_52,code=sm_52 -gencode arch=compute_60,code=sm_60 -gencode arch=compute_61,code=sm_61 -gencode arch=compute_70,code=sm_70 -gencode arch=compute_75,code=sm_75

OPENMP = -fopenmp 

INCLUDES = -I. -I$(CUDA_INSTALL_DIR)/include

ODIR=obj
LDIR =lib

LIBS=-lm #-fopenmp

all:build
build: self_cuda_tile self_cuda_tile_test serial serial_omp

#self_cuda_tile: self_cuda_tile.cu
self_cuda_tile: self_cuda_tile.cu
	$(NVCC) $(STD) $(ARCHFLAG) -o $@ $^ $(CFLAGS) $(LIBS)

#self_cuda_tile_test: self_cuda_tile_test.cu
self_cuda_tile_test: self_cuda_tile_test.cu
	$(NVCC) $(STD) $(ARCHFLAG) -o $@ $^ $(CFLAGS) $(LIBS)

serial: self_serial.cpp
	$(CC) $(STD) $(LIBS) $(INCLUDES) -o $@ $^ $(CFLAGS)

serial_omp: self_serial.cpp
	$(CC) $(STD) $(LIBS) $(OPENMP) $(INCLUDES) -o $@ $^ $(CFLAGS)

test: build
	./self_cuda_tile_test -s 1024 -b 128

clean:
	rm -rf self_cuda_tile self_cuda_tile_test cuda_solution_cuda_seminar serial serial_omp log.dat serial_solution_cuda_seminar




