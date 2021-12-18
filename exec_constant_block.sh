#! /bin/bash
ARCH=$1
# Define variables
RUNS=10
block=128

echo "Tests starting for $ARCH ..." >> log.dat
echo >> log.dat

SIZE=256 

echo "Starting tests for size = $SIZE" >> log.dat
for((k = 0; k < $RUNS; k++))
do
    ./self_cuda_tile -s $SIZE -b $block >> log.dat
done
echo "Done Testing for size $SIZE" >> log.dat
echo >> log.dat

SIZE=$(( $SIZE*2 )) #512
echo "Starting tests for size = $SIZE" >> log.dat
for((k = 0; k < $RUNS; k++))
do
    ./self_cuda_tile -s $SIZE -b $block >> log.dat
done
echo "Done Testing for size $SIZE" >> log.dat
echo >> log.dat

SIZE=$(( $SIZE*2 )) # 1024
echo "Starting tests for size = $SIZE" >> log.dat
for((k = 0; k < $RUNS; k++))
do
    ./self_cuda_tile -s $SIZE -b $block >> log.dat
done
echo "Done Testing for size $SIZE" >> log.dat
echo >> log.dat

SIZE=$(( $SIZE*2 )) # 2048
echo "Starting tests for size = $SIZE" >> log.dat
for((k = 0; k < $RUNS; k++))
do
    ./self_cuda_tile -s $SIZE -b $block >> log.dat
done
echo "Done Testing for size $SIZE" >> log.dat
echo >> log.dat

SIZE=$(( $SIZE*2 )) # 4096
echo "Starting tests for size = $SIZE" >> log.dat
for((k = 0; k < $RUNS; k++))
do
    ./self_cuda_tile -s $SIZE -b $block >> log.dat
done
echo "Done Testing for size $SIZE" >> log.dat
echo >> log.dat

SIZE=$(( $SIZE*2 )) # 8192
echo "Starting tests for size = $SIZE" >> log.dat
for((k = 0; k < $RUNS; k++))
do
    ./self_cuda_tile -s $SIZE -b $block >> log.dat
done
echo "Done Testing for size $SIZE" >> log.dat
echo >> log.dat

SIZE=$(( $SIZE*2 )) # 16384
echo "Starting tests for size = $SIZE" >> log.dat
for((k = 0; k < $RUNS; k++))
do
    ./self_cuda_tile -s $SIZE -b $block >> log.dat
done
echo "Done Testing for size $SIZE" >> log.dat
echo >> log.dat

SIZE=$(( $SIZE*2 )) # 32768
echo "Starting tests for size = $SIZE" >> log.dat
for((k = 0; k < $RUNS; k++))
do
    ./self_cuda_tile -s $SIZE -b $block >> log.dat
done
echo "Done Testing for size $SIZE" >> log.dat
echo >> log.dat

echo "Finished Testing!" >> log.dat