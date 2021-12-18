#! /bin/bash
ARCH=$1
# Define variables
RUNS=10
MAXBLOCK=1024

echo "Tests starting for $ARCH ..." >> log.dat
echo >> log.dat

SIZE=256 
echo "Starting tests for size = $SIZE" >> log.dat
for((block = 16; block < $MAXBLOCK; block=2*block))
do
    for((k = 0; k < $RUNS; k++))
    do
    ./self_cuda_tile -s $SIZE -b $block >> log.dat
    done
    echo "Done Testing for blocksize $block" >> log.dat
    echo >> log.dat
done
echo "Done Testing for size $SIZE" >> log.dat
echo >> log.dat

SIZE=$(( $SIZE*2 )) #512
echo "Starting tests for size = $SIZE" >> log.dat
for((block = 16; block < $MAXBLOCK; block=2*block))
do
    for((k = 0; k < $RUNS; k++))
    do
    ./self_cuda_tile -s $SIZE -b $block >> log.dat
    done
    echo "Done Testing for blocksize $block" >> log.dat
    echo >> log.dat
done
echo "Done Testing for size $SIZE" >> log.dat
echo >> log.dat

SIZE=$(( $SIZE*2 )) #512
echo "Starting tests for size = $SIZE" >> log.dat

for((block = 16; block < $MAXBLOCK; block=2*block))
do
    for((k = 0; k < $RUNS; k++))
    do
    ./self_cuda_tile -s $SIZE -b $block >> log.dat
    done
    echo "Done Testing for blocksize $block" >> log.dat
    echo >> log.dat
done

echo "Done Testing for size $SIZE" >> log.dat
echo >> log.dat

SIZE=$(( $SIZE*2 )) #512
echo "Starting tests for size = $SIZE" >> log.dat

for((block = 16; block < $MAXBLOCK; block=2*block))
do
    for((k = 0; k < $RUNS; k++))
    do
    ./self_cuda_tile -s $SIZE -b $block >> log.dat
    done
    echo "Done Testing for blocksize $block" >> log.dat
    echo >> log.dat
done

echo "Done Testing for size $SIZE" >> log.dat
echo >> log.dat

SIZE=$(( $SIZE*2 )) #512
echo "Starting tests for size = $SIZE" >> log.dat

for((block = 16; block < $MAXBLOCK; block=2*block))
do
    for((k = 0; k < $RUNS; k++))
    do
    ./self_cuda_tile -s $SIZE -b $block >> log.dat
    done
    echo "Done Testing for blocksize $block" >> log.dat
    echo >> log.dat
done

echo "Done Testing for size $SIZE" >> log.dat
echo >> log.dat

SIZE=$(( $SIZE*2 )) #512
echo "Starting tests for size = $SIZE" >> log.dat

for((block = 16; block < $MAXBLOCK; block=2*block))
do
    for((k = 0; k < $RUNS; k++))
    do
    ./self_cuda_tile -s $SIZE -b $block >> log.dat
    done
    echo "Done Testing for blocksize $block" >> log.dat
    echo >> log.dat
done

echo "Done Testing for size $SIZE" >> log.dat
echo >> log.dat

SIZE=$(( $SIZE*2 )) #512
echo "Starting tests for size = $SIZE" >> log.dat

for((block = 16; block < $MAXBLOCK; block=2*block))
do
    for((k = 0; k < $RUNS; k++))
    do
    ./self_cuda_tile -s $SIZE -b $block >> log.dat
    done
    echo "Done Testing for blocksize $block" >> log.dat
    echo >> log.dat
done

echo "Done Testing for size $SIZE" >> log.dat
echo >> log.dat

SIZE=$(( $SIZE*2 )) #512
echo "Starting tests for size = $SIZE" >> log.dat

for((block = 16; block < $MAXBLOCK; block=2*block))
do
    for((k = 0; k < $RUNS; k++))
    do
    ./self_cuda_tile -s $SIZE -b $block >> log.dat
    done
    echo "Done Testing for blocksize $block" >> log.dat
    echo >> log.dat
done

echo "Done Testing for size $SIZE" >> log.dat
echo >> log.dat

SIZE=$(( $SIZE*2 )) #512
echo "Starting tests for size = $SIZE" >> log.dat

for((block = 16; block < $MAXBLOCK; block=2*block))
do
    for((k = 0; k < $RUNS; k++))
    do
    ./self_cuda_tile -s $SIZE -b $block >> log.dat
    done
    echo "Done Testing for blocksize $block" >> log.dat
    echo >> log.dat
done

echo "Done Testing for size $SIZE" >> log.dat
echo >> log.dat

SIZE=$(( $SIZE*2 )) #512
echo "Starting tests for size = $SIZE" >> log.dat

for((block = 16; block < $MAXBLOCK; block=2*block))
do
    for((k = 0; k < $RUNS; k++))
    do
    ./self_cuda_tile -s $SIZE -b $block >> log.dat
    done
    echo "Done Testing for blocksize $block" >> log.dat
    echo >> log.dat
done

echo "Done Testing for size $SIZE" >> log.dat
echo >> log.dat

SIZE=$(( $SIZE*2 )) #512
echo "Starting tests for size = $SIZE" >> log.dat

for((block = 16; block < $MAXBLOCK; block=2*block))
do
    for((k = 0; k < $RUNS; k++))
    do
    ./self_cuda_tile -s $SIZE -b $block >> log.dat
    done
    echo "Done Testing for blocksize $block" >> log.dat
    echo >> log.dat
done

echo "Done Testing for size $SIZE" >> log.dat
echo >> log.dat

echo "Finished Testing!" >> log.dat