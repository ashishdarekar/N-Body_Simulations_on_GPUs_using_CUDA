#! /bin/bash

RUNS=20
block=256

echo "Tests starting for ..." >> log.dat
echo >> log.dat

SIZE=2048

echo "Starting tests for size = $SIZE" >> log.dat
for((k = 0; k < $RUNS; k++))
do
    ./self_cuda_tile -s $SIZE -b $block >> log.dat
done
echo "Done Testing for size $SIZE" >> log.dat
echo >> log.dat

echo "Finished Testing!" >> log.dat