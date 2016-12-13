#!/bin/bash
outputdir="../../../results/pernode_benchs/"
tempfile1="temp1"
tempfile2="temp2"
iterations=200
runs=1
runs2=15

export LD_LIBRARY_PATH=../../lib/current/

make
for ((seed=1;$seed < $runs2;seed=$seed + 1))
do
  ./randompath_bench/randompath_bench ../../../data/404/unrooted.newick ../../../data/404/404.phy $runs $iterations 1 0 95 $seed
./randompath_bench/randompath_bench ../../../data/404/unrooted.newick ../../../data/404/404.phy $runs $iterations 0 0 95 $seed
echo ""
done


