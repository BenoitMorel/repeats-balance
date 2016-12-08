#!/bin/bash
outputdir="../../../results/pernode_benchs/"
tempfile1="temp1"
tempfile2="temp2"
iterations=500
size=10
pathnumber=20

export LD_LIBRARY_PATH=../../lib/current/

make
echo "with sites repeats"
./randompath_bench/randompath_bench ../../../data/404/unrooted.newick ../../../data/404/404.phy $iterations 1 0 $size $pathnumber
echo "no sites repeats"
./randompath_bench/randompath_bench ../../../data/404/unrooted.newick ../../../data/404/404.phy $iterations 0 0 $size $pathnumber



