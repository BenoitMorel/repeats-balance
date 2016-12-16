#!/bin/bash

export LD_LIBRARY_PATH=../../lib/current/

make
echo "with sites repeats"
./rootpath_bench/rootpath_bench ../../../data/404/unrooted.newick ../../../data/404/404.phy 100 1 1 0.94 10
echo "with tip pattern"
./rootpath_bench/rootpath_bench ../../../data/404/unrooted.newick ../../../data/404/404.phy 100 0 0 0.94 10



