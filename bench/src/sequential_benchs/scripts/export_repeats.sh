#!/bin/bash

source scripts/common.sh

#./main export_repeats $path_data/404/unrooted.newick $path_data/404/404.phy $path_data/404/404.part plip
gdb --args  ./main export_repeats $path_data/kyte/unrooted.newick $path_data/kyte/kyte.phy $path_data/kyte/kyte.part /home/benoit/github/phylogenetic_hypergraphs/data/DNA_kyte
#./main export_repeats $path_data/404/unrooted.newick $path_data/404/split-partitions/ 11 plop


