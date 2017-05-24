#!/bin/bash

source scripts/common.sh

data_name=antl_1_1_nt2

states=4
lookupsize=0
cpus=128
iterations=2
use_repeats=1
update_repeats=0
debugger= 
randomized=1
#debugger=valgrind
#debugger="gdb --args "

go="valgrind  ./main kassian_lb_partials $path_data/$data_name/unrooted.newick $path_data/$data_name/$data_name.phy $path_data/$data_name/$data_name.part $states $use_repeats $update_repeats $lookupsize $iterations $cpus $randomized"

echo $debugger $go
$debugger $go
