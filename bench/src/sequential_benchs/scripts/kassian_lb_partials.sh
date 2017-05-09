#!/bin/bash

source scripts/common.sh

data_name=404

states=4
lookupsize=0
cpus=8
iterations=20
use_repeats=1
update_repeats=0
debugger=
#debugger=valgrind
#debugger="gdb --args "

go="./main kassian_lb_partials $path_data/$data_name/unrooted.newick $path_data/$data_name/$data_name.phy $path_data/$data_name/$data_name.part $states $use_repeats $update_repeats $lookupsize $iterations $cpus"

echo $debugger $go
$debugger $go
