#!/bin/bash

source ../sequential_benchs/scripts/common.sh

data_name=1kyte_hyme
#empty string for a random tree
tree=
#tree=$path_data/$data_name/unrooted.newick
states=4

partition_suffix=

go="./main analyse_partitions $path_data/$data_name/$data_name.phy $path_data/$data_name/$data_name$partition_suffix.part $states $tree"

echo $go
$go | grep normalized | awk '{print "("NR"," $3")" }'
#$go | grep plop
#> ${data_name}$partition_suffix


