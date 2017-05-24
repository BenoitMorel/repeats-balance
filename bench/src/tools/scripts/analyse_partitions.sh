#!/bin/bash

source ../sequential_benchs/scripts/common.sh

data_name=1kyte_hyme
partition_suffix=_1000


states=4

go="./main analyse_partitions $path_data/$data_name/$data_name.phy $path_data/$data_name/$data_name$partition_suffix.part $states"

echo $go
$go | grep unique | awk '{print "("NR"," $4")" }'
#> ${data_name}$partition_suffix


