#!/bin/bash

source ../sequential_benchs/scripts/common.sh

data_name=404

states=4

go="./main analyse_partitions $path_data/$data_name/$data_name.phy $path_data/$data_name/$data_name.part $states"

echo $go
$go

