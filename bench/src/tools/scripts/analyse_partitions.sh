#!/bin/bash

source ../sequential_benchs/scripts/common.sh

data_name=pierre


states=4

go="./main analyse_partitions $path_data/$data_name/tax2.phy $path_data/$data_name/tax2.part $states"

echo $go
$go

