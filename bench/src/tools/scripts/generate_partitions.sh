#!/bin/bash

source ../sequential_benchs/scripts/common.sh

data_name=kyte
model=DNA
partitions=5


go="./main generate_partitions $path_data/${data_name}/$data_name.phy $partitions $model $path_data/$data_name/$data_name.part"

echo $go
$go

