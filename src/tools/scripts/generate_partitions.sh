#!/bin/bash

source ../sequential_benchs/scripts/common.sh

data_name=1kite_matc_red
model=DNA
partitions=200

go="./main generate_partitions $path_data/${data_name}/$data_name.phy $partitions $model $path_data/$data_name/${data_name}_$partitions.part"
#$go="./main generate_partitions ../../../../raxml-ng/bin/data/1kyte_hyme/1kyte_hyme.phy $partitions $model ../../../../raxml-ng/bin/data/1kyte_hyme/1kyte_hyme_$partitions.part"

echo $go
$go

