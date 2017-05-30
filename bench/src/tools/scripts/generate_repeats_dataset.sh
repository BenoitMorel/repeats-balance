#!/bin/bash

source ../sequential_benchs/scripts/common.sh

name=../../../data/generated/generated

go="./main generate_repeats_dataset $name.phy $name.part $name.newick "
#$go="./main generate_partitions ../../../../raxml-ng/bin/data/1kyte_hyme/1kyte_hyme.phy $partitions $model ../../../../raxml-ng/bin/data/1kyte_hyme/1kyte_hyme_$partitions.part"

echo $go
$go

