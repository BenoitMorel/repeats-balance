#!/bin/bash

source ../sequential_benchs/scripts/common.sh

data=generated_50_1000
name=../../../data/$data/$data

go="./main generate_repeats_dataset $name.phy $name.part $name.newick "


echo $go
$go

