#!/bin/bash

source ../sequential_benchs/scripts/common.sh

data=generated_terrace_1
states=4
number=10

go="./main random_trees_likelihoods $path_data/$data/$data.phy $path_data/$data/$data.part $states $number"

echo $go
$go 



