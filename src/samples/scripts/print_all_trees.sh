#!/bin/bash

source ../sequential_benchs/scripts/common.sh

data_name=minimal-5
states=4
number=10

go="./main print_all_trees $path_data/$data_name/$data_name.phy $states"

echo $go
$go 


