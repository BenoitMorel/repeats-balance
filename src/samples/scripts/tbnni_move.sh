#!/bin/bash

source ../sequential_benchs/scripts/common.sh

data_name=tbnni_simple
states=4

go="./main tbnni_move $path_data/$data_name/$data_name.phy $path_data/$data_name/$data_name.part  $states $path_data/$data_name/start.newick"

echo $go
$go 


