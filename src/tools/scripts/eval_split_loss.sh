#!/bin/bash

source ../sequential_benchs/scripts/common.sh

data_name=404
#data_name=1kite_2013_10randomtaxa
#empty string for a random tree
states=4
splits=2

go=" ./main eval_split_loss $path_data/$data_name/$data_name.phy $states $splits"

echo $go
$go 
