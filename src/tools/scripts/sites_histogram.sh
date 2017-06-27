#!/bin/bash

source ../sequential_benchs/scripts/common.sh

data_name=kyte
lim=0
data_lim=${data_name}_lim$lim
output=/home/morelbt/github/raxml-ng-alexey/bin/data/${data_lim}/${data_lim}
output=${data_lim}
mkdir -p /home/morelbt/github/raxml-ng-alexey/bin/data/${data_lim}
go=" ./main sites_histogram $path_data/$data_name/$data_name.phy $path_data/$data_name/$data_name.part $output $lim"

echo $go
$go 


