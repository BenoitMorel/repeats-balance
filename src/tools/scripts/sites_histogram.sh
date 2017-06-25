#!/bin/bash

source ../sequential_benchs/scripts/common.sh

data_name=404


go=" ./main sites_histogram $path_data/$data_name/$data_name.phy"

echo $go
$go 


