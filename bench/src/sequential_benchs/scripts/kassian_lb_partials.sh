#!/bin/bash

source scripts/common.sh

data_name=kyte

 ./main kassian_lb_partials $path_data/$data_name/unrooted.newick $path_data/$data_name/$data_name.phy $path_data/$data_name/$data_name.part 4 1 1 0 50 8 


