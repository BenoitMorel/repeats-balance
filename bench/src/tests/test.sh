#!/bin/bash
path_lib=../../lib
path_common=../common
export LD_LIBRARY_PATH=$path_lib/current:$path_lib/common
make
./main
