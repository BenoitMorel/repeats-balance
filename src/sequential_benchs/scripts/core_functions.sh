#!/bin/bash

source scripts/common.sh


outputdir="${path_results}core_functions/"
iterations_likelihood=0
iterations_sumtable=0
iterations_derivatives=1000
dataset_number=4
srlookupsize=2000000

outputfile="core_functions_l${iterations_likelihood}_s${iterations_sumtable}_d${iterations_derivatives}"
outputfile="$outputdir$outputfile"
file=""
filebuffer=buffer.txt

echo $outputfile

write() {
  echo -n "$1"  >> $filebuffer
}
writeln() {
  echo  "$1"  >> $filebuffer
}

launch() {
  echo "./main core_functions $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10}"
  ./main core_functions $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} > temp
  cat temp | grep ms | tr -d '\n'   >> $filebuffer
  
}

bench_dataset() {
  writeln "\\hline"
  echo $8
  write "$8"
  write " & " 
  launch $path_data/59/unrooted.newick $path_data/59/59.phy 4  $1 $2 $3 $4 $5 $6 $7
  write " & " 
  launch $path_data/128/unrooted.newick $path_data/128/128.phy 4  $1 $2 $3 $4 $5 $6 $7
  write " & " 
  launch $path_data/404/unrooted.newick $path_data/404/404.phy 4  $1 $2 $3 $4 $5 $6 $7
  write " & " 
  launch $path_data/kyte/unrooted.newick $path_data/kyte/kyte.phy 4  $1 $2 $3 $4 $5 $6 $7
  writeln "\\\\"
}

bench_arch() {
  extension=".tex"
  file=$outputfile$1$extension
  mkdir -p $outputdir
  rm -r $filebuffer
  touch $filebuffer
  write "\\begin{tabular}{|l|"
  for i in $(seq 1 $dataset_number); do
    write "c|"
  done
  writeln "}"

  writeln "\\hline"
  writeln " & seq59 & seq128 & seq404  & kyte \\\\"

  cp $path_lib/libpll_repeats/* $path_lib/current  
  make clean 
  make 
  bench_dataset 0 0 $srlookupsize $iterations_likelihood $iterations_sumtable $iterations_derivatives $1 "tipinner" 
  bench_dataset 1 0 $srlookupsize $iterations_likelihood $iterations_sumtable $iterations_derivatives $1 "repeats" 

  writeln "\\hline"
  writeln "\\end{tabular}"
  mv $filebuffer $file
}

#bench_arch "cpu"
bench_arch "avx"
#bench_arch "sse"


