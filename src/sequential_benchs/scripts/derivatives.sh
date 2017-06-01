#!/bin/bash

source scripts/common.sh


outputdir="$path_results/latex/derivatives/"
iterations=10000
dataset_number=3
srlookupsize=2000000

outputfile="bench_derivatives_new_${iterations}_"
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
  echo "./main derivatives $1 $2 $3 $4 $5 $6 $7 $8"
  ./main derivatives $1 $2 $3 $4 $5 $6 $7 $8 > temp
  cat temp | grep ms | tr -d '\n'   >> $filebuffer
  
}

bench_dataset() {
  writeln "\\hline"
  echo $6
  write "$6"
  write " & " 
  launch $path_data/59/unrooted.newick $path_data/59/59.phy 4  $1 $2 $3 $4 $5
  write " & " 
  launch $path_data/128/unrooted.newick $path_data/128/128.phy 4  $1 $2 $3 $4 $5
  write " & " 
  launch $path_data/404/unrooted.newick $path_data/404/404.phy 4  $1 $2 $3 $4 $5
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
  writeln " & seq59 & seq128 & seq404  \\\\"

  cp $path_lib/libpll_repeats/* $path_lib/current  
  make clean 
  make 
  bench_dataset 0 0 $srlookupsize $iterations $1 "tipinner benoit" 
  bench_dataset 1 0 $srlookupsize $iterations $1 "repeats benoit" 

  writeln "\\hline"
  writeln "\\end{tabular}"
  mv $filebuffer $file
}

#bench_arch "cpu"
bench_arch "avx"
#bench_arch "sse"


