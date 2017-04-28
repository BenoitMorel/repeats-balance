#!/bin/bash

source scripts/common.sh

outputdir="$path_results/sequential_benchs/"
iterations=500
outputfile="repeats_${iterations}_"
#outputfile="plop_tipinner_hitspc_${iterations}_"
outputfile="$outputdir$outputfile"
dataset_number=4
srlookupsize=2000000
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
  echo "./main full_traversal $1 $2 $3 $4 $5 $6 $7 $8"
  ./main full_traversal $1 $2 $3 $4 $5 $6 $7 $8 > temp
  cat temp | grep ms | tr -d '\n'   >> $filebuffer
  
}

plop_dataset() {
  writeln "\\hline"
  echo $6
  write "$6"
  write " & " 
  launch $path_data/59/unrooted.newick $path_data/59/59.phy 4 $1 $2 $3 $4 $5 
  write " & " 
  launch $path_data/128/unrooted.newick $path_data/128/128.phy 4  $1 $2 $3 $4 $5 
  write " & " 
  launch $path_data/404/unrooted.newick $path_data/404/404.phy 4 $1 $2 $3 $4 $5 
  write " & " 
  launch $path_data/140/unrooted.newick $path_data/140/140_big.phy 20 $1 $2 $3 $4 $5
  writeln "\\\\"
}

plop_arch() {
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
  writeln " & seq59 & seq128 & seq404 & seq140 \\\\"

  cp $path_lib/libpll_benoit_repeats/* $path_lib/current  
  make clean 
  make 
  plop_dataset 0 0 $srlookupsize $iterations $1 "tip pattern" 
  plop_dataset 1 0 $srlookupsize $iterations $1 "repeats" 
  
  
  writeln "\\hline"
  writeln "\\end{tabular}"
  mv $filebuffer $file
  echo "result in $file"
}

#plop_arch "cpu"
plop_arch "avx"
#plop_arch "sse"


