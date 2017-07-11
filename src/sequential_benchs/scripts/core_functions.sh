#!/bin/bash

source scripts/common.sh


outputdir="${path_results}core_functions/"
iterations_partials=25
iterations_likelihood=0
iterations_sumtable=0
iterations_derivatives=0

datasets=(
          #"$path_data/59/unrooted.newick $path_data/59/59.phy 4 4" 
          #"$path_data/128/unrooted.newick $path_data/128/128.phy 4 4"
          #"$path_data/404/unrooted.newick $path_data/404/404.phy 4 4"
          "$path_data/kyte/unrooted.newick $path_data/kyte/kyte.phy 4 1"
          "$path_data/kyte/unrooted.newick $path_data/kyte/kyte.phy 4 2"
          "$path_data/kyte/unrooted.newick $path_data/kyte/kyte.phy 4 4"
          "$path_data/kyte/unrooted.newick $path_data/kyte/kyte.phy 4 8"
          )

datasets_names=(
                #"59"
                #"128"
                #"404"
                "kyte (1 cat)"
                "kyte (2 cat)"
                "kyte (4 cat)"
                "kyte (8 cat)"
                )

outputfile="rates_core_functions_p${iterations_partials}_l${iterations_likelihood}_s${iterations_sumtable}_d${iterations_derivatives}"
outputfile="$outputdir$outputfile"

dataset_number=${#datasets[@]}

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
  echo "./main core_functions $@"
  ./main core_functions $@ > temp
  cat temp | grep ms | tr -d '\n'   >> $filebuffer
  
}

bench_dataset() {
  writeln "\\hline"
  echo ${10}
  write "${10}"
  for ds in "${datasets[@]}"
  do
    write " & " 
    launch $ds  $1 $2 $3 $4 $5 $6 $7 $8 $9
  done
  writeln "\\\\"
}

bench_arch() {
  extension=".tex"
  file=$outputfile$1$extension
  mkdir -p $outputdir
  rm -r $filebuffer
  touch $filebuffer
  echo "output in $file"
  write "\\begin{tabular}{|l|"
  for i in $(seq 1 $dataset_number); do
    write "c|"
  done
  writeln "}"

  writeln "\\hline"
  for dsn in "${datasets_names[@]}"
  do
    write " & $dsn"
  done
  writeln "\\\\"

  cp $path_lib/libpll_repeats/* $path_lib/current  
  make clean 
  make 
  bench_dataset 0 0 0 $srlookupsize $iterations_partials $iterations_likelihood $iterations_sumtable $iterations_derivatives $1 "tipinner" 
  bench_dataset 1 1 0 $srlookupsize $iterations_partials $iterations_likelihood $iterations_sumtable $iterations_derivatives $1 "repeats" 
  bench_dataset 1 0 0 $srlookupsize $iterations_partials $iterations_likelihood $iterations_sumtable $iterations_derivatives $1 "repeats no update" 


  writeln "\\hline"
  writeln "\\end{tabular}"
  mv $filebuffer $file
}

#bench_arch "cpu"
bench_arch "avx"
#bench_arch "sse"


