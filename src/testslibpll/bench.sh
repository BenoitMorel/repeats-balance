#!/bin/bash
outputdir="../../results/sequential_benchs/"
outputfile="bench_500_iterations"
outputfile="$outputdir$outputfile"
iterations=500
dataset_number=3
file=""
export LD_LIBRARY_PATH=.

echo $outputfile

write() {
  echo -n "$1"  >> $file
}
writeln() {
  echo  "$1"  >> $file
}

launch() {
  echo "./test $1 $2 $3 $4 $5 $6 $7"
  ./test $1 $2 $3 $4 $5 $6 $7 > temp
  cat temp | grep ms | tr -d '\n'   >> $file  
  cat temp | grep ll
  
}

bench_dataset() {
  writeln "\\hline"
  write "$5"

  write " & " 
  launch ../../data/59/treechar_samples/tree59. ../../results/exports/59_10  $1 $2 $3 $4 $6
  write " & " 
  launch ../../data/128/treechar_samples/tree128. ../../results/exports/128_10  $1 $2 $3 $4 $6
  write " & " 
  launch ../../data/404/treechar_samples/tree404. ../../results/exports/404_10  $1 $2 $3 $4 $6
  #write " & " 
  #if [ $7 = "1" ]
  #then
  #  launch ../../data/94/tree94. ../../results/exports/94_10  $1 $2 $3 $4 $6
  #else
  #  echo "skip"
  #  write " skip "
  #fi
  writeln "\\\\"
}

bench_arch() {
  extension=".tex"
  file=$outputfile$1$extension
  mkdir -p $outputdir
  rm $file
  touch $file
  write "\\begin{tabular}{|l|"
  for i in $(seq 1 $dataset_number); do
    write "c|"
  done
  writeln "}"

  writeln "\\hline"
  writeln " & seq59 & seq128 & seq404  \\\\"

 # cp tomas/* .
 # make clean && make # header changed
 # bench_dataset $iterations 1 0 0 "xflouris default mode" $1 1
 # cp benoit1000/* .
 # make clean && make # header changed
 # bench_dataset $iterations 1 0 0 "bmorel default mode" $1 1
  
  cp tomas/* .
  make clean && make # header changed
  bench_dataset $iterations 1 0 1 "xflouris tip pattern" $1 1
  cp benoit1000/* .
  make clean && make # header changed
  bench_dataset $iterations 1 0 1 "bmorel tip pattern" $1 1
  
  
  cp benoit1000/* .
  bench_dataset $iterations 1 1 0 "bmorel sites repeats M=1000" $1 0
  cp benoit10000/* .
  bench_dataset $iterations 1 1 0 "bmorel sites repeats M=10000" $1 0
  cp benoit1000000/* . 
  bench_dataset $iterations 1 1 0 "bmorel sites repeats M=1000000" $1 1
  cp benoit1000000/* . 
  bench_dataset $iterations 1 2 0 "bmorel sites repeats no update sr  M=1000000" $1 1
  writeln "\\hline"
  writeln "\\end{tabular}"
}

bench_arch "avx"
bench_arch "cpu"
bench_arch "sse"

