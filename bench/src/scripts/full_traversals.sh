

#!/bin/bash
outputdir="../../results/sequential_benchs/"
outputfile="bench_tipinner_500"
outputfile="$outputdir$outputfile"
iterations=500
dataset_number=3
file=""
filebuffer=buffer.txt
export LD_LIBRARY_PATH=../lib/current

echo $outputfile

write() {
  echo -n "$1"  >> $filebuffer
}
writeln() {
  echo  "$1"  >> $filebuffer
}

launch() {
  echo "./main full_traversal $1 $2 $3 $4 $5 $6 $7"
  ./main full_traversal $1 $2 $3 $4 $5 $6 $7  > temp
  cat temp | grep ms | tr -d '\n'   >> $filebuffer
  
}

bench_dataset() {
  writeln "\\hline"
  echo $6
  write "$6"
  write " & " 
  launch ../../data/59/unrooted.newick ../../data/59/59.phy  $1 $2 $3 $4 $5 
  write " & " 
  launch ../../data/128/unrooted.newick ../../data/128/128.phy  $1 $2 $3 $4 $5
  write " & " 
  launch ../../data/404/unrooted.newick ../../data/404/404.phy  $1 $2 $3 $4 $5
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

  cp ../lib/libpll_benoit_dev/* ../lib/current && make clean && make # header changed
  bench_dataset 1 1 1000000 $iterations $1 "repeats 1000000" 
  cp ../lib/libpll_benoit_tipinner/* ../lib/current && make clean && make # header changed
  bench_dataset 1 1 1000000 $iterations $1 "bmorel sites repeats (ti opt) 1000000" 

  writeln "\\hline"
  writeln "\\end{tabular}"
  mv $filebuffer $file
}

bench_arch "cpu"
#bench_arch "avx"
#bench_arch "sse"


