

#!/bin/bash
outputdir="../../results/update_repeats/"
iterations=10
outputfile="bench_repeats_${iterations}_"
#outputfile="bench_tipinner_hitspc_${iterations}_"
outputfile="$outputdir$outputfile"
dataset_number=3
srlookupsize=2000000
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
  echo "./main full_traversal $1 $2 $3 $4 $5 $6 $7 $8"
  ./main update_repeats $1 $2 $3 $4 $5 $6 $7 $8 > temp
  cat temp | grep ms | tr -d '\n'   >> $filebuffer
  
}

bench_dataset() {
  writeln "\\hline"
  echo $7
  write "$7"
  write " & " 
  launch ../../data/59/unrooted.newick ../../data/59/59.phy  $1 $2 $3 $4 $5 $6
  write " & " 
  launch ../../data/128/unrooted.newick ../../data/128/128.phy  $1 $2 $3 $4 $5 $6
  write " & " 
  launch ../../data/404/unrooted.newick ../../data/404/404.phy  $1 $2 $3 $4 $5 $6
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

  cp ../lib/libpll_benoit_dev/* ../lib/current  
  make clean 
  make 
  bench_dataset 1 1 0 $srlookupsize $iterations $1 "repeats benoit dev" 
  cp ../lib/libpll_benoit_repeats_integration/* ../lib/current  
  make clean 
  make 
  bench_dataset 1 1 0 $srlookupsize $iterations $1 "ri inte" 
  cp ../lib/libpll_benoit_ri_optims/* ../lib/current  
  make clean 
  make 
  bench_dataset 1 1 0 $srlookupsize $iterations $1 "ri optims" 
  cp ../lib/libpll_benoit_ri_optims2/* ../lib/current  
  make clean 
  make 
  bench_dataset 1 1 0 $srlookupsize $iterations $1 "ri optims2" 
  cp ../lib/libpll_benoit_ri_optims2/* ../lib/current  
  make clean 
  make 
  bench_dataset 1 1 0 $srlookupsize $iterations $1 "ri optims3" 
  #cp ../lib/libpll_benoit_tipinner/* ../lib/current 
  #make clean 
  #make 
  #bench_dataset 1 1 0 $srlookupsize $iterations $1 "repeats + bclv 2000000" 
  #bench_dataset 1 1 2048 $srlookupsize $iterations $1 "repeats + bclv2 2000000" 

  writeln "\\hline"
  writeln "\\end{tabular}"
  mv $filebuffer $file
}

#bench_arch "cpu"
bench_arch "avx"
#bench_arch "sse"

