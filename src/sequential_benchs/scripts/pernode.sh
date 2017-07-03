#!/bin/bash

source scripts/common.sh

iterations=20

arch[0]="avx"
#arch[1]="cpu"
#arch[2]="sse"

srlookupsize=2000000


#outputtex="$path_results/pernode/404avx_$iterations.tex"
#dataset[0]="$path_data/404/unrooted.newick $path_data/404/404.phy 4"
outputtex="$path_results/pernode/kyte_$iterations.tex"
dataset[0]="$path_data/kyte/unrooted.newick $path_data/kyte/kyte.phy 4"



runs[0]="$srlookupsize 0 0 0 $iterations"
libs[0]="$path_lib/libpll_repeats/"
runname[0]="tippat"

runs[1]="$srlookupsize 0 1 1 $iterations"
libs[1]="$path_lib/libpll_repeats/"
runname[1]="repeats"


#runs[2]="0 0 0 $srlookupsize $iterations"
#libs[2]="$path_lib/libpll_benoit_dev"
#runname[2]="tip pat"

cp ${libs[0]} $path_lib/current/
make clean
make

echo "" > starttab
echo "" > endtab
echo "" > $outputtex

echo -n "\\begin{tabular}{|l|" >> starttab
for i in "${!runs[@]}"; do
  echo -n "c|" >> starttab
done 
echo -n "c|c|c|" >> starttab
echo "}" >> starttab
echo -n "\\hline node & psrsize & lsrsize & rsrsize  " >> starttab
for i in "${!runs[@]}"; do
  echo -n " & "${runname[$i]}  >> starttab
done
echo "\\\\" >> starttab

echo "\\hline" >> endtab
echo "\\end{tabular} \\" >> endtab
echo "" >> endtab
echo "---" >> endtab
echo "" >> endtab

echo "" > $outputtex
# for all arch
for a in "${!arch[@]}"; do
  cat starttab >> $outputtex
  temps=""
  temp="tempinit"
  temps=$temps" "$temp
  echo -n "" > $temp
  # print nodes info (number of repeats of children etc.)
  mycommand="./main pernode ${dataset[0]} avx -1 $srlookupsize 0"
  echo $mycommand 
  $mycommand > $temp
  for i in "${!runs[@]}"; do
    cp ${libs[$i]}/* $path_lib/current/
    make clean
    make
    temp="temp"$i
    temps=$temps" "$temp
    echo -n "" > $temp
    mycommand="./main pernode ${dataset[0]} ${arch[$a]} -1 ${runs[$i]}"
    echo $mycommand 
    $mycommand > $temp
  done
  paste -d" " $temps | sed 's/ / \& /g' | sed 's/^/    \\hline node/g' | sed 's/$/\\\\/g'  > temp
  lines=0
  while IFS= read -r line
  do
    printf '%s\n' "${line}" >> $outputtex
    lines=$(($lines + 1))
    if [[ $((lines % 20)) -eq 0 ]]; 
    then
      cat endtab >> $outputtex
      cat starttab >> $outputtex
    fi
  done < temp
  
  rm temp*
  
  
  cat endtab >> $outputtex
done

rm starttab
rm endtab
