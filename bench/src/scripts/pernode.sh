#!/bin/bash

outputtex="../../results/pernode/pernode.tex"

iterations=2000
arch[0]="cpu"
arch[1]="avx"
arch[2]="sse"

dataset[0]="../../data/59/unrooted.newick ../../data/59/59.phy"
runs[0]="1 1  1000000 $iterations"
libs[0]="../lib/libpll_benoit_dev"
runname[0]="repeats"
runs[1]="1 1  1000000 $iterations"
libs[1]="../lib/libpll_benoit_tipinner"
runname[1]="repeats tipinner"
runs[2]="0 0  1000000 $iterations"
libs[2]="../lib/libpll_benoit_dev"
runname[2]="tip pat"

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
echo -n "\\hline node & depth & lsize & rsize " >> starttab
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
for a in "${!arch[@]}"; do
  cat starttab >> $outputtex
  temps=""
  temp="tempinit"
  temps=$temps" "$temp
  mycommand="./main pernode ${dataset[0]}"
  echo $mycommand 
  $mycommand > $temp
  for i in "${!runs[@]}"; do
    cp ${libs[$i]}/* ../lib/current/
    make clean
    make
    temp="temp"$i
    temps=$temps" "$temp
    mycommand="./main pernode ${dataset[0]} ${runs[$i]} ${arch[$a]}"
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
  
  
  cat endtab >> $outputtex
done

rm starttab
rm endtab
rm $temps

#./main pernode ${dataset[0]} ${runs[0]}
#./main pernode ${dataset[0]} ${runs[1]}



