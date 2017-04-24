#!/bin/bash


iterations=10000
arch[0]="cpu"
#arch[1]="avx"
#arch[2]="sse"

#nodes[0]=54
nodes[0]=132

srlookupsize=2000000

#outputtex="$path_results/singlenodes/59.tex"
#dataset[0]="$path_data/59/unrooted.newick $path_data/59/59.phy"

outputtex="$path_results/singlenodes/140.tex"
dataset[0]="$path_data/140/unrooted.newick $path_data/140/140.phy"

#outputtex="$path_results/singlenodes/404.tex"
#dataset[0]="$path_data/404/unrooted.newick $path_data/404/404.phy"


runs[0]="1 1  $srlookupsize $iterations"
libs[0]="$path_lib/libpll_benoit_dev"
runname[0]="repeats"

runs[1]="1 0  $srlookupsize $iterations"
libs[1]="$path_lib/libpll_benoit_tipinner"
runname[1]="tipinner opt"

#runs[2]="0 0  $srlookupsize $iterations"
#libs[2]="$path_lib/libpll_benoit_dev"
#runname[2]="tip pattern"

cp $path_lib/libpll_benoit_dev/* $path_lib/current/
make clean 
make

echo "" > starttab
echo "" > endtab
echo "" > $outputtex

echo -n "\\begin{tabular}{|l|" >> starttab
for i in "${!runs[@]}"; do
  echo -n "c|" >> starttab
done 
echo -n "c|c|c|c|c|c|" >> starttab
echo "}" >> starttab
echo -n "\\hline node & depth & lsize & rsize & psrsize & lsrsize & rsrsize  " >> starttab
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
  echo -n "" > $temp
  for n in "${!nodes[@]}"; do
    mycommand="./main pernode ${dataset[0]} ${nodes[$n]} $srlookupsize"
    echo $mycommand 
    $mycommand >> $temp
  done
  for i in "${!runs[@]}"; do
    cp ${libs[$i]}/* $path_lib/current/
    make clean
    make
    temp="temp"$i
    temps=$temps" "$temp
    echo -n "" > $temp
    for n in "${!nodes[@]}"; do
      mycommand="./main pernode ${dataset[0]} ${nodes[$n]} ${runs[$i]} ${arch[$a]}"
      echo $mycommand 
      $mycommand >> $temp
    done
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




