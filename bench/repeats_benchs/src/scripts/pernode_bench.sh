#!/bin/bash
outputdir="../../../results/pernode_benchs/"
tempfile1="temp1"
tempfile2="temp2"
iterations=1000

export LD_LIBRARY_PATH=../lib/current/

write() {
  echo $1 >> $tempfile1
}

launch() {
  echo "./pernode_bench/pernode_bench $1 $2 $3 $4 $5 "
  ./pernode_bench/pernode_bench $1 $2 $3 $4 $5 > $tempfile2

  write '\begin{minipage}{0.49\textwidth}'
  write '\caption*{'
  if [ $4 = "0" ]
  then
   write 'tip pattern mode' 
  else
    write 'sites repeats mode'
  fi
  write '}'
  write '\begin{tikzpicture}[scale=0.75]'
  write '  \begin{axis}[ymin=0, xlabel={Recursive number of children}, ylabel={update partial time}] \addplot[mark=x, color=blue, scatter, only marks] coordinates {'

  awk -F ',' '{print "("$1","$2") "}' $tempfile2 | tr -d '\n' >> $tempfile1

  write '};'
  write '  \end{axis}' 
  write '\end{tikzpicture}'
  write '\end{minipage}'

}

bench() {
  outputfile=$1
  outputfile="$outputdir$outputfile"
  echo "" > $tempfile1
  launch $2 $3 $iterations 0 0
  launch $2 $3 $iterations 1 0
  mv $tempfile1 $outputfile
}

make
mkdir -p $outputdir
bench results59.tex ../../../data/59/unrooted.newick  ../../../data/59/59.phy
bench results128.tex ../../../data/128/unrooted.newick  ../../../data/128/128.phy
bench results404.tex ../../../data/404/unrooted.newick  ../../../data/404/404.phy


