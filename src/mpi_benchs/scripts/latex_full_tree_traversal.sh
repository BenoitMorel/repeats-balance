#!/bin/bash

file=$1


max=$(cat $file | grep "max_time" | awk 'END{print $2}' | sed 's/ms//g' | xargs)
threads=$(cat $file | grep "ms" | awk 'END{print NR}') 
maxy=$(((10 * $max) / 9))
list_times=$(cat $file | grep "ms" | awk '{print "$2"} ' | sed 's/ms//g' )
plop=$(cat $file | grep "ms" | awk '{print "(" NR "," $2 ")"} END{print "(" NR+1 ", 0)"}' | sed 's/ms//g' | xargs)
caption=$(cat $file | grep "caption" | sed 's/caption//g')
echo "\\begin{minipage}{0.24\\textwidth}"
echo "  \\begin{tikzpicture}[scale=0.50]"
echo "    \\begin{axis}[ybar interval, xtick=\\empty, ymax=$maxy,ymin=0, minor y tick num = 3, xlabel={Cores}, ylabel={Elapsed time (ms)}]"
echo "      \\addplot coordinates {$plop};"
echo "        \\draw [red, dashed] ({rel axis cs:0,0}|-{axis cs:0,$max}) -- ({rel axis cs:1,0}|-{axis cs:$threads,$max}) node [pos=0.5, below] {ODDA worst cpu time};
"
echo "     \\end{axis}"
echo "  \\end{tikzpicture}"
echo "  %\\caption*{$caption}"
echo "\\end{minipage}"
echo ""

average_time=$(cat $file | grep "ms" | awk '{print $2} ' | sed 's/ms//g' | awk '{s+=$1}END{print s/NR}') 

echo $max  $average_time

echo "$max $average_time" | awk '{printf "%.2f\n", 1 - $2 / $1}'
