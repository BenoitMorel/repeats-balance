#!/bin/bash

file=$1


max=$(cat $file | grep "ms" | awk 'END{print $2}' | sed 's/ms//g' | xargs)
plop=$(cat $file | grep "ms" | awk '{print "(" NR "," $2 ")"} END{print "(" NR+1 ", 0)"}' | sed 's/ms//g' | xargs)
caption=$(cat $file | grep "caption" | sed 's/caption//g')
echo "\\begin{minipage}{0.49\\textwidth}"
echo "  \\begin{tikzpicture}[scale=0.75]"
echo "    \\begin{axis}[ybar interval, ymax=$max,ymin=0, minor y tick num = 3, xlabel={CPU}, ylabel={Elapsed time (ms)}]"
echo "      \\addplot coordinates {$plop};"
echo "     \\end{axis}"
echo "  \\end{tikzpicture}"
echo "  \\caption*{$caption}"
echo "\\end{minipage}"



