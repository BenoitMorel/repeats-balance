#!/bin/bash
source ~/github/benoit/latex_scripts/histograms.sh

inputfile=$1
file=$2
threads=4

print_wait_thread() {
  inputfile=$1
  threadindex=$2
  cat $inputfile | grep "thread $threadindex" | awk '{print "(" NR "," $4 ")"}' | xargs echo -n
}

print_cpu_thread() {
  inputfile=$1
  threadindex=$2
  cat $inputfile | grep "thread $threadindex" | awk '{print "(" NR "," $7 ")"}' | xargs echo -n
}




echo "" > $file

tempfile="tempplop.txt"
tempfilecaption="tempcaption.txt"
tempfileratio="tempratio.txt"
cat $inputfile | grep cpu_time_ms  | sed 's/cpu_time_ms//g' | sed 's/thread//g' > $tempfile
cat $inputfile | grep stats_step | awk '{print $2}'  | sed 's/_/\\_/g' > $tempfilecaption
cat $inputfile | grep busy_ratio | awk '{print $2}'   > $tempfileratio

n=0
while read -r p && read -r q <&3 && read -r s <&4; do
  if [ $n == 0 ] 
  then 
    write_figure_header $file 
  fi
  write_histo_header $file "thread id" "elapsed CPU time (ms)" "$q ratio:$s"
  write_coordinates_header $file
  echo $p >> $file
  write_coordinates_footer $file
  write_histo_footer $file
  if [ $n == 2 ] 
  then 
    write_figure_footer $file
    n=0
  elif [ $n == 1 ]
  then
    n=2
  else
    n=1
  fi
done <  $tempfile 3< $tempfilecaption 4< $tempfileratio

if [ $n != 0 ]
then
  write_figure_footer $file
fi

rm $tempfile
rm $tempfilecaption
rm $tempfileratio







#write_histo_header $file "core" "elapsed CPU time"
#for i in $(seq 0 $(($threads - 1))); 
#do
#  echo "plop"
#  write_coordinates_header $file
#  print_cpu_thread $inputfile $i >> $file
#  write_coordinates_footer $file
#done
#write_histo_footer $file

#echo "" >> $file

#write_histo_header $file "core" "elapsed CPU time"
#for i in $(seq 0 $(($threads - 1))); 
#do 
#  write_coordinates_header $file
#  print_wait_thread $inputfile $i >> $file
#  write_coordinates_footer $file
#done
#write_histo_footer $file
##
