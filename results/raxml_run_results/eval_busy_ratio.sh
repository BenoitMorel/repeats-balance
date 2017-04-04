cat $1/*.log | grep busy | awk  '{if (NR == 1) { print "("} if (NR != 1) {printf " + ";}  print $2; } END{print ") * 100 /"  NR ".0";}' | xargs | bc
