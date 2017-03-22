cat initial_assign.txt | awk '{print "  M[" $2 "][" $3 "] = "  $5 ";";}'
