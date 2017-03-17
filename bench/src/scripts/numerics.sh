
iterations=200
./main numerics ../../data/404/unrooted.newick ../../data/404/404.phy 4 1 $iterations 42 &> repeats.txt
#./main numerics ../../data/404/unrooted.newick ../../data/404/404.phy 4 2 $iterations 42 &> noopt.txt
./main numerics ../../data/404/unrooted.newick ../../data/404/404.phy 4 0 $iterations 42 &> tipinner.txt


