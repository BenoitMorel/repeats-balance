#!/bin/bash


submit_file=

write()
{
    echo $1 >> $submit_file
}

launch()
{
  threads=$2
  data=$1
  nodes=$(echo "($threads - 1)/16 + 1" | bc)
  randomized=$3
  update_repeats=$4
  iterations=$5
  use_repeats=1
  seed=42
  states=4
  lookupsize=0
  path_data=../../../data/

  submit_file="submit_${data}_${threads}_${randomized}_${update_repeats}_${iterations}"
  rm -f  $submit_file
  write  "#!/bin/bash"
  write "#SBATCH -o ng_$submit_file_%j.out"
  write "#SBATCH -N $nodes"
  write "#SBATCH -n $threads"
  write "#SBATCH -B 2:8:1"
  write "#SBATCH --threads-per-core=1"
  write "#SBATCH --cpus-per-task=1"
  write "#SBATCH -t 24:00:00"

  write "source /etc/profile.d/modules.sh"
  write "module load gompi"
  write "export LD_LIBRARY_PATH=../../lib/current:../common"
  write "path_data=../../../data/"


  preload=/home/morelbt/libjemalloc.so
  write "LD_PRELOAD=$preload mpirun -np $threads ./main full_tree_traversal $path_data/$data/unrooted.newick $path_data/$data/$data.phy $path_data/$data/$data.part $states $use_repeats $update_repeats $lookupsize $iterations $randomized $seed"


  sbatch $submit_file
  rm $submit_file

}


export LD_LIBRARY_PATH=../../lib/current:../common

launch kyte 128 0 1 300
launch kyte 128 1 1 300

#launch kyte 128 0 1 300
#launch kyte 128 1 1 300

#launch 404 16 0 1 3000
#launch 404 16 1 1 3000



