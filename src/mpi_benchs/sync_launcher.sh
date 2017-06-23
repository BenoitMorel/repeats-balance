#!/bin/bash


submit_file=

write()
{
    echo $1 >> $submit_file
}

launch()
{
  data=$1
  threads=$2
  nodes=$(echo "($threads - 1)/16 + 1" | bc)
  randomized=$3
  iterations=$4
  sample_trees=$5
  trees=$6
  use_barrier=$7

  part_suffix=
  results_file=results
  update_repeats=1
  use_repeats=1
  seed=45
  states=4
  lookupsize=0
  path_data=../../data/
  submit_file="submit_${data}_${threads}_${sample_trees}_${trees}_${iterations}"
  rm -f  $submit_file
  write  "#!/bin/bash"
  write "#SBATCH -o ng_$submit_file_%j.out"
  write "#SBATCH -N $nodes"
  write "#SBATCH -n $threads"
  write "#SBATCH -B 2:8:1"
  write "#SBATCH --threads-per-core=1"
  write "#SBATCH --cpus-per-task=1"
  write "#SBATCH -t 1:00:00"

  write "source /etc/profile.d/modules.sh"
  write "module load gompi"
  write "export LD_LIBRARY_PATH=../../lib/current:../common"
  write "path_data=../../data/"

  mkdir -p $results_file
  mkdir -p $results_file/$data$part_suffix
  mkdir -p $results_file/$data$part_suffix/$threads
  outputfile=$results_file/$data$part_suffix/$threads/ng_${sample_trees}_${trees}_${use_barrier}_${randomized}.out

  re='^[0-9]+$'
  if ! [[ $trees = "random" ]] ; then
    trees=$path_data/$data/$trees
  fi
  if ! [[ $sample_trees =~ $re ]] ; then
    sample_trees=$path_data/$data/$sample_trees
  fi

  echo $trees
  echo $sample_trees
  preload=/home/morelbt/libjemalloc.so
  write "LD_PRELOAD=$preload mpirun -np $threads ./main synchronized_ftt $path_data/$data/$data.phy $path_data/$data/${data}${part_suffix}.part $states $use_repeats $update_repeats $lookupsize $iterations $randomized $seed $sample_trees $trees $use_barrier  &> $outputfile"

  sbatch $submit_file
  rm $submit_file
}


export LD_LIBRARY_PATH=../../lib/current:../common

#launch kyte 128 0 1 300
#launch kyte 128 1 1 300

data=kyte
#data=1kite_2013_10randomtaxa
threads=256
iterations=10000
#sample_tree=worsttree.newick
sample_tree=10
tree=random

#launch $data $threads 1 $use_repeats $iterations 20 0 1 0

echo $tree
echo $sample_tree

# always recreate tree, no barriers
#launch $data $threads 0 $iterations $sample_tree $tree 0 
#launch $data $threads 1 $iterations $sample_tree $tree 0 

# always recreate tree, barriers
launch $data $threads 0 $iterations $sample_tree random 1 
launch $data $threads 1 $iterations $sample_tree random 1 



