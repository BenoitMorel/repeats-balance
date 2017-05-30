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
  update_repeats=$4
  iterations=$5
  trees_number=$6
  use_randomized_tree=$7
  use_barrier=$8
  use_update_operations=$9

  part_suffix=

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
  write "#SBATCH -t 1:00:00"

  write "source /etc/profile.d/modules.sh"
  write "module load gompi"
  write "export LD_LIBRARY_PATH=../../lib/current:../common"
  write "path_data=../../../data/"

  mkdir -p results
  mkdir -p results/$data$part_suffix
  mkdir -p results/$data$part_suffix/$threads
  outputfile=results/$data$part_suffix/$threads/ng_${randomized}_${trees_number}_${use_randomized_tree}_${use_barrier}.out

  preload=/home/morelbt/libjemalloc.so
  write "LD_PRELOAD=$preload mpirun -np $threads ./main synchronized_ftt $path_data/$data/$data.phy $path_data/$data/${data}${part_suffix}.part $states $use_repeats $update_repeats $lookupsize $iterations $randomized $seed $trees_number $use_randomized_tree $use_barrier $use_update_operations &> $outputfile"

  sbatch $submit_file
  rm $submit_file
}


export LD_LIBRARY_PATH=../../lib/current:../common

#launch kyte 128 0 1 300
#launch kyte 128 1 1 300

data=generate_50_1000
threads=1024
iterations=10000
use_repeats=1

#launch $data $threads 1 $use_repeats $iterations 20 0 1 0

# always recreate tree, no barriers
launch $data $threads 0 $use_repeats $iterations 1 1 0 1
launch $data $threads 1 $use_repeats $iterations 1 1 0 1



# always recreate tree,  barriers
#launch $data $threads 0 $use_repeats $iterations 1 1 1 1
#launch $data $threads 1 $use_repeats $iterations 1 1 1 1




