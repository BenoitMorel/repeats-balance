
filename=$1
repeatsfilename=$(echo $filename | sed 's/tipinner/repeats/')
tipinnerfilename=$(echo $filename | sed 's/repeats/tipinner/')





threads=$(cat $filename | awk '{if (NR!=1) {print}}' | grep "parallelization: MPI (" | awk '{print $3}' | sed 's/(//')
sites=$(cat $filename | awk '{if (NR!=1) {print}}' | grep "Alignment comprises " | awk '{print $6}')
elapsed_ti=$(cat $tipinnerfilename | awk '{if (NR!=1) {print}}' | grep "Elapsed time:" | awk '{print $3}' | sed 's/\..*//')
elapsed_rep=$(cat $repeatsfilename | awk '{if (NR!=1) {print}}' | grep "Elapsed time:" | awk '{print $3}' | sed 's/\..*//')
losttime_ti=$(cat $tipinnerfilename | awk '{if (NR!=1) {print}}' | grep "average possible speedup " | awk 'END{print $4}')
losttime_rep=$(cat $repeatsfilename | awk '{if (NR!=1) {print}}' | grep "average possible speedup " | awk 'END{print $4}')

losttime_ti=${losttime_ti:0:4}\\\%
losttime_rep=${losttime_rep:0:4}\\\%
sites_per_thread=$(($sites / $threads))
elapsed_ti_per_thread=$(($elapsed_ti * $threads))
elapsed_rep_per_thread=$(($elapsed_rep * $threads))
speedup=$(($elapsed_ti_per_thread * 100 / $elapsed_rep_per_thread))
speedup=$(echo $speedup | sed 's/\([0-9]\)/\1./')

echo "\\hline" $threads "&" $sites_per_thread "&" ${elapsed_ti_per_thread}s "&" ${elapsed_rep_per_thread}s "&" $losttime_ti "&" $losttime_rep "&" $speedup "\\\\"


