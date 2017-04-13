


print_escaped()
{
  echo $1 | sed 's/_/\\_/g'
}

analyse_run() 
{

  filename=$1
  repeatsfilename=$(echo $filename | sed 's/tipinner/repeats/')
  tipinnerfilename=$(echo $filename | sed 's/repeats/tipinner/')

  repeatsfinished=$(cat $repeatsfilename | awk '{if (NR!=1) {print}}' | grep "Elapsed time:")
  tipinnerfinished=$(cat $tipinnerfilename | awk '{if (NR!=1) {print}}' | grep "Elapsed time:")
  if [[ -z "${repeatsfinished// }" ]] 
  then
    return
  fi
  if [[ -z "${tipinnerfinished// }" ]] 
  then
    return
  fi


  threads=$(cat $filename | awk '{if (NR!=1) {print}}' | sed 's|NONE/sequential|MPI (1 ranks)|g' | grep "parallelization: MPI (" | awk '{print $3}' | sed 's/(//')
  sites=$(cat $filename | awk '{if (NR!=1) {print}}' | grep "Alignment comprises " | awk '{print $6}')
  elapsed_ti=$(cat $tipinnerfilename | awk '{if (NR!=1) {print}}' | grep "Elapsed time:" | awk '{print $3}' | sed 's/\..*//')
  elapsed_rep=$(cat $repeatsfilename | awk '{if (NR!=1) {print}}' | grep "Elapsed time:" | awk '{print $3}' | sed 's/\..*//')

  sites_per_thread=$(($sites / $threads))
  elapsed_ti_per_thread=$(($elapsed_ti * $threads))
  elapsed_rep_per_thread=$(($elapsed_rep * $threads))
  speedup=$( echo $elapsed_ti_per_thread "/" $elapsed_rep_per_thread | awk '$0="scale=2;"$0' | bc  )

  echo "\\hline" $threads "&" $sites_per_thread "&" ${elapsed_ti}s "&" ${elapsed_rep}s "&" ${elapsed_ti_per_thread}s "&" ${elapsed_rep_per_thread}s "&" $speedup "\\\\"

}

analyse_data()
{
  print_escaped $2
  print_escaped ""
  print_escaped "\\begin{tabular}{|c|c|c|c|c|c|c|}"
  print_escaped "\\hline threads  &  sites/thread & time & time repeats & time*threads & time*threads repeats & speedup \\\\"
  for rundir in `ls -v $1/$2 | grep repeats `; do
    analyse_run  $1/$2/$rundir/$2.raxml.log
  done 
  echo "\\hline"
  echo "\\end{tabular}"
  echo "\\newline"
  echo ""
}



dir=${1%/}
for datadir in `ls $dir`; do
  analyse_data  $dir $datadir
done
