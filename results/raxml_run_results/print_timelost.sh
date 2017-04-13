


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
  losttime_ti=$(cat $tipinnerfilename | awk '{if (NR!=1) {print}}' | grep "average possible speedup " | awk 'END{print $4}')
  losttime_rep=$(cat $repeatsfilename | awk '{if (NR!=1) {print}}' | grep "average possible speedup " | awk 'END{print $4}')
  reducetime_ti=$(cat $tipinnerfilename | awk '{if (NR!=1) {print}}' | grep "average unavoidable wait" | awk 'END{print $4}')
  reducetime_rep=$(cat $repeatsfilename | awk '{if (NR!=1) {print}}' | grep "average unavoidable wait" | awk 'END{print $4}')
  reducetime_ti=$(echo $reducetime_ti | awk '$0="scale=2;100*"$0' | bc | awk '$0=$0"\\%"')
  reducetime_rep=$(echo $reducetime_rep | awk '$0="scale=2;100*"$0' | bc | awk '$0=$0"\\%"')

  losttime_ti=${losttime_ti:0:4}\\\%
  losttime_rep=${losttime_rep:0:4}\\\%
  reducetime_ti=${reducetime_ti:0:4}\\\%
  reducetime_rep=${reducetime_rep:0:4}\\\%
  sites_per_thread=$(($sites / $threads))
  elapsed_ti_per_thread=$(($elapsed_ti * $threads))
  elapsed_rep_per_thread=$(($elapsed_rep * $threads))
  speedup=$( echo $elapsed_ti_per_thread "/" $elapsed_rep_per_thread | awk '$0="scale=2;"$0' | bc  )

  echo "\\hline" $threads "&" $sites_per_thread "&" $losttime_ti "&" $losttime_rep "&" $reducetime_ti "&" $reducetime_rep "\\\\"

}

analyse_data()
{
  print_escaped $2
  print_escaped ""
  print_escaped "\\begin{tabular}{|c|c|c|c|c|c|}"
  print_escaped "\\hline threads  &  sites/thread & bad LB & bad LB repeats & reduce & reduce repeats \\\\"
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


#echo "0.005" |  awk '$0="scale=2;100*"$0' | bc | awk '$0=$0"%"'

