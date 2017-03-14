for file in ./*.out
do
  ./out_to_tex.sh "$file" "${file%.out}.tex"
done
