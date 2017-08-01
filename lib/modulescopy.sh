if [ "$#" -ne 2 ]; then
  echo "usage : ./copy.sh pllmodules_rep_path dest_path"
  exit
fi

modulespath=$1
libpllpath=$modulespath/libs/libpll
dest=$2


mkdir -p $dest
cp $libpllpath/src/.libs/*.so* $dest
cp $libpllpath/src/pll.h $dest
cp $modulespath/src/*/.libs/*so* $dest
cp $modulespath/src/*/*.h $dest
cp $dest/* current
echo "done"
