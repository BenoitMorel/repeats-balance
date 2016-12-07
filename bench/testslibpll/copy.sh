
dest=$1
mkdir -p $dest
rm $dest/*
cp ../../../libpll/src/.libs/*.so* $dest
cp ../../../libpll/src/pll.h $dest
cp $dest/* .
