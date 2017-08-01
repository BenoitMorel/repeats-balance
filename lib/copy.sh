
from=$1
to=$2
echo "usage : ./copy.sh libpll_rep_path dest_path"
mkdir -p $to
touch $to/plop
cp $from/src/.libs/*.so* $to
cp $from/src/pll.h $to
cp $to/* current
echo "done"
