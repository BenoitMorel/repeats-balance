cd ../../../libpll
git checkout dev
make clean
make
cp src/pll.h src/.libs/libpll.so* ../repeats-balance/bench/lib/lipll_benoit_dev
git checkout tipinner
make clean
make
cp src/pll.h src/.libs/libpll.so* ../repeats-balance/bench/lib/lipll_benoit_tipinner
cd ../repeats-balance/bench/lib/

