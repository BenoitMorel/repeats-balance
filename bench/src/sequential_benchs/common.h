#ifndef _LIBPLLBENCH_COMMON_H
#define _LIBPLLBENCH_COMMON_H

#include "../common/Tree.hpp"
#include "../common/Partition.hpp"
#include "../common/Model.hpp"
#include "../common/LikelihoodEngine.hpp"
#include <time.h>
#include <vector>
#include "depreciated/PLLHelper.hpp"

#ifndef PLL_ATTRIB_SITES_REPEATS // old pll version
  #define PLL_ATTRIB_SITES_REPEATS 0
#else
  #define HAS_REPEATS
#endif

#define INVALID_ATTRIBUTE ((unsigned int)1)

class Timer {
public:
  Timer();
  // ms
  long get_time();

private:
  timespec start;
};


void full_traversal(int argc, char *params[]);
void derivatives(int argc, char *params[]);
void sumtables(int argc, char *params[]);
void partitioned_full_traversal(int argc, char *params[]);
void pernode(int argc, char *params[]);
void repeats_rates(int argc, char *params[]);
void export_repeats(int argc, char *params[]);

#endif
