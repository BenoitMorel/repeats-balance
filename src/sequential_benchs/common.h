#ifndef _LIBPLLBENCH_COMMON_H
#define _LIBPLLBENCH_COMMON_H

#include "../common/repeatsbalance.hpp"
#include <vector>
#include "depreciated/PLLHelper.hpp"

#ifndef PLL_ATTRIB_SITES_REPEATS // old pll version
  #define PLL_ATTRIB_SITES_REPEATS 0
#else
  #define HAS_REPEATS
#endif

#define INVALID_ATTRIBUTE ((unsigned int)1)


void full_traversal(int argc, char *params[]);
void core_functions(int argc, char *params[]);
void sumtables(int argc, char *params[]);
void partitioned_full_traversal(int argc, char *params[]);
void pernode(int argc, char *params[]);
void repeats_rates(int argc, char *params[]);
void export_repeats(int argc, char *params[]);
void kassian_lb_partials(int argc, char *params[]);

#endif
