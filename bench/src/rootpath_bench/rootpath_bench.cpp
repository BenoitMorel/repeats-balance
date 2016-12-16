#include "../common/common.h"
#include <iostream>
#include <time.h>

int main(int argc, char *params[])
{
  srand(42);
  if (argc != 8) {
    std::cout << "Expected syntax : ./rootpath_bench newick phy "
      "iterations use_repeats update_repeats p paths_number" << std::endl;
    return 1;
  }
  const char *tree = params[1];
  const char *seq = params[2];
  unsigned int iterations = atoi(params[3]);
  bool use_repeats = atoi(params[4]);
  bool update_repeats = atoi(params[5]) && use_repeats;
  double p = atof(params[6]);
  unsigned int paths_number = atoi(params[7]);
  unsigned int attribute = 0;
  if (use_repeats) {
    std::cout << "sites repeats" << std::endl;
    attribute |= PLL_ATTRIB_SITES_REPEATS;
  } else {
    attribute |= PLL_ATTRIB_PATTERN_TIP;	
  }
  PLLHelper d(tree, seq, attribute);
  d.update_all_partials();

  std::vector<pll_operation_t> ops;

  for (unsigned int i = 0; i < paths_number; ++i) {
    do {
      d.build_random_path(p, ops);
    } while (ops.size() < 11);
    Timer t;
    d.update_partials(ops, iterations, update_repeats);
    unsigned int elapsed = t.get_time();
    std::cout  <<  ops.size() << "," << elapsed << std::endl;
  }
  return 1;
}
