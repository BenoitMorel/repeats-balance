#include "../common/common.h"
#include <iostream>

int main(int argc, char *params[])
{
  if (argc != 6) {
    std::cout << "Expected syntax : ./update_per_level newick phy "
      "iterations use_repeats update_repeats" << std::endl;
    return 1;
  }
  const char *tree = params[1];
  const char *seq = params[2];
  unsigned int iterations = atoi(params[3]);
  bool use_repeats = atoi(params[4]);
  bool update_repeats = atoi(params[5]) && use_repeats;
  unsigned int attribute = 0;
  if (use_repeats) {
    attribute |= PLL_ATTRIB_SITES_REPEATS;
  } else {
    attribute |= PLL_ATTRIB_PATTERN_TIP;	
  }
  PLLHelper d(tree, seq, attribute);
  d.update_all_partials();

  std::vector<unsigned int> children_number(d.nodes_count * 2);
  d.fill_children_number(children_number);

  for (unsigned int i = 0; i < d.inner_nodes_count; ++i) {
    Timer t;
    d.update_partial(&d.operations[i], iterations, update_repeats);
    unsigned int elapsed = t.get_time();
    std::cout  <<  children_number[d.operations[i].parent_clv_index] << "," << elapsed << std::endl;
  }
  return 1;
}
