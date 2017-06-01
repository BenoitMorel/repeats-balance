#include <iostream>
#include "common.h"

/*
 *  For each node, print : left and right children number, and clv_computation time
 */
void pernode(int argc, char *params[])
{
  if (argc != 10 && argc != 9 && argc != 5 && argc != 4) {
    std::cerr << "Error : syntax is" << std::endl;
    std::cerr 
      << "newick sequence states use_repeats update_repeats additional_attr repeats_lookup_size iterations arch" 
      << std::endl;
    std::cerr 
      << "newick sequence states node_index use_repeats update_repeats additional_attr repeats_lookup_size iterations arch" 
      << std::endl;
    std::cerr 
      << "newick sequence states repeats_lookup_size" 
      << std::endl;
    std::cerr 
      << "newick sequence states node_index repeats_lookup_size" 
      << std::endl;
    return ;
  }
  unsigned int i = 0;
  const char *newick = params[i++];
  const char *seq = params[i++];
  unsigned int states = atoi(params[i++]);
  int node_index = -1;
  if (argc == 10 || argc == 5) {
    node_index = atoi(params[i++]);
  }
  // just print nodes stats
  if (argc == 4 || argc == 5) {
    PLLHelper helper(newick, seq, PLL_ATTRIB_SITES_REPEATS | PLL_ATTRIB_ARCH_AVX); 
    unsigned int repeats_lookup_size = atoi(params[i++]);
    helper.set_srlookup_size(repeats_lookup_size);
    helper.update_all_partials();
    helper.init_tree_stats();
    if (argc == 4) {
      for (unsigned int op = 0; op < helper.ops_count; ++op) {
        std::cout << op << " ";
        helper.print_op_stats(helper.operations[op]);
      }
    } else {
      std::cout << node_index << " ";
      helper.print_op_stats(helper.operations[node_index]);

    }
    return;
  }

  unsigned int use_repeats = atoi(params[i++]);
  unsigned int update_repeats = atoi(params[i++]);
  unsigned int additional_attr = atoi(params[i++]);
  unsigned int repeats_lookup_size = atoi(params[i++]);
  unsigned int iterations = atoi(params[i++]);
  const char *arch = params[i++];
  

  unsigned int attribute = PLLHelper::compute_attribute(use_repeats, additional_attr, arch);
  if (INVALID_ATTRIBUTE == attribute) {
    return;
  }
  PLLHelper helper(newick, seq, attribute, states); 
  helper.set_srlookup_size(repeats_lookup_size);
  helper.update_all_partials();
 
  if (node_index == -1) {
    for (unsigned int op = 0; op < helper.ops_count; ++op) {
      Timer t;
      helper.update_partial(helper.operations[op], iterations, update_repeats);
      std::cout << t.get_time() << "ms" << std::endl; 
    }
  } else {
    Timer t;
    helper.update_partial(helper.operations[node_index], iterations, update_repeats);
    std::cout << t.get_time() << "ms" << std::endl; 

  }
}



