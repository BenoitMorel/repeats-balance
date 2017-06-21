#include "samples.hpp"
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "../common/repeatsbalance.hpp"

void print_random_trees(int argc, char *params[])
{
  if (argc != 3) {
    std::cerr << "Error : syntax is" << std::endl;
    std::cerr 
      << "phylip states tree_numbers" 
      << std::endl;
    return ;
  }
  unsigned int i = 0;
  const char *phy_filename = params[i++];
  unsigned int states_number = atoi(params[i++]);
  unsigned int trees_number = atoi(params[i++]);


  srand(time(0));
  MSA msa(phy_filename, states_number);
  Tree tree(&msa);

  tree.update_operations(Tree::traverser_full);
  
  for (unsigned int i = 0; i < trees_number; ++i) {
    tree.randomize_pll_utree(&msa);
    tree.print();
  }

}

