#include "samples.hpp"
#include <iostream>
#include <stdlib.h>
#include "../common/repeatsbalance.hpp"
#include <time.h>

void random_trees_likelihoods(int argc, char *params[])
{
  if (argc != 4) {
    std::cerr << "Error : syntax is" << std::endl;
    std::cerr 
      << "phylip part_file states trees_number" 
      << std::endl;
    return ;
  }
  unsigned int i = 0;
  const char *phy_filename = params[i++];
  const char *part_filename = params[i++];
  unsigned int states_number = atoi(params[i++]);
  unsigned int trees_number = atoi(params[i++]);

  unsigned int attribute = Partition::compute_attribute(true, 
      0, 
		  "avx");

  srand(time(0));
  MSA msa(phy_filename, states_number);
 
  Tree tree(&msa);
  LikelihoodEngine engine(&tree, &msa, part_filename, attribute, states_number, 4, 0);  
  
  for (unsigned int i = 0; i < trees_number; ++i) {
    tree.randomize_pll_utree(&msa);
    engine.update_operations();
    engine.update_matrices();
    engine.update_partials();
    std::cout << engine.compute_likelihood() << std::endl;
    //engine.get_partitions()[0]->print_partials();
  }
}

