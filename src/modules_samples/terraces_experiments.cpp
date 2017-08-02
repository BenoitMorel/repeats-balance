


#include <iostream>



extern "C" {
#include "pll_tree.h"
#include "pllmod_algorithm.h"
}
#include "../common/repeatsbalance.hpp"


pllmod_treeinfo_t *load_treeinfo(const char *phy_file,
        const char *part_file)
{
  const unsigned int states_number = 4;
  const unsigned int rates = 4;
  unsigned int attributes = PLL_ATTRIB_SITE_REPEATS | PLL_ATTRIB_ARCH_AVX;
  MSA full_msa(phy_file, states_number); 
  Tree tree(&full_msa);
  std::vector<PartitionIntervals> partitionning;
  PartitionIntervals::parse(part_file, partitionning);
  std::vector<Partition *> partitions; 
  for (unsigned int i = 0; i < partitionning.size(); ++i) 
  {
    MSA submsa(&full_msa, partitionning[i], i);
    submsa.compress();
    partitions.push_back(new Partition(&submsa,
        attributes,
        states_number,
        rates,
        0));
  }
  
  // pllmod stuff
  pllmod_treeinfo_t * treeinfo = pllmod_treeinfo_create(tree.get_pll_root(),
                                                        tree.get_tips_number(),
                                                        partitions.size(),
                                                        PLLMOD_TREE_BRLEN_LINKED);
  for (unsigned int i = 0; i < partitionning.size(); ++i) 
  {
    unsigned int params_indices[4] = {0,0,0,0};
    const unsigned int params_to_optimize = 0;
    pllmod_treeinfo_init_partition(treeinfo,
                                  i,
                                  partitions[i]->get_partition(),
                                  params_to_optimize,
                                  0, 
                                  1.0, 
                                  params_indices,
                                  NULL  
                                  );
  }
  return treeinfo;
    
}

void terraces_experiments(int argc, char *argv[])
{
  srand(42);
  pllmod_treeinfo_t *treeinfo = load_treeinfo("../../data/404/404.phy", "../../data/404/404.part");
  std::cout << "likelihood: " << pllmod_treeinfo_compute_loglh(treeinfo, 0) << std::endl;


}

