#include "tools.hpp"
#include <iostream>
#include <stdlib.h>
#include "../common/repeatsbalance.hpp"
#include <time.h>

void analyse_partitions(int argc, char *params[])
{
  if (argc != 3) {
    std::cerr << "Error : syntax is" << std::endl;
    std::cerr 
      << "phylip partition states" 
      << std::endl;
    return ;
  }
  unsigned int i = 0;
  const char *phy_filename = params[i++];
  const char *part_filename = params[i++];
  unsigned int states_number = atoi(params[i++]);

  unsigned int attribute = Partition::compute_attribute(true, 
      0, 
		  "avx");

  srand(time(0));
  MSA msa(phy_filename, states_number);
  std::vector<PartitionIntervals> partitions_intervals;
  PartitionIntervals::parse(part_filename, partitions_intervals);
  Tree tree(&msa);
  tree.update_operations(Tree::traverser_full);
  for (unsigned int i = 0; i < partitions_intervals.size(); ++i) {
    MSA submsa(&msa, partitions_intervals[i], i);
    unsigned int initial_sites = submsa.get_pll_msa()->length;
    submsa.compress();
    Partition partition(&submsa, attribute, states_number, 4, 0);
    partition.update_matrices(tree);
    partition.update_partials(tree);
    std::cout << "Partition " << i << ":" << std::endl;
    std::cout << "  pattern compression " << double(partition.get_sites_number())/double(initial_sites) << std::endl;
    std::cout << "  compressed sites: " << partition.get_sites_number() << std::endl;
    std::cout << "  unique pattern ratio: " << partition.get_unique_repeats_pattern_ratio() << std::endl;
  }
}

